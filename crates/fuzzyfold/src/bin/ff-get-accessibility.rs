use rustc_hash::FxHashMap;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use std::fs::{File, read_to_string};

use clap::Parser;
use anyhow::{bail, Result};

use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::Base;

#[derive(Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT")]
    input: Vec<PathBuf>,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,
}

#[derive(Default)]
struct FileData {
    sequence: Option<NucleotideVec>,
    times: Vec<f64>,
    t_ext: f64,
    t_end: usize,
    t_sep: f64,
    num_sims: usize,
    structures: Vec<FxHashMap<DotBracketVec, usize>>,
    accessibility: Vec<Vec<usize>>,
}

fn main() -> Result<()> {

    let cli = Cli::parse();
    
    let mut data = Vec::new();

    for i in cli.input.iter() {
        let file = parse_file(i.clone())?;
        data.push(file);
    }
    
    if data.len() <= 1 {
        bail!("At least two files have to be given as input!")
    }

    for file in data.iter_mut() {
        let sequence = file.sequence.as_ref().ok_or_else(|| anyhow::anyhow!("missing sequence in input file"))?;
        file.accessibility = get_accessibility(&file.structures)?;
    }

    let csv_path = cli.output.with_extension("csv");
    
    let mut data_iter = data.into_iter();
    
    let mut combined_data = data_iter.next().unwrap();

    for file in data_iter {
        combined_data = accumulate(combined_data, file)?;
    }

    let normalized_accessibility: Vec<Vec<f64>> = combined_data.accessibility.iter()
        .map(|row| {
            row.iter().map(|&count| count as f64 / combined_data.num_sims as f64).collect()
        })
        .collect();
    let sequence = combined_data.sequence.unwrap();
   
    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
    writeln!(writer, "# sequence={}", sequence)?;
    writeln!(writer, "# t_ext={}", combined_data.t_ext)?;
    writeln!(writer, "# t_end={}", combined_data.t_end)?;
    writeln!(writer, "# t_sep={}", combined_data.t_sep)?;
    writeln!(writer, "# num_sims={}", combined_data.num_sims)?;
    write_accessibility(&mut writer, &combined_data.times, &normalized_accessibility)?;

    Ok(())
}


fn parse_file(path: PathBuf) -> Result<FileData> {

    let content = read_to_string(path)?;
    let mut data = FileData::default();

    for line in content.lines() {
        if let Some(rest) = line.strip_prefix("# t_ext=") {
            let inner = rest.strip_prefix("Some(").and_then(|s| s.strip_suffix(")")).unwrap_or(rest);
            data.t_ext = inner.parse()?; 
        }
        if let Some(rest) = line.strip_prefix("# t_end=") {
            data.t_end = rest.parse()?;
        }   
        if let Some(rest) = line.strip_prefix("# t_sep=") {
            data.t_sep = rest.parse()?;
        }
        if let Some(rest) = line.strip_prefix("# num_sims=") {
            data.num_sims = rest.parse()?; 
        }
        if let Some(rest) = line.strip_prefix("# sequence=") {
            data.sequence = Some(NucleotideVec::try_from(rest)?); 
        }    
    }

    for line in content.lines() {
        if line.starts_with('#') {
            continue;
        }
        if line == "time,structure,count" {
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        let t:f64 = cols[0].parse()?;

        if data.times.last() != Some(&t) {
            data.times.push(t);
            data.structures.push(FxHashMap::default());
        }

        let t_idx = data.times.len() - 1; 

        let structure = DotBracketVec::try_from(cols[1])?;
        let count: usize = cols[2].parse()?;

        *data.structures[t_idx].entry(structure).or_insert(0) += count;
    }
    Ok(data)
}

fn get_accessibility(structures: &Vec<FxHashMap<DotBracketVec, usize>>) -> Result<Vec<Vec<usize>>> {
    
    let mut result = Vec::new();
    for (i, counts) in structures.iter().enumerate() {
            let mut entries: Vec<(&DotBracketVec, &usize)> = counts.iter().collect();
            entries.sort_by(|a, b| b.1.cmp(a.1)); // descending by count
            let seq_len = entries.first().map(|(db, _)| db.0.len()).unwrap_or(0);
            for (db, _) in entries.iter() {
                if db.0.len() > seq_len + 1 {
                    bail!("timepoint {}: structure length {} exceeds expected {} by more than 1", i, db.0.len(), seq_len);
                }
            }
            let mut unpaired = vec![0usize; seq_len + 1];
            for (db, count) in entries {
                for (pos, state) in db.iter().enumerate() {
                    if matches!(state, DotBracket::Unpaired) {
                        unpaired[pos] += count;
                    }
                }
            }
            result.push(unpaired);
    }
    Ok(result)
}

fn accumulate(data1: FileData, data2: FileData) -> Result<FileData> {

    let mut combined_data = FileData::default();

    if data1.t_ext != data2.t_ext {
        bail!("t_ext isn't matching")
    } else {
        combined_data.t_ext = data1.t_ext;
    }

    if data1.t_end != data2.t_end {
        bail!("t_end isn't matching")
    } else {
        combined_data.t_end = data1.t_end;
    }

    if data1.t_sep != data2.t_sep {
        bail!("t_sep isn't matching")
    } else {
        combined_data.t_sep = data1.t_sep; 
    }

    if data1.sequence != data2.sequence {
        bail!("sequences aren't matching")
    } else {
        combined_data.sequence = data1.sequence;
    }

    const EPS: f64 = 1e-6; 

    if data1.times.len() != data2.times.len() {
        bail!("times aren't matching")
    }

    for (t1, t2) in data1.times.iter().zip(data2.times.iter()) {
        if (t1 - t2).abs() > EPS {
            bail!("times aren't matching")
        }
    }
    combined_data.times = data1.times.clone(); 
    combined_data.accessibility = data1.accessibility.clone();
   
    for (t_idx, pos) in data2.accessibility.into_iter().enumerate() {
        for (i, unpaired) in pos.into_iter().enumerate() {
            combined_data.accessibility[t_idx][i] += unpaired
        }
    }

    combined_data.num_sims = data1.num_sims + data2.num_sims;
    
    Ok(combined_data)
}


fn write_accessibility( 
    writer: &mut impl Write,
    times: & [f64],
    accessibility: &Vec<Vec<f64>>,
) -> Result<()> {
    write!(writer, "time,")?;
    for (t_idx, t) in times.iter().enumerate() {
        let row = accessibility[t_idx].iter()
            .map(|val| val.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(writer, "{},{}", t, row)?;
    }
    Ok(())
}