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
    #[arg(value_name = "INPUT", default_value = "-")]
    input: Vec<PathBuf>,

    #[arg(short, long, value_name = "POSITION", num_args = 1..)]
    positions: Vec<usize>,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,
}

#[derive(Default)]
struct FileData {
    sequence: Option<NucleotideVec>,
    times: Vec<f64>,
    t_ext: usize,
    t_end: usize,
    t_sep: usize,
    num_sims: usize,
    structures: Vec<FxHashMap<DotBracketVec, usize>>,
    counts: Vec<Vec<usize>>,
}

fn main() -> Result<()> {

    let cli = Cli::parse();

    let zero_indexed_positions: Vec<usize> = cli.positions
        .iter()
        .map(|&position| {
            if position == 0 {
                Err(anyhow::anyhow!("Positions are 1-indexed, can't be zero!"))
            } else {
                Ok(position - 1)
            }
        })
        .collect::<Result<Vec<usize>,_>>()?;
    
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
        file.counts = edit_count(&sequence, &file.structures, &zero_indexed_positions);
    }

    let csv_path = cli.output.with_extension("csv");
    
    let mut data_iter = data.into_iter();
    
    let mut combined_data = data_iter.next().unwrap();

    for file in data_iter {
        combined_data = accumulate(combined_data, file)?;
    }
   
    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
    write_counts(&mut writer, &combined_data.times, &zero_indexed_positions, &combined_data.counts)?;

    Ok(())
}


fn parse_file(path: PathBuf) -> Result<FileData> {

    let content = read_to_string(path)?;
    let mut data = FileData::default();

    for line in content.lines() {
        if let Some(rest) = line.strip_prefix("# t_ext=") {
            data.t_ext = rest.parse()?; 
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

fn edit_count (sequence: &NucleotideVec, structures: &Vec<FxHashMap<DotBracketVec, usize>> , positions: &Vec<usize>) -> Vec<Vec<usize>> {

    structures.iter().map(|counts| {
        let mut edit_counts = vec![0usize; positions.len()];
        for (db, count) in counts.iter() { 
            let result = gets_edited(sequence, db, positions);
            for (i, edited) in result.iter().enumerate() {
                if *edited {
                    edit_counts[i] += count;
                }
            }     
        }
        edit_counts

    })
    .collect()
}

// checks whether position is unpaired and across from a 'C', and whether the 
// 5' and 3' neighbor are paired
fn gets_edited (sequence: &NucleotideVec, structure: &DotBracketVec, positions: &Vec<usize>) -> Vec<bool> {

    let mut result = Vec::new();
    let pt = PairTable::try_from(structure).unwrap();

    for position in positions.iter() {
       
        if structure[*position] != DotBracket::Unpaired {
            result.push(false);
            continue;
        }

        let outer = *position - 1;
        let inner = *position + 1;

        if !pt[outer].is_some() {
            result.push(false );
            continue; 
        } 
        
        if !pt[inner].is_some() {
            result.push(false);
            continue; 
        } 

        let outer_partner = pt[outer].unwrap() as usize;
        let inner_partner = pt[inner].unwrap() as usize;

        if outer_partner > inner_partner {
            if outer_partner - inner_partner == 2 && sequence[outer_partner - 1] == Base::C {
                result.push(true);
                continue;
            }

        } else if inner_partner - outer_partner == 2 && sequence[inner_partner - 1] == Base::C {
            result.push(true);
            continue;
        }
        result.push(false);
    }
        
    return result  
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
    combined_data.counts = data1.counts.clone();

   
    for (t_idx, pos_counts) in data2.counts.into_iter().enumerate() {
        for (i, count) in pos_counts.into_iter().enumerate() {
            combined_data.counts[t_idx][i] += count
        }
    }

    combined_data.num_sims = data1.num_sims + data2.num_sims;
    
    Ok(combined_data)
}


fn write_counts( 
    writer: &mut impl Write,
    times: & [f64],
    positions: &Vec<usize>,
    counts: &Vec<Vec<usize>>,
) -> Result<()> {
    let header = positions.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
    writeln!(writer, "# positions={}", header)?;
    write!(writer, "time,")?;
    writeln!(writer, "{}", header)?;
    for (t_idx, t) in times.iter().enumerate() {
        let row = positions.iter().enumerate()
            .map(|(i, _pos)| counts[t_idx][i].to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(writer, "{},{}", t, row)?;
    }
    Ok(())
}