use rustc_hash::FxHashMap;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use std::fs::{File, read_to_string};

use clap::Parser;
use anyhow::{bail, Result};

use ff_structure::DotBracketVec;
use ff_structure::PairList;
use ff_energy::NucleotideVec;

use fuzzyfold::category::Category;

#[derive(Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: Vec<PathBuf>,

    #[arg(long, value_name = "START-END", num_args = 1..)]
    regions: Vec<String>,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,
}

#[derive(Default)]
struct FileData {
    sequence: Option<NucleotideVec>,
    times: Vec<f64>,
    t_ext: f64,
    t_end: f64,
    t_sep: f64,
    num_sims: usize,
    structures: Vec<FxHashMap<DotBracketVec, usize>>,
    categories: Vec<FxHashMap<Category, usize>>,
}

fn main() -> Result<()> {

    let cli = Cli::parse();

    let regions: Vec<(usize, usize)> = cli.regions.iter()
    .map(|r| {
        let (a, b) = r.split_once('-').expect("region must be in format START-END");
        let start = a.parse::<usize>().expect("invalid start");
        let end = b.parse::<usize>().expect("invalid end");
        (start - 1, end - 1)
    })
    .collect();

    let mut data = Vec::new();

    for i in cli.input.iter() {
        let file = parse_file(i.clone())?;
        data.push(file);
    }

    if data.len() <= 1 {
        bail!("At least two files have to be given as input!")
    }
    
    for file in data.iter_mut() {
        file.categories = count_categories(&file.structures, &regions)?;
    }

    let csv_path = cli.output.with_extension("csv");
    
    let mut data_iter = data.into_iter();
    
    let mut combined_data = data_iter.next().unwrap();

    for file in data_iter {
        combined_data = accumulate(combined_data, file)?;
    }

    let sequence = combined_data.sequence.unwrap();
    let co_time = combined_data.t_ext * (sequence.len() - 1) as f64;
    let post_time = co_time + combined_data.t_end;
   
    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
    writeln!(writer, "# t_split={}", co_time)?;
    writeln!(writer, "# t_ext={:?}", combined_data.t_ext)?;
    writeln!(writer, "# t_end={}", post_time)?;
    let regions_str = regions.iter()
        .map(|(a, b)| format!("{}-{}", a, b))
        .collect::<Vec<_>>()
        .join(",");
    writeln!(writer, "# regions={}", regions_str)?;
    writeln!(writer, "# num_sims={}", combined_data.num_sims)?;
    write_categories(&mut writer, &combined_data.times, &combined_data.categories)?;

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


fn count_categories(structures: &Vec<FxHashMap<DotBracketVec, usize>>, regions: &Vec<(usize, usize)>,) -> Result<Vec<FxHashMap<Category, usize>>>{

    let mut result = Vec::with_capacity(structures.len());
   
    for timepoint_structures in structures.iter() {

        let mut counts: FxHashMap<Category, usize> = FxHashMap::default();

        for (structure, count) in timepoint_structures.iter() {

            let pl = PairList::try_from(structure)?;

            for (i, j) in pl.iter() {
                let region_i = regions.iter().position(|&(start, end) | (*i as usize) >= start && (*i as usize) <= end);
                let region_j = regions.iter().position(|&(start, end) | (*j as usize) >= start && (*j as usize) <= end);

                let category = match (region_i, region_j) {
                    (Some(a), Some(b)) if a == b => Category::Within(a),
                    (Some(a), Some(b)) => Category::Between(a.min(b), a.max(b)),
                    (Some(a), None) => Category::WithRest(a),
                    (None, Some(b)) => Category::WithRest(b),
                    (None, None) => continue,
                };

                *counts.entry(category).or_insert(0) += count;
            }
        }

        result.push(counts); 

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
    combined_data.categories = data1.categories.clone();

    for (t_idx, categories) in data2.categories.into_iter().enumerate() {
        for (cat, count) in categories.into_iter() {
            *combined_data.categories[t_idx].entry(cat).or_insert(0) += count;
        }
    }

    combined_data.num_sims = data1.num_sims + data2.num_sims;

    Ok(combined_data)
}
   

fn write_categories( 
    writer: &mut impl Write,
    times: &[f64],
    counts: &[FxHashMap<Category, usize>],
) -> Result<()> {
    writeln!(writer, "time,category,count")?;
    for (t_idx, t) in times.iter().enumerate() {
        for (cat, count) in counts[t_idx].iter() {
            writeln!(writer, "{},{},{}", t, cat.to_key(), count)?;
        }
    }
    Ok(())
}