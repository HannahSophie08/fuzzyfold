use rustc_hash::FxHashMap;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use std::fs::{File, read_to_string};

use clap::Parser;
use anyhow::{bail, Result};

use ff_energy::NucleotideVec;
use ff_energy::Base;

#[derive(Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like)
    #[arg(value_name = "INPUT")]
    input: Vec<PathBuf>,

    #[arg(long, default_value_t = 20)]
    num_sims: usize,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,
}

#[derive(Default)]
struct FileData {
    sequences: FxHashMap<NucleotideVec, f64>,
    num_sims: usize,
}


fn main() -> Result<()> {

    let cli = Cli::parse();


    let mut data = Vec::new();

    for i in cli.input.iter() {
        let file = parse_drf_file(i.clone(), cli.num_sims)?;
        data.push(file);
    }

    let csv_path = cli.output.with_extension("csv");
    
    let mut data_iter = data.into_iter();
    
    let mut combined_data = data_iter.next().unwrap();

    for file in data_iter {
        combined_data = accumulate(combined_data, file)?;
    }

    let rate_per_position = editing_rate(&combined_data)?;

    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
    writeln!(writer, "# num_sims={}", combined_data.num_sims)?;
    writeln!(writer, "position,rate")?;
    for (pos, rate) in rate_per_position.iter().enumerate() {
        writeln!(writer, "{}, {}", pos + 1, rate)?;
    }
    writer.flush()?;
    Ok(())
}


fn parse_drf_file(path: PathBuf, num_sims: usize) -> Result<FileData> {

    let content = read_to_string(path)?;
    let mut data = FileData::default();
    let mut last_t = 0.0;
    data.num_sims = num_sims;

    for line in content.lines() {
        if line.starts_with('#') {
            continue;
        }
        if line == "id time occupancy structure energy sequence" {
            continue;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        let t:f64 = cols[1].parse()?;

        last_t = t;
    }

    if last_t != 0.0 {
        for line in content.lines() {
            if line.starts_with('#') {
                continue;
            }
            if line == "id time occupancy structure energy sequence" {
                continue;
            }
            let cols: Vec<&str> = line.split_whitespace().collect();
            let t:f64 = cols[1].parse()?;
            if t == last_t {
                let sequence = NucleotideVec::try_from(cols[5])?;
                let count: f64 = cols[2].parse()?;
                *data.sequences.entry(sequence).or_insert(0.0) += count * num_sims as f64;
            }

        }
    }

    Ok(data)
}


fn accumulate(data1: FileData, data2: FileData) -> Result<FileData> {

    let mut combined_data = FileData::default(); 

    for (seq, count) in data1.sequences.into_iter() {
        *combined_data.sequences.entry(seq).or_insert(0.0) += count;
    }

    for (seq, count) in data2.sequences.into_iter() {
        *combined_data.sequences.entry(seq).or_insert(0.0) += count;
    }

    combined_data.num_sims = data1.num_sims + data2.num_sims;

    Ok(combined_data)
}
   

fn editing_rate(combined_data: &FileData) -> Result<Vec<f64>> {

    let num_sims = combined_data.num_sims;
    let Some(first) = combined_data.sequences.keys().next() else {
        bail!("no sequences found");
    };
    let len = first.len();
    let mut rate_per_position = vec![0.0; len];

    let mut total_count = 0.0;

    for (seq, count) in combined_data.sequences.iter() {
        for (idx, n) in seq.iter().enumerate() {
            if *n == Base::IN {
                rate_per_position[idx] +=*count; 
            }
        }
        total_count += count;
    }

    if (total_count - num_sims as f64).abs() > 1e-6 {
        eprintln!("Total count not the same as num_sims!");
    }

    for rate in rate_per_position.iter_mut() {
        *rate /= num_sims as f64;
    }

    Ok(rate_per_position)
}

