use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::path::PathBuf;

use clap::Parser;
use anyhow::bail;
use anyhow::Result;
use rustc_hash::FxHashMap;

use ff_alu_analysis::category::Category;

#[derive(Debug, Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    
    #[arg(value_name = "FILE")]
    input: Vec<PathBuf>,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,
}

#[derive(Default)]
struct FileData {
    t_ext: f64,
    t_end: f64,
    t_split: f64,
    regions: Vec<(usize, usize)>,
    num_sims: usize,
    times: Vec<f64>,
    counts: FxHashMap<(usize, Category), usize>,
}


fn main() -> Result<()> {
    let cli = Cli::parse();
    
    let csv_path = cli.output.with_extension("csv");

    let mut data = Vec::new();

    for i in cli.input.iter() {
        let file = parse_file(i.clone())?;
        data.push(file);
    }

    if data.len() <= 1 {
        bail!("At least two files have to be given as input!")
    }
    
    let mut data_iter = data.into_iter();
    
    let mut combined_data = data_iter.next().unwrap();

    for file in data_iter {

        combined_data = accumulate(combined_data, file)?;
    }
   
    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
        writeln!(writer, "# t_split={}", combined_data.t_split)?;
        writeln!(writer, "# t_ext={:?}", combined_data.t_ext)?;
        writeln!(writer, "# t_end={}", combined_data.t_end)?;
        let regions_str = combined_data.regions.iter()
            .map(|(a, b)| format!("{}-{}", a, b))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(writer, "# regions={}", regions_str)?;
        writeln!(writer, "# num_sims={}", combined_data.num_sims)?;
        write_categories(&mut writer, &combined_data.times, &combined_data.counts)?;

    Ok(())
}

fn parse_file(path: PathBuf) -> Result<FileData> {

    
    let content = fs::read_to_string(path)?;
    let mut data = FileData::default();

    for line in content.lines() {
        if let Some(rest) = line.strip_prefix("# t_ext=") {
            data.t_ext = rest.parse()?; 
        }
        if let Some(rest) = line.strip_prefix("# t_end=") {
            data.t_end = rest.parse()?;
        }   
        if let Some(rest) = line.strip_prefix("# t_split=") {
            data.t_split = rest.parse()?;
        }
        if let Some(rest) = line.strip_prefix("# regions=") {
            for piece in rest.split(',') {
                let (a, b) = piece.split_once('-').expect("region must be in format START-END");
                let a: usize = a.parse()?;
                let b: usize = b.parse()?;
                data.regions.push((a, b));
            }
        }
        if let Some(rest) = line.strip_prefix("# num_sims=") {
            data.num_sims = rest.parse()?; 
        }  
    }

    for line in content.lines() {
        if line.starts_with('#') {
            continue;
        }
        if line == "time,category,count" {
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        let t:f64 = cols[0].parse()?;

        if data.times.last() != Some(&t) {
            data.times.push(t);
        }

        let t_idx = data.times.len() - 1; 

        let category = Category::from_key(cols[1])?;
        let count: usize = cols[2].parse()?;

        data.counts.insert((t_idx, category), count);
    }
    Ok(data)
} 



fn accumulate(data1: FileData, data2: FileData) -> Result<FileData> {

    let mut combined_data = FileData::default();

    if data1.t_ext != data2.t_ext {
        bail!("t-ext isn't matching")
    } else {
        combined_data.t_ext = data1.t_ext
    }

    if data1.t_end != data2.t_end {
        bail!("t-end isn't matching")
    } else {
        combined_data.t_end = data1.t_end
    }

    if data1.t_split != data2.t_split {
        bail!("t-split isn't matching")
    } else {
        combined_data.t_split = data1.t_split
    }

    if data1.regions != data2.regions {
        bail!("regions aren't matching")
    } else {
        combined_data.regions = data1.regions 
    }

    combined_data.times = data1.times.clone();
    combined_data.counts = data1.counts.clone();

    const EPS: f64 = 1e-6;

    for (&(t_idx, ref cat), &count) in data2.counts.iter() {
        let mut found = false;

        for (&(comb_t, ref comb_cat), comb_count) in combined_data.counts.iter_mut() {
            if (data2.times[t_idx] - combined_data.times[comb_t]).abs() < EPS && cat == comb_cat {
                *comb_count += count;
                found = true;
                break;
            }
        }

        if !found {
            let idx = combined_data.times.iter() 
                .position(|&t| (t - data2.times[t_idx]).abs() < EPS)
                .unwrap_or_else(|| {
                    combined_data.times.push(data2.times[t_idx]);
                    combined_data.times.len() -1
                });
            
            combined_data.counts.insert((idx, cat.clone()), count);
        }
    }

    combined_data.num_sims = data1.num_sims + data2.num_sims;
    
    Ok(combined_data)
}


fn write_categories( 
    writer: &mut impl Write,
    times: & [f64],
    counts: &FxHashMap<(usize, Category), usize>,
) -> Result<()> {
    writeln!(writer, "time,category,count")?;
    let mut order: Vec<usize> = (0..times.len()).collect();
    order.sort_by(|&a, &b| times[a].partial_cmp(&times[b]).unwrap());

    for &t_idx in &order {
        let t= times[t_idx];
        for ((entry_t_idx, cat), count) in counts.iter() {
            if entry_t_idx == &t_idx {
                writeln!(writer, "{},{},{}", t, cat.to_key(), count)?;
            }
        }
    }
    Ok(())
}