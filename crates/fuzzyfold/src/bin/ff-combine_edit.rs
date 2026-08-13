use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::path::PathBuf;

use clap::Parser;
use anyhow::bail;
use anyhow::Result;
use rustc_hash::FxHashMap;


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
    positions: Vec<usize>,
    times: Vec<f64>,
    counts: FxHashMap<(usize, usize), usize>,
    num_sims: usize,
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
    let positions_str = combined_data.positions.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",");
    //writeln!(writer, "# num_sims={}", combined_data.num_sims)?;
    write_counts(&mut writer, &combined_data.times, &combined_data.positions, &combined_data.counts)?;

    Ok(())
}

fn parse_file(path: PathBuf) -> Result<FileData> {

    
    let content = fs::read_to_string(path)?;
    let mut data = FileData::default();

    //for line in content.lines() {
    //  if let Some(rest) = line.strip_prefix("# num_sims=") {
    //        data.num_sims = rest.parse()?; 
    //    }

    for line in content.lines() {
        if line.starts_with('#') {
            continue;
        }
        if line.starts_with("time") {
            let col: Vec<&str> = line.split(',').collect();
            for p in &col[1..] {
                data.positions.push(p.parse()?);
            }
            continue;
        }
        let cols: Vec<&str> = line.split(',').collect();
        let t:f64 = cols[0].parse()?;

        if data.times.last() != Some(&t) {
            data.times.push(t);
        }

        let t_idx = data.times.len() - 1; 

        for (i, &position) in data.positions.clone().iter().enumerate() {
            let count: usize = cols[i + 1].parse()?;
            data.counts.insert((t_idx, position), count);
        }
    }
    Ok(data)
} 


fn accumulate(data1: FileData, data2: FileData) -> Result<FileData> {

    let mut combined_data = FileData::default();

    if data1.positions != data2.positions {
        bail!("positions aren't matching")
    } else {
        combined_data.positions = data1.positions
    }

    if data1.times != data2.times {
        bail!("times aren't matching")
    } else {
        combined_data.times = data1.times 
    }

    combined_data.counts = data1.counts.clone();

    const EPS: f64 = 1e-6; 

    for (&(t_idx, pos), &count) in data2.counts.iter() {
        let mut found = false;

        for (&(comb_t, comb_pos), comb_count) in combined_data.counts.iter_mut() {
            if (data2.times[t_idx] - combined_data.times[comb_t]).abs() < EPS && pos == comb_pos {
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
            
            combined_data.counts.insert((idx, pos), count);
        }
    }

    //combined_data.num_sims = data1.num_sims + data2.num_sims;
    
    Ok(combined_data)
}


fn write_counts( 
    writer: &mut impl Write,
    times: & [f64],
    positions: &Vec<usize>,
    counts: &FxHashMap<(usize, usize), usize>,
) -> Result<()> {
    let header = positions.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
    writeln!(writer, "# positions={}", header)?;
    write!(writer, "time,")?;
    writeln!(writer, "{}", header)?;
    for (t_idx, t) in times.iter().enumerate() {
        let row = positions.iter()
            .map(|&pos| counts.get(&(t_idx, pos)).copied().unwrap_or(0).to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(writer, "{},{}", t, row)?;
    }
    Ok(())
}