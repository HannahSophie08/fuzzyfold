use rustc_hash::FxHashMap;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use std::fs::{File, read_to_string};

use clap::Parser;
use anyhow::{bail, Result};

use ff_structure::{DotBracketVec, MultiPairTable, MultiStruct};
use ff_energy::NucleotideVec;

#[derive(Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT")]
    input: Vec<PathBuf>,

    #[arg(long, value_name = "START-END", num_args = 1..)]
    regions: Vec<String>,

    #[arg(long)]
    outlier: usize,

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
    outlier_pairs: Vec<FxHashMap<OutlierPairs, usize>>,
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub enum OutlierPairs {
    Alu1,
    Alu2,
    Both,
    Within,
    WithRest,
    Unpaired,
    None,
}

impl OutlierPairs {
    fn to_key(&self) -> &'static str {
        match self {
            OutlierPairs::Alu1 => "alu1",
            OutlierPairs::Alu2 => "alu2",
            OutlierPairs::Both => "both",
            OutlierPairs::Within => "within",
            OutlierPairs::WithRest => "withrest",
            OutlierPairs::Unpaired => "unpaired",
            OutlierPairs::None => "none",
        }
    }
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

    if regions.len() != 3 {
        bail!("exactly three regions are required");
    }

    if !(1..=3).contains(&cli.outlier) {
        bail!("--outlier must be 1, 2 or 3");
    }

    let outlier = cli.outlier - 1;

    let mut data = Vec::new();

    for i in cli.input.iter() {
        let file = parse_file(i.clone())?;
        data.push(file);
    }
    
    for file in data.iter_mut() {
        file.outlier_pairs = count_outlier_pairs(&file.structures, &regions, outlier)?;
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
    writeln!(writer, "# outlier={}", outlier+1)?;
    write_outlier_pairs(&mut writer, &combined_data.times, &combined_data.outlier_pairs)?;
    writer.flush()?;

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


fn count_outlier_pairs(structures: &Vec<FxHashMap<DotBracketVec, usize>>, regions: &Vec<(usize, usize)>, outlier: usize) -> Result<Vec<FxHashMap<OutlierPairs, usize>>>{

    let mut result = Vec::with_capacity(structures.len());

    let others: Vec<usize> = (0..3).filter(|&r| r != outlier).collect();
    let (alu_1, alu_2) = (others[0], others[1]);

    for timepoint_structures in structures.iter() {

        let mut counts: FxHashMap<OutlierPairs, usize> = FxHashMap::default();

        for (structure, count) in timepoint_structures.iter() {

            let mpt = MultiPairTable::try_from(structure).unwrap();
            let transcript_len = structure.len();
            let (mut with_alu1, mut with_alu2, mut within, mut with_rest, mut unpaired) = (0usize, 0usize, 0usize, 0usize, 0usize);

            for i in 0..mpt.len() {
                
                let j = match mpt.get(i) {
                    MultiStruct::Paired(j) if i < j => j,
                    _ => continue,
                };

                let region_i = regions.iter().position(|&(start, end) | i  >= start && i <= end);
                let region_j = regions.iter().position(|&(start, end) | j >= start && j <= end);
                
                match (region_i, region_j) {
                    (Some(a), Some(b)) if a == outlier && b == outlier => within += 2,
                    (Some(a), Some(b)) if (a == outlier && b == alu_1 || b == outlier && a == alu_1) => with_alu1 += 1,
                    (Some(a), Some(b)) if (a == outlier && b == alu_2 || b == outlier && a == alu_2) => with_alu2 += 1,
                    (Some(a), None) if (a == outlier) => with_rest += 1,
                    (None, Some(b)) if (b == outlier) => with_rest += 1,
                    _ => {}
                }
            }
            let (out_start, out_end) = regions[outlier];
            if transcript_len <= out_start {
                continue
            } else if transcript_len <= out_end {
                unpaired = (transcript_len - out_start) - (within + with_alu1 + with_alu2 + with_rest);
            } else {
                unpaired = (out_end - out_start + 1) - (within + with_alu1 + with_alu2 + with_rest);
            }

        
            // Minimum of 25% have to be paired
            // Ties: Both/Alu1/Alu2 > Within > WithRest
            if 4* ((within + with_alu1 + with_alu2 + with_rest)) < ((within + with_alu1 + with_alu2 + with_rest + unpaired)) {
              *counts.entry(OutlierPairs::Unpaired).or_insert(0) += *count;
            } else if within != 0 && with_alu1 == 0 && with_alu2 == 0 {
                if within >= with_rest {
                     *counts.entry(OutlierPairs::Within).or_insert(0) += *count;
                } else if with_rest > within {
                    *counts.entry(OutlierPairs::WithRest).or_insert(0) += *count;
                }
            } else if with_alu1 != 0 || with_alu2 != 0 { 
                let alus = with_alu1 + with_alu2;
                if alus >= within && alus >= with_rest {
                    if with_alu1 != 0 && with_alu2 != 0 {
                        *counts.entry(OutlierPairs::Both).or_insert(0) += *count;
                    } else if with_alu1 != 0 && with_alu2 == 0 {
                        *counts.entry(OutlierPairs::Alu1).or_insert(0) += *count;
                    } else if with_alu1 == 0 && with_alu2 != 0 {
                        *counts.entry(OutlierPairs::Alu2).or_insert(0) += *count;
                    } 
                } else if within >= with_rest {
                     *counts.entry(OutlierPairs::Within).or_insert(0) += *count;
                } else if with_rest > within {
                    *counts.entry(OutlierPairs::WithRest).or_insert(0) += *count;
                }
            } else if within == 0 && with_alu1 == 0 && with_alu2 == 0 {
                if with_rest != 0 {
                    *counts.entry(OutlierPairs::WithRest).or_insert(0) += *count;
                }
            }

            /* 
            // Ties: Both/Alu1/Alu2 > Within > WithRest > Unpaired
            if within != 0 && with_alu1 == 0 && with_alu2 == 0 {
                if within >= with_rest && within >= unpaired {
                     *counts.entry(OutlierPairs::Within).or_insert(0) += *count;
                } else if with_rest > within && with_rest >= unpaired {
                    *counts.entry(OutlierPairs::WithRest).or_insert(0) += *count;
                } else if unpaired > within && unpaired > with_rest {
                    *counts.entry(OutlierPairs::Unpaired).or_insert(0) += *count;
                }
            } else if with_alu1 != 0 || with_alu2 != 0 { 
                let alus = with_alu1 + with_alu2;
                if alus >= within && alus >= with_rest && alus >= unpaired {
                    if with_alu1 != 0 && with_alu2 != 0 {
                        *counts.entry(OutlierPairs::Both).or_insert(0) += *count;
                    } else if with_alu1 != 0 && with_alu2 == 0 {
                        *counts.entry(OutlierPairs::Alu1).or_insert(0) += *count;
                    } else if with_alu1 == 0 && with_alu2 != 0 {
                        *counts.entry(OutlierPairs::Alu2).or_insert(0) += *count;
                    } 
                } else if within >= with_rest && within >= unpaired {
                     *counts.entry(OutlierPairs::Within).or_insert(0) += *count;
                } else if with_rest > within && with_rest >= unpaired {
                    *counts.entry(OutlierPairs::WithRest).or_insert(0) += *count;
                } else if unpaired > within && unpaired > with_rest {
                    *counts.entry(OutlierPairs::Unpaired).or_insert(0) += *count;
                }
            } else if within == 0 && with_alu1 == 0 && with_alu2 == 0 {
                if with_rest >= unpaired {
                    *counts.entry(OutlierPairs::WithRest).or_insert(0) += *count;
                } else if unpaired > with_rest {
                    *counts.entry(OutlierPairs::Unpaired).or_insert(0) += *count;
                }
            }

            
            if within != 0 && with_alu1 == 0 && with_alu2 == 0{
                *counts.entry(OutlierPairs::Within).or_insert(0) += *count;
            } else if with_alu1 != 0 && with_alu2 != 0 { 
                if within > (with_alu1 + with_alu2) {
                    *counts.entry(OutlierPairs::Mixed).or_insert(0) += *count;
                } else {
                    *counts.entry(OutlierPairs::Both).or_insert(0) += *count;
                }
            } else if with_alu1 != 0 && with_alu2 == 0 { 
                if within > with_alu1 {
                    *counts.entry(OutlierPairs::Mixed).or_insert(0) += *count;
                } else {
                    *counts.entry(OutlierPairs::Alu1).or_insert(0) += *count;
                }
            } else if with_alu1 == 0 && with_alu2 != 0 { 
                if within > with_alu2 {
                    *counts.entry(OutlierPairs::Mixed).or_insert(0) += *count;
                } else {
                    *counts.entry(OutlierPairs::Alu2).or_insert(0) += *count;
                }
            } else if within == 0 && with_alu1 == 0 && with_alu2 == 0 {
                *counts.entry(OutlierPairs::None).or_insert(0) += *count;
            }
            

            // Present
            if within != 0 && with_alu1 == 0 && with_alu2 == 0{
                *counts.entry(OutlierPairs::Within).or_insert(0) += *count;
            } else if with_alu1 != 0 && with_alu2 != 0 { 
                *counts.entry(OutlierPairs::Both).or_insert(0) += *count;
            } else if with_alu1 != 0 && with_alu2 == 0 { 
                *counts.entry(OutlierPairs::Alu1).or_insert(0) += *count;
            } else if with_alu1 == 0 && with_alu2 != 0 { 
                *counts.entry(OutlierPairs::Alu2).or_insert(0) += *count;
            } else if within == 0 && with_alu1 == 0 && with_alu2 == 0 {
                *counts.entry(OutlierPairs::None).or_insert(0) += *count;
            }
            */
            
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
    combined_data.outlier_pairs = data1.outlier_pairs.clone();

    for (t_idx, outlier_pairs) in data2.outlier_pairs.into_iter().enumerate() {
        for (cat, count) in outlier_pairs.into_iter() {
            *combined_data.outlier_pairs[t_idx].entry(cat).or_insert(0) += count;
        }
    }

    combined_data.num_sims = data1.num_sims + data2.num_sims;

    Ok(combined_data)
}
   

fn write_outlier_pairs( 
    writer: &mut impl Write,
    times: &[f64],
    counts: &[FxHashMap<OutlierPairs, usize>],
) -> Result<()> {
    writeln!(writer, "time,category,count")?;
    for (t_idx, t) in times.iter().enumerate() {
        for (cat, count) in counts[t_idx].iter() {
            writeln!(writer, "{},{},{}", t, cat.to_key(), count)?;
        }
    }
    Ok(())
}