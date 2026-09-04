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

    #[arg(short, long, value_name = "POSITION", num_args = 1..)]
    positions: Vec<usize>,

    #[arg(short = '5', long, default_value_t = 1)]
    duplex_5: usize,

    #[arg(short = '3', long, default_value_t = 1)]
    duplex_3: usize,

    #[arg(short, long)]
    all_adenosins: bool,

    #[arg(long, value_name = "START-END", num_args = 1..)]
    regions: Vec<String>,
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

    let mut positions: Vec<usize> = if cli.all_adenosins {
        let sequence = data[0].sequence.as_ref()
            .ok_or_else(|| anyhow::anyhow!("missing sequence in input file"))?;
        let adenines = sequence.iter().enumerate()
            .filter(|(_, base)| **base == Base::A)
            .map(|(i, _)| i);
            adenines.collect()
    } else {
        zero_indexed_positions
    };

    let sequence = data[0].sequence.as_ref()
            .ok_or_else(|| anyhow::anyhow!("missing sequence in input file"))?;
    let alu_as: Vec<usize> = sequence.iter().enumerate()
        .filter(|(_, base)| **base == Base::A)
        .map(|(i, _)| i)
        .filter(|&i| regions.iter().any(|&(start, end)| i >= start && i <= end))
        .collect();

    let a_count = positions.len();

    if positions.is_empty() {
        bail!("No positions given: use --positions or --all-adenosins");
    }

    for file in data.iter_mut() {
        let sequence = file.sequence.as_ref().ok_or_else(|| anyhow::anyhow!("missing sequence in input file"))?;
        file.counts = edit_count(sequence, &file.structures, &positions, cli.duplex_5, cli.duplex_3);
    }
    
    let mut data_iter = data.into_iter();
    
    let mut combined_data = data_iter.next().unwrap();

    for file in data_iter {
        combined_data = accumulate(combined_data, file)?;
    }

    let mut edited_pos = vec![false; positions.len()];

    for row in combined_data.counts.iter() {
        for (i, &count) in row.iter().enumerate() {
            if count > 0 {
                edited_pos[i] = true; 
            }
        }
    }

    let mut edited_count = 0;
    let mut alu_count = 0;
    let mut alu_positions = Vec::new();
    let mut alu_specific = vec![0; regions.len()];

    for (i, &is_edited) in edited_pos.iter().enumerate() {
        if is_edited {
            edited_count += 1;
            if regions.iter().position(|&(start, end)| positions[i] >= start && positions[i] <= end).is_some() {
                alu_count += 1;
                alu_positions.push(positions[i]);
            } 
            for (r, &(start, end)) in regions.iter().enumerate() {
                if positions[i] >= start && positions[i] <= end {
                    alu_specific[r] += 1;
                }
            }
        } 
    }

    println!("Total As: {}", a_count);
    println!("As in Alus: {}", alu_as.len());
    println!("Edited As: {}", edited_count);
    println!("Edited As in Alus: {}", alu_count);
    println!("Edited As in Alu_1: {}", alu_specific[0]);
    println!("Edited As in Alu_2: {}", alu_specific[1]);
    println!("Edited As in Alu_3: {}", alu_specific[2]);
    println!("Edited positions in Alus: {:?}", alu_positions);
    let sequence = combined_data.sequence.unwrap();
   

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

fn edit_count (sequence: &NucleotideVec, structures: &Vec<FxHashMap<DotBracketVec, usize>> , positions: &Vec<usize>, duplex_5: usize, duplex_3: usize) -> Vec<Vec<usize>> {

    structures.iter().map(|counts| {
        let mut edit_counts = vec![0usize; positions.len()];
        for (db, count) in counts.iter() { 
            let result = gets_edited(sequence, db, positions, duplex_5, duplex_3);
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
// 5' and 3' neighbors are paired
fn gets_edited (sequence: &NucleotideVec, structure: &DotBracketVec, positions: &Vec<usize>, duplex_5: usize, duplex_3: usize) -> Vec<bool> {

    let mut result = Vec::new();
    let pt = PairTable::try_from(structure).unwrap();

    'outer_pos: for position in positions.iter() {

        if *position < duplex_5 || *position + duplex_3 >= structure.len() {
            result.push(false);
            continue;
        }
            
        if structure[*position] != DotBracket::Unpaired {
            result.push(false);
            continue;
        }

        let mut outer = Vec::new();
        for i in 1..=duplex_5 {
            outer.push((*position - i) as usize);
        }

        let mut inner = Vec::new();
        for i in 1..=duplex_3 {
            inner.push((*position + i) as usize);
        }

        let mut outer_partner = Vec::new();
        for o in outer {
            if !pt[o].is_some() {
                result.push(false );
                continue 'outer_pos; 
            } else {
                outer_partner.push(pt[o].unwrap() as usize);
            }
        } 

        let mut inner_partner = Vec::new();
        for i in inner {
            if !pt[i].is_some() {
                result.push(false );
                continue 'outer_pos; 
            } else {
                inner_partner.push(pt[i].unwrap() as usize);
            }
        }
        
        for idx in 0..outer_partner.len() {
            if idx + 1 < outer_partner.len() {
                if outer_partner[idx] +  1 != outer_partner[idx + 1] {
                    result.push(false );
                    continue 'outer_pos; 
                }
            }
        }

        for idx in 0..inner_partner.len() {
            if idx + 1 < inner_partner.len() {
                if inner_partner[idx+1] + 1 != inner_partner[idx] {
                    result.push(false );
                    continue 'outer_pos; 
                }
            }
        }

        let duplex_len = duplex_5 + duplex_3;

        if outer_partner[duplex_5 -1] > inner_partner[duplex_3 - 1] {
            if outer_partner[duplex_5 -1] - inner_partner[duplex_3 - 1] == duplex_len && sequence[outer_partner[0] - 1] == Base::C {
                result.push(true);
                continue;
            }
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
