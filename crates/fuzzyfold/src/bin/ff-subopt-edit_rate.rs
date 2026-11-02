use rustc_hash::FxHashMap;
use std::path::PathBuf;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::fs::{File};

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
    input: PathBuf,

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
    
    let sequence = read_full_sequence(&cli.input)?;

    let positions: Vec<usize> = if cli.all_adenosins {
        sequence.iter().enumerate()
            .filter(|(_, base)| **base == Base::A)
            .map(|(i, _)| i)
            .collect()
    } else {
        zero_indexed_positions
    };

    let alu_as: Vec<usize> = sequence.iter().enumerate()
        .filter(|(_, base)| **base == Base::A)
        .map(|(i, _)| i)
        .filter(|&i| regions.iter().any(|&(start, end)| i >= start && i <= end))
        .collect();

    let a_count = positions.len();

    if positions.is_empty() {
        bail!("No positions given: use --positions or --all-adenosins");
    }

    let (lengths, samples, counts) =
        count_per_length(&cli.input, &sequence, &positions, cli.duplex_5, cli.duplex_3)?;

    let mut edited_pos = vec![false; positions.len()];

    for row in counts.iter() {
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
    
    Ok(())
}

// first pass: the sequence line of the last (longest) block is the full sequence
fn read_full_sequence(path: &PathBuf) -> Result<NucleotideVec> {
    let reader = BufReader::new(File::open(path)?);
    let mut sequence: Option<String> = None;
    let mut next_line_is_sequence = false;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            next_line_is_sequence = true;
        } else if next_line_is_sequence {
            sequence = Some(line);
            next_line_is_sequence = false;
        }
    }
    let sequence = sequence.ok_or_else(|| anyhow::anyhow!("no sequence found in input file"))?;
    Ok(NucleotideVec::try_from(sequence.as_str())?)
}

// second pass: count one length block at a time, keep only the counts
fn count_per_length(
    path: &PathBuf,
    sequence: &NucleotideVec,
    positions: &Vec<usize>,
    duplex_5: usize,
    duplex_3: usize,
) -> Result<(Vec<usize>, Vec<usize>, Vec<Vec<usize>>)> {
    let reader = BufReader::new(File::open(path)?);
    let mut lengths: Vec<usize> = Vec::new();
    let mut samples: Vec<usize> = Vec::new();
    let mut counts: Vec<Vec<usize>> = Vec::new();
    let mut block: FxHashMap<DotBracketVec, usize> = FxHashMap::default();
    let mut next_line_is_sequence = false;

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            // a new header closes the previous block
            if !lengths.is_empty() {
                samples.push(block.values().sum());
                counts.push(edit_count(sequence, &block, positions, duplex_5, duplex_3));
                block.clear();
            }
            let len = rest.split('_').nth(1)
                .ok_or_else(|| anyhow::anyhow!("header without length: {}", line))?
                .parse::<usize>()?;
            lengths.push(len);
            next_line_is_sequence = true;
        } else if next_line_is_sequence {
            next_line_is_sequence = false;
        } else {
            if lengths.is_empty() {
                bail!("structure found before the first header");
            }
            let structure = DotBracketVec::try_from(line.as_str())?;
            *block.entry(structure).or_insert(0) += 1;
        }
    }

    // the last block has no following header
    if !lengths.is_empty() {
        samples.push(block.values().sum());
        counts.push(edit_count(sequence, &block, positions, duplex_5, duplex_3));
    }
    Ok((lengths, samples, counts))
}

// counts for one length: in how many samples each position is editable
fn edit_count(
    sequence: &NucleotideVec,
    block: &FxHashMap<DotBracketVec, usize>,
    positions: &Vec<usize>,
    duplex_5: usize,
    duplex_3: usize,
) -> Vec<usize> {
    let mut edit_counts = vec![0usize; positions.len()];
    for (db, count) in block.iter() {
        let result = gets_edited(sequence, db, positions, duplex_5, duplex_3);
        for (i, edited) in result.iter().enumerate() {
            if *edited {
                edit_counts[i] += count;
            }
        }
    }
    edit_counts
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
            
        if structure[*position] == DotBracket::Unpaired {
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
            if outer_partner[duplex_5 -1] - inner_partner[duplex_3 - 1] == duplex_len && sequence[outer_partner[0] - 1] == Base::U {
                result.push(true);
                continue;
            }
        }
     
        result.push(false);
    }
        
    return result  
}

