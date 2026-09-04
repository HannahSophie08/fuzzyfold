use rustc_hash::FxHashMap;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use std::fs::{File, read_to_string};

use clap::Parser;
use anyhow::{bail, Result};

use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_structure::PairList;
use ff_energy::NucleotideVec;
use ff_energy::Base;

use fuzzyfold::category::Category;

#[derive(Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    #[arg(short, long, value_name = "POSITION", num_args = 1..)]
    positions: Vec<usize>,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[arg(long, default_value_t = 1)]
    num_samples: usize,

    #[arg(short = '5', long, default_value_t = 1)]
    duplex_5: usize,

    #[arg(short = '3', long, default_value_t = 1)]
    duplex_3: usize,

    #[arg(short, long)]
    all_adenosins: bool,

    #[arg(long, value_name = "START-END", num_args = 1..)]
    regions: Vec<String>,

    #[arg(short, long)]
    interactions: bool,
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
    
    let (sequence, lengths, structures) = parse_file(cli.input)?;
    
    if cli.interactions {
        let categories = count_categories(&structures, &regions)?;
    
        let csv_path = cli.output.with_extension("csv");
        
        let mut writer = BufWriter::new(File::create(csv_path.clone())?);
        write_categories(&mut writer, &lengths, &categories)?;
        
    } else {
        let mut positions: Vec<usize> = if cli.all_adenosins {
            let adenines = sequence.iter().enumerate()
                .filter(|(_, base)| **base == Base::A)
                .map(|(i, _)| i);

                if regions.is_empty() {
                    adenines.collect()
                } else {
                    adenines.filter(|&i| regions.iter().any(|&(start, end)| i >= start && i <= end))
                    .collect()
                }
        } else {
            zero_indexed_positions
        };

        if positions.is_empty() {
            bail!("No positions given: use --positions or --all-adenosins");
        }

        let mut counts = edit_count(&sequence, &structures, &positions, cli.duplex_5, cli.duplex_3);
        if cli.all_adenosins {
            if regions.is_empty() {
                counts = counts.iter().map(|inner| vec![inner.iter().sum()]).collect();
            } else {
                counts = counts.iter().map(|inner| {
                    let mut region_counts = vec![0usize; regions.len()];
                    for (i, &count) in inner.iter().enumerate() {
                        let r = regions.iter().position(|&(start, end)| positions[i] >= start && positions[i] <= end).unwrap(); 
                        region_counts[r] += count; 
                    }
                    region_counts
                }).collect();
            }
        }

        if cli.all_adenosins {
            if regions.is_empty() {
                let len_positions = positions.len();
                positions = Vec::new();
                positions.push(len_positions -1);
            } else {
                positions = Vec::new();
                for (idx, _region) in regions.iter().enumerate() {
                    positions.push(idx);
                }
                println!("regions, positions: {:?}", positions);
            }
        }

        let csv_path = cli.output.with_extension("csv");
    
        let mut writer = BufWriter::new(File::create(csv_path.clone())?);
        writeln!(writer, "# sequence={}", sequence)?;
        writeln!(writer, "# num_samples={}", cli.num_samples)?;
        writeln!(writer, "# duplex_5={}", cli.duplex_5)?;
        writeln!(writer, "# duplex_3={}", cli.duplex_3)?;
        let one_indexed_positions: Vec<usize> = positions.iter().map(|p| p + 1).collect();
        write_counts(&mut writer, &lengths, &one_indexed_positions, &counts)?;

    }
    Ok(())
}

fn parse_file(path: PathBuf) -> Result<(NucleotideVec, Vec<usize>, Vec<FxHashMap<DotBracketVec, usize>>)> {

    let content = read_to_string(path)?;
    let mut lengths = Vec::new();
    let mut structures = Vec::new();
    let mut sequence: Option<NucleotideVec> = None;
    let mut next_line_is_sequence = false;


    for line in content.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            let len = rest.split('_').nth(1).unwrap().parse()?;
            lengths.push(len);
            structures.push(FxHashMap::default());
            next_line_is_sequence = true;
        } else if next_line_is_sequence {
            sequence = Some(NucleotideVec::try_from(line)?);
            next_line_is_sequence = false;
        } else {
            let structure = DotBracketVec::try_from(line)?;
            *structures.last_mut().unwrap().entry(structure).or_insert(0) += 1;
        }
    }
    Ok((sequence.ok_or_else(|| anyhow::anyhow!("no sequence found in input file"))?, lengths, structures))
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

fn write_counts( 
    writer: &mut impl Write,
    lengths: &Vec<usize>,
    positions: &Vec<usize>,
    counts: &Vec<Vec<usize>>,
) -> Result<()> {
    let header = positions.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
    writeln!(writer, "# positions={}", header)?;
    write!(writer, "time,")?;
    writeln!(writer, "{}", header)?;
    for (t_idx, t) in lengths.iter().enumerate() {
        let row = positions.iter().enumerate()
            .map(|(i, _pos)| counts[t_idx][i].to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(writer, "{},{}", t, row)?;
    }
    Ok(())
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


fn write_categories( 
    writer: &mut impl Write,
    lengths: &Vec<usize>,
    counts: &[FxHashMap<Category, usize>],
) -> Result<()> {
    writeln!(writer, "length,category,count")?;
    for (idx, l) in lengths.iter().enumerate() {
        for (cat, count) in counts[idx].iter() {
            writeln!(writer, "{},{},{}", l, cat.to_key(), count)?;
        }
    }
    Ok(())
}

