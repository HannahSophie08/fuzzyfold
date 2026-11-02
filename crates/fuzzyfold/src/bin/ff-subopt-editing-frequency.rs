use rustc_hash::FxHashMap;
use std::path::PathBuf;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::fs::File;

use clap::Parser;
use anyhow::{bail, Result};

use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::Base;

#[derive(Parser)]
#[command(version, about = "Mean equilibrium editable fraction per A over all transcript lengths (RNAsubopt samples).")]
pub struct Cli {
    /// RNAsubopt sample file (one '>len_N' block per transcript length)
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,
}

pub struct Pattern {
    name: &'static str,
    paired: bool,
    opposite: Option<Base>,
    duplex_3: usize,
    duplex_5: usize,
}

fn main() -> Result<()> {

    let cli = Cli::parse();

    let sequence = read_full_sequence(&cli.input)?;

    let adenosines: Vec<usize> = sequence.iter().enumerate()
        .filter(|(_, base)| **base == Base::A)
        .map(|(i, _)| i)
        .collect();

    let ac_3: Pattern = Pattern {
        name: "ac_3",
        paired: false,
        opposite: Some(Base::C),
        duplex_3: 1,
        duplex_5: 1,
    };

    let ac_5: Pattern = Pattern {
        name: "ac_5",
        paired: false,
        opposite: Some(Base::C),
        duplex_3: 2,
        duplex_5: 2,
    };

    let ac_7: Pattern = Pattern {
        name: "ac_7",
        paired: false,
        opposite: Some(Base::C),
        duplex_3: 3,
        duplex_5: 3,
    };
    
    let au_3: Pattern = Pattern {
        name: "au_3",
        paired: true,
        opposite: Some(Base::U),
        duplex_3: 1,
        duplex_5: 1,
    };

    let au_5: Pattern = Pattern {
        name: "au_5",
        paired: true,
        opposite: Some(Base::U),
        duplex_3: 2,
        duplex_5: 2,
    };

    let au_7: Pattern = Pattern {
        name: "au_7",
        paired: true,
        opposite: Some(Base::U),
        duplex_3: 3,
        duplex_5: 3,
    };

    let ag_3: Pattern = Pattern {
        name: "ag_3",
        paired: false,
        opposite: Some(Base::G),
        duplex_3: 1,
        duplex_5: 1,
    };

    let ag_5: Pattern = Pattern {
        name: "ag_5",
        paired: false,
        opposite: Some(Base::G),
        duplex_3: 2,
        duplex_5: 2,
    };

    let ag_7: Pattern = Pattern {
        name: "ag_7",
        paired: false,
        opposite: Some(Base::G),
        duplex_3: 3,
        duplex_5: 3,
    };

    let bulge_3: Pattern = Pattern {
        name: "bulge_3",
        paired: false,
        opposite: None,
        duplex_3: 1,
        duplex_5: 1,
    };

    let bulge_5: Pattern = Pattern {
        name: "bulge_5",
        paired: false,
        opposite: None,
        duplex_3: 2,
        duplex_5: 2,
    };

    let bulge_7: Pattern = Pattern {
        name: "bulge_7",
        paired: false,
        opposite: None,
        duplex_3: 3,
        duplex_5: 3,
    };

    let patterns: Vec<Pattern> = vec![ac_3, ac_5, ac_7, au_3, au_5, au_7, ag_3, ag_5, ag_7, bulge_3, bulge_5, bulge_7];

    let (lengths, samples, counts) =
        count_per_length(&cli.input, &sequence, &adenosines, &patterns)?;

    if lengths.is_empty() {
        bail!("no length blocks found in input file");
    }

    let seq_str = sequence.to_string();

    for (k, pattern) in patterns.iter().enumerate() {
        let mut mean_frac = vec![0.0f64; adenosines.len()];
        for (length_counts, &n) in counts.iter().zip(&samples) {
            for (i, &c) in length_counts[k].iter().enumerate() {
                mean_frac[i] += c as f64 / n as f64;
            }
        }
        for f in mean_frac.iter_mut() {
            *f /= lengths.len() as f64;
        }
        let mut value: Vec<Option<f64>> = vec![None; sequence.len()];
        for (j, &p) in adenosines.iter().enumerate() {
            value[p] = Some(mean_frac[j]);
        }

        let path = cli.output.with_file_name(format!("{}_{}.csv", cli.output.file_stem().unwrap().to_string_lossy(), pattern.name));
         let mut writer = BufWriter::new(File::create(&path)?);
        writeln!(writer, "# sequence={}", seq_str)?;
        writeln!(writer, "# num_lengths={}", lengths.len())?;
        writeln!(writer, "# duplex_5={}", pattern.duplex_5)?;
        writeln!(writer, "# duplex_3={}", pattern.duplex_3)?;
        writeln!(writer, "# opposite={:?}", pattern.opposite)?;
        writeln!(writer, "# paired={}", pattern.paired)?;
        writeln!(writer, "position\tnucleotide\tmean_editable_fraction")?;

        for (i, nt) in seq_str.chars().enumerate() {
            match value[i] {
                Some(v) => writeln!(writer, "{}\t{}\t{}", i + 1, nt, v)?,
                None => writeln!(writer, "{}\t{}\tNA", i + 1, nt)?,
            }
        }
        writer.flush()?;
        println!("Wrote {}", path.display());
    }
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
    patterns: &[Pattern],
) -> Result<(Vec<usize>, Vec<usize>, Vec<Vec<Vec<usize>>>)> {

    let reader = BufReader::new(File::open(path)?);
    let mut lengths: Vec<usize> = Vec::new();
    let mut samples: Vec<usize> = Vec::new();
    let mut counts: Vec<Vec<Vec<usize>>> = Vec::new();
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
                counts.push(edit_count(sequence, &block, positions, patterns));
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
        counts.push(edit_count(sequence, &block, positions, patterns));
    }
    Ok((lengths, samples, counts))
}

// counts for one length: in how many samples each position is editable
fn edit_count(
    sequence: &NucleotideVec,
    block: &FxHashMap<DotBracketVec, usize>,
    positions: &Vec<usize>,
    patterns: &[Pattern],
) -> Vec<Vec<usize>> {
    let mut edit_counts = vec![vec![0usize; positions.len()]; patterns.len()];
    for (db, count) in block.iter() {
        let pt = PairTable::try_from(db).expect("invalud structure in input file");
        for (k, pattern) in patterns.iter().enumerate() {
            let result = position_check(sequence, db, &pt, positions, pattern);
            for (i, edited) in result.iter().enumerate() {
                if *edited {
                    edit_counts[k][i] += count;
                }
            }
        }
    }
    edit_counts
}

fn position_check (sequence: &NucleotideVec, structure: &DotBracketVec, pt: &PairTable, positions: &Vec<usize>, pattern: &Pattern) -> Vec<bool> {

    let mut result = Vec::new();

    'outer_pos: for position in positions.iter() {

        if *position < pattern.duplex_5 || *position + pattern.duplex_3 >= structure.len() {
            result.push(false);
            continue;
        }
        
        if !pattern.paired {
            if structure[*position] != DotBracket::Unpaired {
                result.push(false);
                continue;
            }
        } else {
            if structure[*position] == DotBracket::Unpaired {
                result.push(false);
                continue;
            }
        }

        let mut outer = Vec::new();
        for i in 1..=pattern.duplex_5 {
            outer.push((*position - i) as usize);
        }

        let mut inner = Vec::new();
        for i in 1..=pattern.duplex_3 {
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

        let duplex_len = pattern.duplex_5 + pattern.duplex_3;

        if !pattern.opposite.is_none() {
            if outer_partner[pattern.duplex_5 -1] > inner_partner[pattern.duplex_3 - 1] {
            
                if outer_partner[pattern.duplex_5 -1] - inner_partner[pattern.duplex_3 - 1] == duplex_len && sequence[outer_partner[0] - 1] == pattern.opposite.unwrap() {
                    result.push(true);
                    continue;
                }
            }
        } else {
            if outer_partner[pattern.duplex_5 -1] > inner_partner[pattern.duplex_3 - 1] {
                if outer_partner[pattern.duplex_5 -1] - inner_partner[pattern.duplex_3 - 1] == duplex_len - 1{
                    result.push(true);
                    continue;
                }
            } 
        }
     
        result.push(false);
    }
        
    return result  
}