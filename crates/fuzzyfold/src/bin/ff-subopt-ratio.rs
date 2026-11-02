use rustc_hash::FxHashMap;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use std::fs::{File, read_to_string};

use clap::Parser;
use anyhow::{bail, Result};

use ff_structure::DotBracketVec;
use ff_structure::PairList;
use ff_energy::NucleotideVec;

#[derive(Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT")]
    input: PathBuf,

    #[arg(long, value_name = "START-END", num_args = 1..)]
    regions: Vec<String>,

    #[arg(long)]
    outlier: usize,

    #[arg(long, default_value_t = 1)]
    num_samples: usize,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,
}

#[derive(Clone, Eq, Hash, PartialEq)]
pub enum OutlierPairs {
    Alu1,
    Alu2,
    Both,
    Within,
    WithRest,
    Unpaired,
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
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct AluPairs {
    alu1: usize,
    alu2: usize,
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

    let (sequence, lengths, structures) = parse_file(cli.input)?;
    
    let alu_pairs = count_outlier_pairs(&structures, &regions, outlier)?;

    let csv_path = cli.output.with_extension("csv");
    
    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
    let regions_str = regions.iter()
        .map(|(a, b)| format!("{}-{}", a, b))
        .collect::<Vec<_>>()
        .join(",");
    writeln!(writer, "# num_samples={}", cli.num_samples)?;
    writeln!(writer, "# regions={}", regions_str)?;
    writeln!(writer, "# outlier={}", outlier+1)?;
    write_alu_pairs(&mut writer, &lengths, &alu_pairs)?;
    writer.flush()?;

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



fn count_outlier_pairs(structures: &Vec<FxHashMap<DotBracketVec, usize>>, regions: &Vec<(usize, usize)>, outlier: usize) -> Result<Vec<FxHashMap<AluPairs, usize>>> {

    let mut result = Vec::with_capacity(structures.len());

    let others: Vec<usize> = (0..3).filter(|&r| r != outlier).collect();
    let (alu_1, alu_2) = (others[0], others[1]);

    for timepoint_structures in structures.iter() {

        let mut counts: FxHashMap<AluPairs, usize> = FxHashMap::default();

        for (structure, count) in timepoint_structures.iter() {

            let pl = PairList::try_from(structure)?;
            let transcript_len = structure.len();
            let (mut with_alu1, mut with_alu2, mut within, mut with_rest, mut unpaired) = (0usize, 0usize, 0usize, 0usize, 0usize);
            for (i, j) in pl.iter() {
                let region_i = regions.iter().position(|&(start, end) | (*i as usize) >= start && (*i as usize) <= end);
                let region_j = regions.iter().position(|&(start, end) | (*j as usize) >= start && (*j as usize) <= end);
                
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
        
            if 4 * ((with_alu1 + with_alu2)) < 3 * ((within + with_alu1 + with_alu2 + with_rest + unpaired)) {
              continue
            }
            
            if with_alu1 != 0 || with_alu2 != 0 { 
                let alus = with_alu1 + with_alu2;
                if alus >= within && alus >= with_rest {
                    let key = AluPairs { alu1: with_alu1, alu2: with_alu2 };
                    *counts.entry(key).or_insert(0) += *count;
                }
            }
        }
        result.push(counts); 

    }
    Ok(result)
}

fn write_alu_pairs( 
    writer: &mut impl Write,
    lengths: &[usize],
    counts: &[FxHashMap<AluPairs, usize>],
) -> Result<()> {
    writeln!(writer, "length,alu1,alu2,count")?;
    for (l_idx, l) in lengths.iter().enumerate() {
        for (key, count) in counts[l_idx].iter() {
            writeln!(writer, "{},{},{},{}", l, key.alu1, key.alu2, count)?;
        }
    }
    Ok(())
}