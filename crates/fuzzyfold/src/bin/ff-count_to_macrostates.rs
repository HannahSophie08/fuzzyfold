//! Converts output from ff-co_stucture_count into macrostates.
//! Selects the most frequent structure for each timepoint and checks
//! whether each structure is unique. If there are structures more than 
//! once only the first is kept. Then the pairlists are converted to 
//! dotbracket vectors and each structure is safed in a seperate file with
//! the sequence and a header.
//! 
//! Input: csv file from ff-co_structure_count and fasta file that was used as
//! input for ff-co_structure_count. Optional: path to directory, where macrostate
//! files are stored. 
//! 
//! Output: macrostate files (header: macrostate index, sequence, structure in dbv)
//! 


use clap::Parser;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use ff_structure::{Pair, PairList};
use ff_structure::DotBracketVec;

use fuzzyfold::input_parsers::read_fasta_like_input;

#[derive(Debug, Parser)]
#[command(version, about = "Record structure frequencies during cotranscriptional folding")]

pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// CSV file output from ff-co_structure_count
    #[arg(short, long, value_name = "COUNTS")]
    csv: PathBuf,
 
    /// Output directory for macrostates
    #[arg(short, long, default_value = "macrostates")]
    output_dir: PathBuf,
}

struct Row {
    timepoint_idx: usize,
    time: f64,
    structure: PairList,
    count: usize,
}


fn parse_pairlist(s: &str) -> Result<PairList> {

    let mut pl = PairList::new();
    let s = s.trim();
    if s.is_empty() {
        return Ok(pl);
    }

    for part in s.split("),(") {
        let part = part.trim().trim_matches('(').trim_matches(')');
        let mut nums = part.splitn(2, ',');
        let i: u16 = nums.next().context("missing i in pair")?.trim().parse()?;
        let j: u16 = nums.next().context("missing j in pair")?.trim().parse()?;
        pl.insert(Pair::new(i, j));
    }
    Ok(pl)
}


fn parse_csv(path: &PathBuf) -> Result<Vec<Row>> {
    
    let file = File::open(path).with_context(|| format!("Could not open CSV file!"))?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();

    for (line_idx, line) in reader.lines().enumerate() {

        let line = line?;

        if line_idx == 0 {
            continue;
        }

        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let last_comma = line.rfind(',').context("No comma found in line")?;
        let count: usize = line[last_comma + 1..].trim().parse().with_context(|| format!("Could not parse count!"))?;
        let rest = &line[..last_comma];

        let mut iter = rest.splitn(3, ',');
   
        let timepoint_idx: usize = iter.next().context("Missing timepoint idx")?.trim().parse().with_context(|| format!("Could not parse timepoint_idx!"))?;
        let time: f64 = iter.next().context("Missing time")?.trim().parse().with_context(|| format!("Could not parse time!"))?;
        let structure_str = iter.next().context("Missing structure")?.trim();
        let structure = parse_pairlist(structure_str)?; 
        

        rows.push(Row {timepoint_idx, time, structure, count})

    }

    Ok(rows)

}



fn find_most_frequent_structure (rows: &[Row]) -> Vec<&Row> {

    let mut best: HashMap<usize, &Row> = HashMap::new();
    for row in rows {
        best.entry(row.timepoint_idx).or_insert(row);
    }

    let mut sorted_keys: Vec<usize> = best.keys().cloned().collect();
    sorted_keys.sort();

    sorted_keys.iter().map(|k| best[k]).collect()
    
    
}


fn remove_duplicates<'a>(rows: &[&'a Row]) -> Vec<&'a Row> {

    let mut seen: Vec<&PairList> = Vec::new();
    let mut unique = Vec::new();

    for row in rows {
        if !seen.contains(&&row.structure) {
            seen.push(&row.structure);
            unique.push(*row);
        }
    }

    unique

}


fn main() -> Result<()> {

    let cli = Cli::parse();

    let (_header, sequence, _structure) = read_fasta_like_input(&cli.input)?;

    let sequence_len = sequence.len();

    let rows = parse_csv(&cli.csv)?;

    let most_frequent = find_most_frequent_structure(&rows);

    let structures = remove_duplicates(&most_frequent);

    let mut dbvs = Vec::new();

    for row in &structures {
        let dbv = DotBracketVec::from(&row.structure);
        if dbv.is_empty() {
            continue;
        }
        dbvs.push(dbv);
    }

    std::fs::create_dir_all(&cli.output_dir)?;

    for (macro_idx, row) in structures.iter().enumerate() {
        let dbv = &dbvs[macro_idx];
        let filename = cli.output_dir.join(format!("macrostate_{}.txt", macro_idx + 1));
        let file = File::create(&filename)?;
        let mut w = BufWriter::new(file);

        writeln!(w, ">Macrostate_{}", macro_idx+1)?;
        writeln!(w, "{}", sequence)?;
        writeln!(w, "{}", dbv)?;
    }

    println!("Wrote {} mscrostate files to {:?}", structures.len(), cli.output_dir);

    Ok(())



}