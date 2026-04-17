use clap::Parser;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use ff_structure::{PairList, PairTable};
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

        let parts: Vec<&str> = line.splitn(4, ',').collect();
        if parts.len() < 4 {
            continue;
        }

        let timepoint_idx: usize = parts[0].trim().parse().with_context(|| format!("Could not parse timepoint_idx!"))?;
        let time: f64 = parts[1].trim().parse().with_context(|| format!("Could not parse time!"))?;
        let dbv = DotBracketVec::try_from(parts[2].trim()).with_context(|| "Invalid dot-bracket")?;
        let pt = PairTable::try_from(&dbv)?;
        let structure = PairList::from(&pt);
        let count: usize = parts[3].trim().parse().with_context(|| format!("Could not parse count!"))?;

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