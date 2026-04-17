use clap::Parser;
use anyhow::{Context, Result};
use clap::builder::Str;
use std::collections::HashMap;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use ff_structure::{Pair, PairList};
use ff_structure::DotBracketVec;

use fuzzyfold::input_parsers::read_fasta_like_input;

#[derive(Debug, Default, Parser)]
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


struct StructureCounter {
    counts: HashMap<PairList, usize>,
}


impl StructureCounter {
    fn new() -> Self {
        Self {
            counts: HashMap::new(),
        }
    }

    fn add_structure(&mut self, structure: PairList, count: usize) {
        self.counts.insert(structure, count);
    }


    fn increase(&mut self, structure: PairList, amount: usize) {
        *self.counts.entry(structure).or_insert(0) += amount;
    }

    fn get_count(&self, structure: &PairList) -> usize {
        *self.counts.get(structure).unwrap_or(&0)
    }

    fn max_structure(&self) -> Option<(&PairList, &usize)> {
        self.counts
            .iter()
            .max_by_key(|(_, count)| *count)
    }
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



fn find_most_frequent_structure (rows: &[Row]) -> (PairList, usize) {

    let mut counter = StructureCounter::new();

    
    for row in rows {
        if counter.counts.contains_key(&row.structure) {
            counter.increase(row.structure.clone(), row.count)
        } else {
            counter.add_structure(row.structure.clone(), row.count);
        }
    }
    
    let (best_structure, best_count) = counter.max_structure().unwrap();

    if best_structure.is_empty() {
        let best_structure = best_structure.clone();
        counter.counts.remove(&best_structure);
        let (new_best_structure, new_best_count) = counter.max_structure().unwrap();
        return (new_best_structure.clone(), *new_best_count);
    }

    return (best_structure.clone(), *best_count)
}




fn main() -> Result<()> {

    let cli = Cli::parse();

    let (_header, sequence, _structure) = read_fasta_like_input(&cli.input)?;

    let rows = parse_csv(&cli.csv)?;

    let (best_structure, best_count) = find_most_frequent_structure(&rows);

    let dbv = DotBracketVec::from(&best_structure);

    std::fs::create_dir_all(&cli.output_dir)?;

    let filename = cli.output_dir.join(&format!("macrostate.txt"));
    let file = File::create(&filename)?;
    let mut w = BufWriter::new(file);

    writeln!(w, ">Macrostate with count {}", best_count,)?;
    writeln!(w, "{}", sequence)?;
    writeln!(w, "{}", dbv)?;
    
    println!("Wrote macrostate with count {} to {:?}", best_count, cli.output_dir);

    Ok(())
}