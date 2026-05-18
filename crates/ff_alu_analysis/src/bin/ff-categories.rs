use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::path::PathBuf;

use ff_structure::DotBracketVec;
use ff_structure::PairList;
use ff_structure::PairTable;
use clap::Parser;
use anyhow::Result;
use rustc_hash::FxHashMap;

use ff_alu_analysis::category::Category;

#[derive(Debug, Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[arg(long, value_name = "START-END", num_args = 1..)]
    regions: Vec<String>
}



fn main() -> Result<()> {
    let cli = Cli::parse();
   
    let csv_path = cli.output.with_extension("csv");

    // Parse regions
    let regions: Vec<(usize, usize)> = cli.regions.iter()
        .map(|r| {
            let (a, b) = r.split_once('-').expect("region must be in format START-END");
            (a.parse::<usize>().expect("invalid start"), b.parse::<usize>().expect("invalid end"))
        })
        .collect();

    let content = std::fs::read_to_string(&cli.input)?;
    let mut lines = content.lines().filter(|l| !l.is_empty());

    // skip header and sequence
    lines.next(); // >STK4_IVT
    lines.next(); // GGGCGAA...

    // all remaining lines are structures
    let dot_bracket_structures: Vec<DotBracketVec> = lines
    .map(|l| {
        let pt = PairTable::try_from(l).expect("invalid dot-bracket line");
        DotBracketVec::from(&pt)
    })
    .collect();

println!("Loaded {} structures", dot_bracket_structures.len());
    let categories =  count_categories(&dot_bracket_structures, &regions)?;

    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
    
    write_categories(&mut writer, &categories)?;


    Ok(())
}


fn count_categories(structures: &[DotBracketVec], regions: &Vec<(usize, usize)>) -> Result<FxHashMap<Category, f64>>{

    let mut counts: FxHashMap<Category, usize> = FxHashMap::default();
    let mut total_pairs = 0;

    let mut skipped = 0;
    for structure in structures.iter() {

        let pl = PairList::try_from(structure)?;

        for (i, j) in pl.iter() {
            let region_i = regions.iter().position(|&(start, end) | (*i as usize) >= start && (*i as usize) <= end);
            let region_j = regions.iter().position(|&(start, end) | (*j as usize) >= start && (*j as usize) <= end);

            let category = match (region_i, region_j) {
                (Some(a), Some(b)) if a == b => Category::Within(a),
                (Some(a), Some(b)) => Category::Between(a.min(b), a.max(b)),
                (Some(a), None) => Category::WithRest(a),
                (None, Some(b)) => Category::WithRest(b),
                (None, None) => { skipped += 1; continue },
            };

            *counts.entry(category).or_insert(0) += 1;
            total_pairs += 1;
        
        }
    }
    eprintln!("Skipped {} pairs (both outside regions)", skipped);
    eprintln!("Counted {} pairs", total_pairs);


    let percentages = counts.into_iter()
        .map(|(cat, count)| (cat, count as f64 / total_pairs as f64 * 100.0))
        .collect();
    
    Ok(percentages)
}


fn write_categories( 
    writer: &mut impl Write,
    categories: &FxHashMap<Category, f64>,
) -> Result<()> {
    writeln!(writer, "category,value")?;
    for (cat, val) in categories {
        writeln!(writer, "{},{}", cat.to_key(), val)?;
    }
    Ok(())
}