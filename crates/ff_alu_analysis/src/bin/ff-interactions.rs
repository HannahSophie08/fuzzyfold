use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::sync::Arc;
use std::path::Path;
use std::path::PathBuf;

use ff_structure::DotBracketVec;
use ff_structure::PairList;
use rayon::prelude::*;
use rand::rng;
use colored::*;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use rustc_hash::FxHashMap;

use ff_structure::PairTable;
use ff_energy::EnergyModel;
use ff_kinetics::RateModel;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::shift_policy::*;
use ff_kinetics::SSA;

use ff_alu_analysis::input_parsers::read_cotr_input;
use ff_alu_analysis::energy_parsers::EnergyModelArguments;
use ff_alu_analysis::kinetics_parsers::RateModelArguments;
use ff_alu_analysis::kinetics_parsers::TimelineParameters;
use ff_alu_analysis::category::Category;

#[derive(Debug, Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[arg(long, value_name = "START-END", num_args = 1..)]
    regions: Vec<String>,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimelineParameters,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}


fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.simulation.validate()?;

    // --- Build simulator ---
    let emodel = Arc::new(cli.energy.build_model());
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) = read_cotr_input(&cli.input, is_rna)?;
    let pairings = PairTable::try_from(&structure)?;

    let regions: Vec<(usize, usize)> = cli.regions.iter()
    .map(|r| {
        let (a, b) = r.split_once('-').expect("region must be in format START-END");
        (a.parse::<usize>().expect("invalid start"), b.parse::<usize>().expect("invalid end"))
    })
    .collect();

    let csv_path = cli.output.with_extension("csv");

    // --- Check input ---
    if cli.simulation.t_ext.is_none() {
        panic!("Error:'t-ext' (extension time) missing for cotranscriptional simulation!"); 
    }

    if structure.len() == sequence.len() {
        panic!("Error: Full length structure given, for cotranscriptional simulation give either 
        shorter start structure or no structure!");
    }

    if let Some(h) = header {
        println!("{}", h.yellow());
    }

    println!("{}", sequence);

    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    let times = cli.simulation.get_output_times(&sequence, Some(&structure)); // build times vector for output times 

    let mut sim_times = vec![cli.simulation.t_ext.unwrap(); sequence.len() - structure.len()];
    sim_times.push(cli.simulation.t_end);

    let all_structures: Vec<Vec<DotBracketVec>> =
        match (rmodel.k3ws().is_some(), rmodel.k4ws().is_some()) {
            (false, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, NoShift))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
            (true, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
            (false, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, FourWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
            (true, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeAndFour))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
        };

    
    let sorted_structures = sort_by_timepoint(&times, &all_structures, cli.num_sims);

    let categories =  count_categories(&sorted_structures, &regions)?;

    let co_time = cli.simulation.t_ext.unwrap() * (sequence.len() - structure.len()) as f64;
    let post_time = co_time + cli.simulation.t_end;

    println!("{}", "Finished simulations!".red());

    let mut existing_num_sims: usize = 0; 
    let mut existing_counts: FxHashMap<(usize, Category), usize> = FxHashMap::default();

    let accumulate = if Path::new(&csv_path).exists() {
        
        let exisiting_content = fs::read_to_string(&csv_path)?;
        let mut matches = true; 

        for line in exisiting_content.lines() {
            if let Some(rest) = line.strip_prefix("# t_ext=") {
                let existing_t_ext: f64 = rest.parse()?; 
                if existing_t_ext != cli.simulation.t_ext.unwrap() {
                    println!("Error: t-ext of existing file doesn't match! Simualtions can not be accumulated");
                    matches = false; 
                    break 
                }         
            }
            if let Some(rest) = line.strip_prefix("# t_end=") {
                let existing_t_end: f64 = rest.parse()?;
                if existing_t_end != post_time {
                    println!("Error: t-end of existing file doesn't match! Simulations can not be accumulated");
                    matches = false; 
                    break
                }
            }   
            if let Some(rest) = line.strip_prefix("# t_split=") {
                let existing_t_split: f64 = rest.parse()?;
                if existing_t_split != co_time {
                    println!("Error: t-end of existing file doesn't match! Simulations can not be accumulated");
                    matches = false; 
                    break
                } 
            }
            if let Some(rest) = line.strip_prefix("# regions=") {
                let mut existing_regions: Vec<(usize, usize)> = Vec::new();
                for piece in rest.split(',') {
                    let (a, b) = piece.split_once('-').expect("region must be in format START-END");
                    let a: usize = a.parse()?;
                    let b: usize = b.parse()?;
                    existing_regions.push((a, b));
                }
                if existing_regions != regions {
                    println!("Error: regions of existing file don't match! Simulations can not be accumulated");
                    matches = false; 
                    break
                } 
            }

            if let Some(rest) = line.strip_prefix("# num_sims=") {
                existing_num_sims = rest.parse()?; 
            }  
        }

        if matches {

            let mut existing_times: Vec<f64> = Vec::new();
        
            for line in exisiting_content.lines() {
                if line.starts_with('#') {
                    continue;
                }

                if line == "time,category,count" {
                    continue;
                }
        
                let cols: Vec<&str> = line.split(',').collect();
                let t:f64 = cols[0].parse()?;

                if existing_times.last() != Some(&t) {
                    existing_times.push(t);
                }

                let category = Category::from_key(cols[1])?;
                let count: usize = cols[2].parse()?;

                let t_idx = times.iter().position(|x| *x == t). expect("time in existing file should be one of the known timepoints");

                existing_counts.insert((t_idx, category), count);
            }

            if existing_times != times {
                    println!("Error: times of existing file don't match! Simulations can not be accumulated");
                    matches = false;
            }
        }

        matches 
    } else {
        false 
    };

    if !accumulate {
        existing_counts.clear();
        existing_num_sims = 0;
    }

    for (t_idx, map) in categories.iter().enumerate() {
            for (cat, count) in map {
                *existing_counts.entry((t_idx, cat.clone())).or_insert(0) += count;
            }
    }

    let total_num_sims = existing_num_sims + cli.num_sims;
   
    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
        writeln!(writer, "# t_split={}", co_time)?;
        writeln!(writer, "# t_ext={:?}", cli.simulation.t_ext.unwrap())?;
        writeln!(writer, "# t_end={}", post_time)?;
        let regions_str = regions.iter()
            .map(|(a, b)| format!("{}-{}", a, b))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(writer, "# regions={}", regions_str)?;
        writeln!(writer, "# num_sims={}", total_num_sims)?;
        write_categories(&mut writer, &times, &existing_counts)?;


    Ok(())
}




/// Cotranscriptional folding simulation 
fn run_cotimecourse<W, K>(
    moves: W,
    rmodel: K,
    sim_times: &[f64],
    num_sims: u64,
    times: &[f64],
) -> impl ParallelIterator<Item = Vec<DotBracketVec>>
where
    W: Walker + Clone + Send + Sync,
    K: RateModel + Clone + Send + Sync,
    SSA<W, K>: From<(W, K)>,
{
    let pb = ProgressBar::new(num_sims);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    (0..num_sims)
    .into_par_iter()
    .map_init(
        move || pb.clone(), // each thread gets a clone
        move |pb, _| {
           
            let mut structures: Vec<DotBracketVec> = Vec::with_capacity(times.len());

            let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
            let mut t_idx = 0;
            simulator.co_simulate(
                &mut rng(),
                sim_times,
                |t, tinc, _, w| {
                    while t_idx < times.len() && t + tinc >= times[t_idx] {
                        structures.push(w.current_structure());
                        t_idx += 1;
                    }
                    true
                },
            );

            pb.inc(1);
            structures
        },
    )

}

fn sort_by_timepoint(times: &Vec<f64>, structures: &Vec<Vec<DotBracketVec>>, num_sims: usize) -> Vec<Vec<DotBracketVec>> {

    let num_times = times.len();
    let mut sorted: Vec<Vec<DotBracketVec>> = vec![Vec::with_capacity(num_sims); num_times];

    for s in structures.iter() {
        for (t_idx, structure) in s.iter().enumerate() {
            sorted[t_idx].push(structure.clone());
        }
    }
    
    sorted
}

fn count_categories(structures: &[Vec<DotBracketVec>], regions: &Vec<(usize, usize)>,) -> Result<Vec<FxHashMap<Category, usize>>>{

    let mut result = Vec::with_capacity(structures.len());
   
    for timepoint_structures in structures.iter() {

        let mut counts: FxHashMap<Category, usize> = FxHashMap::default();

        for structure in timepoint_structures.iter() {

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

                *counts.entry(category).or_insert(0) += 1;
            
            }
        }

        result.push(counts); 

    }
    Ok(result)
}


fn write_categories( 
    writer: &mut impl Write,
    times: & [f64],
    counts: &FxHashMap<(usize, Category), usize>,
) -> Result<()> {
    writeln!(writer, "time,category,count")?;
    for (t_idx, t) in times.iter().enumerate() {
        for ((entry_t_idx, cat), count) in counts.iter() {
            if entry_t_idx == &t_idx {
                writeln!(writer, "{},{},{}", t, cat.to_key(), count)?;
            }
        }
    }
    Ok(())
}