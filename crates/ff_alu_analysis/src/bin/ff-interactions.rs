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
use plotters::prelude::*;
use plotters::style::Palette;
use plotters::style::IntoFont;
use plotters::prelude::IntoLogRange;
use plotters::style::Color;

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

    let svg_path = cli.output.with_extension("svg");
    let csv_path = cli.output.with_extension("csv");

    // --- Check input ---
    if structure.len() < sequence.len() && cli.simulation.t_ext.is_none() {
        panic!("Error:'t-ext' (extension time) missing for cotranscriptional simulation!"); 
    }

    if structure.len() == sequence.len() && cli.simulation.t_ext.is_some() {
        panic!("Error: 't-ext' (extension time) given in combination full length structure. For cotranscriptional simulation give either 
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

    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
    writeln!(writer, "# t_split={}", co_time)?;
    writeln!(writer, "# t_end={}", post_time)?;
    write_categories(&mut writer, &times, &categories)?;
   
    
        
    println!("Plotted svg file: {}", svg_path.display());


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

fn count_categories(structures: &[Vec<DotBracketVec>], regions: &Vec<(usize, usize)>,) -> Result<Vec<FxHashMap<Category, f64>>>{

    let mut result = Vec::with_capacity(structures.len());
   
    for timepoint_structures in structures.iter() {

        let mut counts: FxHashMap<Category, usize> = FxHashMap::default();
        let mut total_pairs = 0; // total pairs in region 


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
                total_pairs += 1;
            
            }
        }

        let percentages = counts.into_iter()
            .map(|(cat, count)| (cat, count as f64 / total_pairs as f64 * 100.0))
            .collect();

        result.push(percentages); 
    }
    Ok(result)
}


fn write_categories( 
    writer: &mut impl Write,
    times: & [f64],
    categories: &[FxHashMap<Category, f64>],
) -> Result<()> {
    writeln!(writer, "time,category,value")?;
    for (t, map) in times.iter().zip(categories.iter()) {
        for (cat, val) in map {
            writeln!(writer, "{},{},{}", t, cat.to_key(), val)?;
        }
    }
    Ok(())
}