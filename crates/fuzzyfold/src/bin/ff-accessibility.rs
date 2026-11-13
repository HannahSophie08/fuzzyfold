use rustc_hash::FxHashMap;
use std::sync::Arc;

use rayon::prelude::*;
use rand::rng;
use colored::*;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;


use ff_structure::{DotBracket, DotBracketVec, MultiPairTable, MultiStruct};
use ff_energy::Base;
use ff_energy::EnergyModel;
use fuzzyfold::kinetics::Elongation;
use fuzzyfold::kinetics::LoopNeighbors;
use fuzzyfold::kinetics::SSA;
use fuzzyfold::kinetics::RateModel; 

use fuzzyfold::input_parsers::read_cotr_input;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
use fuzzyfold::kinetics_parsers::TimeArguments;

#[derive(Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(short, long)]
    title: Option<String>,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimeArguments,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

fn main() -> Result<()> {
    let mut cli = Cli::parse();
    let emodel = Arc::new(cli.energy.build_model()?);
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) =
        if cli.simulation.t_ext.is_some() {
            read_cotr_input(&cli.input, is_rna)?
        } else {
            match read_eval_input(&cli.input, is_rna) {
                Ok(v) => v,
                Err(e) => return Err(anyhow::anyhow!("{e} (or use --t-ext?)")),
            }
        };

    if let Some(h) = header {
        println!("{}", h.yellow());
    }
    println!("{}", sequence);
 
    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    let num_ext = sequence.len() - structure.len();
    let pairings = MultiPairTable::try_from(&structure)?;

    let output = cli.simulation.output_times(num_ext)?;
    let schedule = cli.simulation.schedule(num_ext, &output)?;
    let times = output.times();

            
 
    let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, rmodel.move_types()))
        .map_err(|e| anyhow::anyhow!("failed to build the move set: {e}"))?;
    let timelines: Vec<Vec<DotBracketVec>> = run_timecourse(moves, rmodel, &schedule, cli.num_sims as u64, &times).collect();


    let mut master: Vec<FxHashMap<DotBracketVec, usize>> = vec![FxHashMap::default(); times.len()];    
    for timeline in timelines {        
        for (i, db) in timeline.into_iter().enumerate() {           
             *master[i].entry(db).or_insert(0) += 1;        
        }    
    } 

    for (i, counts) in master.iter().enumerate() {
        let mut entries: Vec<(&DotBracketVec, &usize)> = counts.iter().collect();
        entries.sort_by(|a, b| b.1.cmp(a.1)); // descending by count
        let seq_len = entries.first().map(|(db, _)| db.0.len()).unwrap_or(0);
        let mut unpaired = vec![0usize; seq_len];
        let mut total = 0;
        for (db, count) in entries {
            println!("  {:>6}  {}", count, db);
            for (pos, c) in db.to_string().chars().enumerate() {
                if c == '.' {
                    unpaired[pos] += count;
                }
            }
            total += count;
        }
        let normalized: Vec<f64> = unpaired.iter().map(|&x| x as f64 / total as f64).collect();
        println!("{:>8.5} {:?}", times[i], normalized);
    }
    /*
    let line: Vec<String> = paired
        .iter()
        .map(|&c| format!("{:.3}", c as f64 / num_sims as f64))
        .collect();
    println!("p_paired: [{}]", line.join(", "));
    */

    Ok(())
}

fn run_timecourse<W, K>(
    moves: W,
    rmodel: K,
    schedule: &[f64],
    num_sims: u64,
    times: &[f64],
) -> impl ParallelIterator<Item = Vec<DotBracketVec>>
where
    W: Elongation + Clone + Send + Sync,
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
    let done = std::sync::atomic::AtomicU64::new(0);

    (0..num_sims)
        .into_par_iter()
        .map_init(
            move || pb.clone(), // each thread gets a clone
            move |pb, _
            | {
                let mut timeline: Vec<DotBracketVec> = Vec::with_capacity(times.len());
                let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
                let mut t_idx = 0;
                simulator.co_simulate(
                    &mut rng(),
                    schedule,
                    |t, tinc, _, w| {
                        while t_idx < times.len() && t + tinc >= times[t_idx] {
                            timeline.push(w.current_structure());
                            t_idx += 1;
                        }
                        true
                    },
                );
                let n = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                if n % 10 == 0 || n == num_sims {
                    eprintln!("progress: {}/{}", n, num_sims);
                }

                pb.inc(1);
                timeline
            },
        )
}