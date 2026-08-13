use rustc_hash::FxHashMap;
use std::sync::Arc;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use std::fs::File;

use rayon::prelude::*;
use rand::rng;
use colored::*;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::Base;
use ff_energy::EnergyModel;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::shift_policy::*;
use ff_kinetics::SSA;
use ff_kinetics::RateModel; 

use fuzzyfold::input_parsers::read_cotr_input;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
use fuzzyfold::kinetics_parsers::TimelineParameters;

#[derive(Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimelineParameters,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

fn main() -> Result<()> {
    let mut cli = Cli::parse();
    let emodel = Arc::new(cli.energy.build_model());
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
    let num_ext = sequence.len() - structure.len();
    cli.simulation.validate(cli.kinetics.k0, num_ext)?;
    let _t_sep = cli.simulation.t_sep.expect("t-sep must exist after validation!");
    let pairings = PairTable::try_from(&structure)?;

    if let Some(h) = header {
        println!("{}", h.yellow());
    }
    println!("{}", sequence);
 
    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    let times = cli.simulation.get_output_times(num_ext)?;

    let (sim_times, _t_fin) = if num_ext > 0 {
        let t_ext = cli.simulation.t_ext.unwrap();
        let t_end = cli.simulation.t_end;
        let mut a = vec![t_ext; num_ext];
        a.push(t_end);
        (a, t_ext * (num_ext as f64) + t_end)
    } else { 
        (vec![cli.simulation.t_end], cli.simulation.t_end)
    };
    
    let timelines: Vec<Vec<DotBracketVec>> =
        match (rmodel.k3ws().is_some(), rmodel.k4ws().is_some()) {
            (false, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, NoShift))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
            (true, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
            (false, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, FourWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
            (true, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeAndFour))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
        };

    let mut master: Vec<FxHashMap<DotBracketVec, usize>> = vec![FxHashMap::default(); times.len()];    
    for timeline in timelines {        
        for (i, db) in timeline.into_iter().enumerate() {           
             *master[i].entry(db).or_insert(0) += 1;        
        }    
    } 

    let csv_path = cli.output.with_extension("csv");

    println!("{}", "Finished simulations!".red());

    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
    writeln!(writer, "# num_sims={}", cli.num_sims)?;
    writeln!(writer, "# t_sep={:?}", cli.simulation.t_sep.unwrap())?;
    writeln!(writer, "# t_ext={:?}",  cli.simulation.t_ext.unwrap())?;
    writeln!(writer, "# t_end={}", cli.simulation.t_end)?;
    writeln!(writer, "time,structure,count")?;

    for (t_idx, structures) in master.iter().enumerate() {
        let t = times[t_idx];
        for (structures, count) in structures.iter() {
            writeln!(writer, "{},{},{}", t, structure, count)?;
        }
    }
    println!("Wrote csv file: {}", format!("{}",csv_path.display()).green()); 

    Ok(())
}


fn run_timecourse<W, K>(
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
                let mut timeline: Vec<DotBracketVec> = Vec::with_capacity(times.len());
                let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
                let mut t_idx = 0;
                simulator.co_simulate(
                    &mut rng(),
                    sim_times,
                    |t, tinc, _, w| {
                        while t_idx < times.len() && t + tinc >= times[t_idx] {
                            timeline.push(w.current_structure());
                            t_idx += 1;
                        }
                        true
                    },
                );

                pb.inc(1);
                timeline
            },
        )
}