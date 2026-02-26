//! Stochastic cotranscriptional folding simulator 
//! 
//! This binary performs stochastic cotranscriptional folding simulations, running multiple simulations in parallel
//! and merging them into a master timeline. This timeline is used to produce an occupancy plot.
//! The plot shows the occupancy of predefined macrostates in average over time.
//! 
//! Input & start state
//!-  If no structure is provided, in the input fasta file, the simulation starts by default 
//!   at transcript length 1 with ".".
//! - If an initial dot-bracket structure is provided, its length defines the starting transcript lenth.
//!   The structure is used as the initial structure.
//! - If the provided structure is full length, it is ignored and the simzlation starts at length 1.  
//! 
//! Parameters
//! - t-ext: extension time (simulation time per transcript length)
//! - t-end: time of postranscriptional simulation
//! - t-lin: recorded time steps per transcript length
//! - t-log: time points recorded for the posttranscriptional simulation on a logarithmic timescale 
//! - p-pos: pausing sites positions
//! - t-pos: duration of pausing sites corresponding to `--p-pos`` 
//! - num-sims: number of simulations performed
//! - k0: kinetic rate constant 

use std::fs;
use std::sync::Arc;
use std::path::Path;
use std::path::PathBuf;
use rayon::prelude::*;
use rand::rng;
use clap::Parser;
use anyhow::Result;
use colored::*;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use serde_json::to_string_pretty;

use ff_structure::PairTable;
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;
use ff_kinetics::LoopStructure;
use ff_kinetics::LoopStructureSSA;
use ff_kinetics::timeline_pairlists::Timeline;
use ff_kinetics::timeline_plotting_pairlists::plot_occupancy_over_time;
use ff_kinetics::MacrostateRegistryPL;
//use ff_kinetics::timeline::Timeline;
//use ff_kinetics::timeline_plotting::plot_occupancy_over_time;
//use ff_kinetics::MacrostateRegistry;

use fuzzyfold::input_parsers::read_fasta_like_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
use fuzzyfold::kinetics_parsers::TimelineParameters;

#[derive(Debug, Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(long, value_name = "FILE", num_args = 1.., required = false)]
    macrostates: Vec<PathBuf>,

    /// Backup/Store timeline in this file.
    #[arg(long, value_name = "FILE")]
    timeline: Option<PathBuf>,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimelineParameters,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,

    /// Pausing sites 
    #[arg(long, value_delimiter = ',')]
    t_p: Option<Vec<f64>>,

    /// Pausing sites positions  
    #[arg(long, value_delimiter = ',')]
    p_pos: Option<Vec<usize>>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.simulation.validate()?;

    // --- Build simulator ---
    let emodel = cli.energy.build_model();
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let (header, sequence, mut structure) = read_fasta_like_input(&cli.input)?;
    if structure.len() >= sequence.len() { // if the given structure is full length, start from transcript length one 
        structure = DotBracketVec::try_from(".")?;
        println!("Input structure is full-length");
    }
    let pairings = PairTable::try_from(&structure)?;

    let name = if let Some(h) = header {
        println!("{}", h.yellow());
        h.strip_prefix('>')
            .and_then(|s| s.split_whitespace().next())
            .unwrap_or("anonymous")
            .to_string()
    } else {
        println!("{}", ">anonymous".yellow());
        "anonymous".to_string()
    };
    println!("{}", sequence);

    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    // --Output times--

    let sequence_len = sequence.len();

    let t_lin = cli.simulation.t_lin; //timepoints per nucleotide
    let t_ext = cli.simulation.t_ext; //extension time (simulation time per nucleotide)
    let t_log = cli.simulation.t_log; //timepoints for posttranscriptional folding 
    let t_end = cli.simulation.t_end; // posttranscriptional simulation time
    
    let mut times_tl = vec![0.0];

    // Build times_tl vector for timeline
    // Linear segments: append `t_lin` evenly spaced points
    let mut len = 1;

    while len < sequence_len {
        let start = *times_tl.last().unwrap();
        let step = t_ext / t_lin as f64;
        for i in 1..= t_lin {
            times_tl.push(start + i as f64 * step);
        }
        len += 1;
    }

    // Logarithmic tail
    let start = *times_tl.last().unwrap();
    let log_start = start.ln();
    let log_end = (start + t_end).ln();
    for i in 1..t_log {
        let frac = i as f64 / t_log as f64;
        let value = (log_start + frac * (log_end - log_start)).exp();
        times_tl.push(value);
    }
    times_tl.push(start + t_end);

    let mut registry = MacrostateRegistryPL::from((&sequence, &emodel));
    registry.insert_files(&cli.macrostates)?;

    println!("Macrostates:\n{}", registry.iter()
        .map(|(_, m)| format!(" - {} {:6.2}", m.name(), m.ensemble_energy().unwrap_or(0.0)))
        .collect::<Vec<_>>().join("\n"));

    let shared_registry = Arc::new(registry);
    
    // Create or load master timeline
    // If timeline.json exists, reload instead of starting empty
    let mut master = if let Some(path) = &cli.timeline {
        if Path::new(path).exists() {
            println!("Loading existing timeline from: {}", path.display());
            Timeline::from_file(path, &times_tl, Arc::clone(&shared_registry))?
        } else {
            println!("A new timeline file will be created: {}", 
                path.display());
            Timeline::new(&times_tl, Arc::clone(&shared_registry))
        }
    } else {
        Timeline::new(&times_tl, Arc::clone(&shared_registry))
    };

   
    //---build times vector--- 
    let mut times: Vec<f64> = Vec::new();

    let start_len = structure.len(); //starting transcript length
    let mut idx = 0;

    // pausing sites
    if let Some(pause_times) = &cli.t_p {
       if let Some(pause_positions) = &cli.p_pos {
            if pause_times.len() != pause_positions.len() { //Check that pausing positions and times match
                panic!("p_pos and t_pos don't have the same length")
            }
            let mut pos = start_len;
            times.push(t_ext);
            pos += 1;
            for &p in pause_positions {
                while pos < p {
                    times.push(times.last().unwrap() + t_ext);
                    pos += 1;
                }
                if pos == p {
                    times.push(times.last().unwrap() + pause_times[idx]);
                    idx += 1;
                    pos += 1;
                }
            }
            while pos < (sequence.len()) {
                times.push(times.last().unwrap() + t_ext);
                pos += 1;
            }
            // Add final postranscriptional folding time 
            times.push(times.last().unwrap() + t_end);
            
       }
    } else { //no pausing sites

        times.push(t_ext);
        for _ in (start_len + 1)..(sequence.len()) {
            times.push(times.last().unwrap() + t_ext);
        }
        // Add final posttranscriptional folding time 
        times.push(times.last().unwrap() + t_end);
    }

    // ---Progress bar---
    println!("Simulation progress:");
    let pb = ProgressBar::new(cli.num_sims as u64);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    // --- Run simulations in parallel ---
    // Each simulation produces its own timeline
    let timelines: Vec<_> = (0..cli.num_sims)
        .into_par_iter()
        .map_init(
            || pb.clone(), // each thread gets a clone
            |pb, _| {
                let registry = Arc::clone(&shared_registry);
                let mut timeline = Timeline::new(&times_tl, registry);
                // Build loop structure and simulator 
                let loops = LoopStructure::try_from((&sequence[..], &pairings, &emodel)).unwrap();
                let mut simulator = LoopStructureSSA::from((loops, &rmodel));
                
                let mut t_idx = 0;
                // Simulation 
                simulator.co_simulate(
                    &mut rng(),
                    times.clone(),
                    |t, tinc, _, ls| {
                        while t_idx < times_tl.len() && t + tinc >= times_tl[t_idx] {
                            let structure = DotBracketVec::from(ls);
                            timeline.assign_structure(t_idx, &structure);
                            t_idx += 1;
                        }
                        true
                    },
                );

                pb.inc(1);
                timeline
            },
        ).collect();
    pb.finish_with_message("All simulations complete!");

    // Merge timelines to master timeline 
    for timeline in timelines {
        master.merge(timeline);
    }

    let cotrans_time = times[times.len() - 2]; 

    println!("Final Timeline:\n{}", master);
    // Generate occupancy plot
    plot_occupancy_over_time(&master, &format!("ff_{}.svg", name), cotrans_time, *times.last().unwrap());

    if let Some(path) = cli.timeline {
        let serial = master.to_serializable();
        let json = to_string_pretty(&serial)?;
        fs::write(path, json)?;
    }

    Ok(())
}



