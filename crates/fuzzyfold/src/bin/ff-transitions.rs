use clap::Parser;
use colored::*;
use anyhow::Result;
use ff_energy::NucleotideVec;
use std::sync::Arc;
use std::path::PathBuf;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io;



use std::fs;
//use std::io::{self, BufRead, BufReader};
use std::path::{Path};

use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde_json::to_string_pretty;


use ff_energy::EnergyModel;
use ff_kinetics::Metropolis;
use ff_kinetics::MacrostateRegistry;

use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelParams;
use ff_kinetics::commit_and_delay::ExitMacrostateRegistry;
use ff_kinetics::commit_and_delay::CommitAndDelay;

#[derive(Debug, Parser)]
#[command(name = "ff-transitions")]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
pub struct Cli {
    #[arg(long, value_name = "FILE", num_args = 1.., required = true)]
    macrostates: Vec<PathBuf>,

    /// Where to store / reload the simulation database
    #[arg(long, value_name = "FILE")]
    database: Option<PathBuf>,
    
    #[arg(short, long, default_value_t = 1_000)]
    num_sims: usize,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelParams,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,
}

fn get_sequence(msfile: &PathBuf) -> io::Result<NucleotideVec> {
    let file = File::open(msfile).expect("Failed to open first macrostate file");
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let header_line = lines
        .next()
        .ok_or_else(|| io_err("Missing header line", &msfile.display().to_string()))??
        .trim()
        .to_string();

    let _ = header_line
        .strip_prefix('>')
        .ok_or_else(|| io_err("First line must start with '>'", &msfile.display().to_string()))?
        .trim()
        .to_string();

    let seq_line = lines
        .next()
        .ok_or_else(|| io_err("Missing sequence line", &msfile.display().to_string()))??
        .trim()
        .to_string();
    let sequence = NucleotideVec::from_lossy(&seq_line);
    Ok(sequence)
}

fn io_err(msg: &str, src: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("{} in {}", msg, src))
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // --- Build simulator ---
    let emodel = cli.energy.build_model();
    let rmodel = Metropolis::new(emodel.temperature(), cli.kinetics.k0);

    println!("Parameters:\n - {:?}\n - {:?}",
        cli.kinetics, cli.energy);

    let sequence = get_sequence(&cli.macrostates[0]).unwrap();
    let mut registry = MacrostateRegistry::from((&sequence, &emodel));
    let _ = registry.insert_files(&cli.macrostates);

    println!("Macrostates:\n{}", registry.iter()
        .map(|(_, m)| format!(" - {} {:6.2}", m.name(), m.ensemble_energy().unwrap_or(0.0)))
        .collect::<Vec<_>>().join("\n"));

    let exitreg = ExitMacrostateRegistry::from((&registry, &rmodel));
        let exitreg = Arc::new(exitreg);
    
    // --- Load or initialize CommitAndDelay ---
    let mut cad = if let Some(path) = &cli.database {
        if Path::new(path).exists() {
            println!("{} {}", "Loading existing database from:".yellow(), path.display());
            CommitAndDelay::load_json(path, Arc::clone(&exitreg))?
        } else {
            println!("{} {}", "Creating new database:".yellow(), path.display());
            CommitAndDelay::from(Arc::clone(&exitreg))
        }
    } else {
        CommitAndDelay::from(Arc::clone(&exitreg))
    };

    // --- Simulations with progress bar ---
    println!("{}\n", "Simulating transitions:".green());
    let pb = ProgressBar::new(cli.num_sims as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );
    let results: Vec<_> = (0..cli.num_sims)
        .into_par_iter()
        .map_init(
            || pb.clone(),
            |pb, _| {
                let mut local_cad = CommitAndDelay::from(Arc::clone(&exitreg));

                for (id, _ms) in registry.iter() {
                    if id == 0 { continue; }
                    local_cad.simulate_from(id);
                }

                pb.inc(1);
                local_cad
            },
        )
        .collect();

    // merge results
    for local in results {
        cad.merge(local);
    }
    pb.finish_with_message("All simulations complete!");

    cad.gather_data();

    println!("{}", "Aggregated transition matrix:".green());
    for row in cad.trajectories().rows() {
        let line = row
            .iter()
            .map(|el| el.as_ref().map_or("0".into(), |ens| ens.len().to_string()))
            .collect::<Vec<_>>()
            .join(" ");
        println!("{line}");
    }
  
    // --- Save checkpoint ---
    if let Some(path) = cli.database {
        let serial = cad.to_serializable();
        let json = to_string_pretty(&serial)?;
        fs::write(&path, json)?;
        println!("{} {}", "Database saved to".green(), path.display());
    }

    Ok(())
}


