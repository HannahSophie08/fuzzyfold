use std::fs;
use std::fs::File;
use std::path::Path;
use std::path::PathBuf;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::sync::Arc;

use colored::*;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use rayon::prelude::*;
use serde_json::to_string_pretty;

use ff_energy::EnergyModel;
use ff_energy::NucleotideVec;
use ff_kinetics::MacrostateRegistry;
use ff_kinetics::commit_and_delay::ExitMacrostateRegistry;
use ff_kinetics::commit_and_delay::CommitAndDelay;

use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;

#[derive(Debug, Parser)]
#[command(name = "ff-transitions")]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
pub struct Cli {
    #[arg(long, value_name = "FILE", num_args = 1.., required = true)]
    macrostates: Vec<PathBuf>,

    #[arg(short, long, default_value_t = 1_000)]
    num_sims: usize,

    #[arg(long, default_value_t = 0)]
    split_trajectories: usize,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
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
    let sequence = NucleotideVec::try_from_rna(&seq_line).unwrap();
    Ok(sequence)
}

fn io_err(msg: &str, src: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("{} in {}", msg, src))
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // --- Build simulator ---
    let emodel = Arc::new(cli.energy.build_model());
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    println!("Parameters:\n - {:?}\n - {:?}",
        cli.kinetics, cli.energy);

    let sequence = get_sequence(&cli.macrostates[0]).unwrap();
    let mut macrostates = MacrostateRegistry::from((sequence.clone(), emodel.clone()));
    let _ = macrostates.insert_files(&cli.macrostates, false);

    println!("{:>4} {:<10} {} {:>5} {:>8}",
        "ID",
        "Macrostate".cyan(), format!("{}", sequence).yellow(), "Size", "Energy");

    for (id, (_, m)) in macrostates.iter() {
        if m.name() == "Unassigned" {
            continue
        }
        println!("{:4} {:<10} {:<} {:>5} {:>8.2}",
            id, 
            m.name(),
            m.get_lowest_microstate().unwrap(),
            m.len(),
            m.ensemble_energy().unwrap());
    }
    let exitreg = Arc::new(ExitMacrostateRegistry::from((&macrostates, &rmodel)));
    
    // Output paths 
    let dat_path = cli.output.with_extension("dat");
    let crn_path = cli.output.with_extension("crn");

    // --- Load or initialize CommitAndDelay ---
    let mut cad = if Path::new(&dat_path).exists() {
        println!("{} {}", "Loading existing trajectory data from:".yellow(), dat_path.display());
        CommitAndDelay::load_json(&dat_path, Arc::clone(&exitreg))?
    } else {
        println!("{} {}", "Creating new trajectory database:".yellow(), dat_path.display());
        CommitAndDelay::from(Arc::clone(&exitreg))
    };

    // --- Simulations with progress bar ---
    println!("{}", "Simulating transitions:".green());
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

                for (id, _ms) in macrostates.iter() {
                    if id == 0 { continue; }
                    local_cad.simulate_from(id);
                }

                pb.inc(1);
                local_cad
            },
        )
        .collect();
    pb.finish_with_message("All simulations complete!");

    // merge results
    for local in results {
        cad.merge(local);
    }
    cad.recompute_marginals();

    println!("{}", "Number of successful transitions:".green());
    let rows: Vec<Vec<String>> = cad.trajectories()
        .rows()
        .into_iter()
        .map(|row| {
            row.iter()
                .map(|el| el.as_ref().map_or("0".into(), |ens| ens.len().to_string()))
                .collect()
        })
    .collect();

    // max width per column
    let ncols = rows.first().map_or(0, |r| r.len());
    let col_widths: Vec<usize> = (0..ncols)
        .map(|c| rows.iter().map(|r| r[c].len()).max().unwrap_or(0))
        .collect();

    let names: Vec<_> = macrostates.iter().map(|(_, (_, n))| n.name()).collect();
    let header = names.iter()
        .enumerate()
        .map(|(c, name)| {
            let name = if *name == "Unassigned" { "" } else { name };
            format!("{:>width$}", name, width = col_widths[c])
        })
    .collect::<Vec<_>>()
        .join(" ");
    println!("{:5} {header}", "");
    for (name, row) in names.iter().zip(&rows) {
        let name = if *name == "Unassigned" { "" } else { name };
        let line = row
            .iter()
            .enumerate()
            .map(|(c, cell)| format!("{:>width$}", cell, width = col_widths[c]))
            .collect::<Vec<_>>()
            .join(" ");
        println!("{:5} {line}", name);
    }

    std::fs::write(&crn_path, cad.to_crn_string(cli.split_trajectories))?;  
    println!("{} {}", "CRN written to".green(), crn_path.display());
  
    // --- Save checkpoint ---
    let serial = cad.to_serializable();
    let json = to_string_pretty(&serial)?;
    fs::write(&dat_path, json)?;
    println!("{} {}", "Trajectory data saved to".green(), dat_path.display());

    Ok(())
}


