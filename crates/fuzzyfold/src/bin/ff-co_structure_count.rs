//! A command line tool to run multiple co-transcriptional folding simulations and identify the structures 
//! that appear most often throughout the simulation at defined timepoints. 
//!
//! This tool runs stochastic cotranscriptional folding simulations...
//! 
//! # Parameters
//! -'num-sims': number of simulations
//! -'t-ext': extensions time 
//! -'output': output file name
use rayon::prelude::*;
use rand::rng;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use ff_structure::PairTable;
use ff_structure::PairList;
use ff_structure::DotBracketVec;
use ff_energy::Base;
use ff_energy::EnergyModel;
use ff_kinetics::LoopStructure;
use ff_kinetics::LoopStructureSSA;

use fuzzyfold::input_parsers::read_fasta_like_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;


#[derive(Debug, Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time during transcription")]

pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Number of simulations
    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    /// Extension time 
    #[arg(long, default_value_t = 0.02)]
    t_ext: f64,

    /// Structures recorded per nucleotid given as frequency (1.0 = 1/nt, 0.5 = 0.5/nt, 2.0 = 2/nt)
    #[arg(long, default_value_t = 1.0)]
    frequency: f64,

    /// Output file name 
    #[arg(short, long, default_value = "accessibility.png")]
    output: String,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}


/// Vector to store structures in a single simulation run 
struct SimulationRun  {
    structures: Vec<PairList>, 
}

impl SimulationRun {

    // new simulation run
    pub fn new() -> Self {
        Self {
            structures: Vec::new()
        }
    }


    // Insert a new structure
    pub fn insert(&mut self, pairlist: PairList) {
        self.structures.push(pairlist);
    }

}

/// Build times_tl vector for timeline
fn build_timeline (
    frequency: f64,
    sequence: &[Base],
    t_ext: f64,
) -> Vec<f64> {

    if frequency <= 0.0 || sequence.len() == 0 {
        return Vec::new();
    }

    let mut times_tl = vec![0.0];
    let mut len = 1;

    let total_steps = sequence_len * freq;
    let total_time = sequence_len * t_ext;

    let step = total_time/ total_steps;

    for i in 1..total_steps {
        times_tl.push(i *)
    }

    if freq >= 1.0 {
        while len < sequence_len {
            let start = *times_tl.last().unwrap();
            let step = t_ext / freq as f64;
            for i in 1..= freq {
                times_tl.push(start + i as f64 * step);
            }
            len += 1;
        }
    } else {
        let count = sequence_len / freq;
        let mut c = 1;
        for i in 1..sequence_len {
            if c < count {
                c += 1;
            } else {
                let start = *times_tl.last().unwrap();
                times_tl.push(start + t_ext * count);
                c = 1;
            }
        }
    }
}

// Runs a single cotranscriptional folding simulation
fn run_simulation(
    sequence: &[Base],
    pairings: &PairTable,
    times: &Vec<f64>,
    frequency: f64,
    emodel: &impl EnergyModel,
    rmodel: &impl ff_kinetics::RateModel,
) -> Result<Vec<Vec<bool>>> {
    let loops = LoopStructure::try_from((&sequence[..], pairings, emodel)).unwrap();
    let mut simulator = LoopStructureSSA::from((loops, rmodel));

    let mut simulation_data SimulationRun = structures.new();
    let mut t_idx = 0; 


    simulator.co_simulate(
        &mut rng(),
        times.clone(),
        |t, tinc, _, ls| {
            
            if (tinc - times[0]) < 1e-10 && t_idx < sequence.len() { //record one structure per transcript length (Here structure before extension/ final structure)
                     
                    let structure = DotBracketVec::from(ls);
                    let pt = PairTable::try_from(&structure).unwrap();
                 
                    let row = t_idx + 1;

                    if row < sequence.len() {
                        for pos in 0..=row {
                            simulation_data[row][pos] = match pt[pos] {
                                None => true,      // Unpaired (accessible)
                                Some(_) => false,  // Paired (unaccessible)          
                            };
                        }  
                    } 
                t_idx += 1;
                }
        true
        }
    );
    Ok(simulation_data)
}





fn main() -> Result<()> {
    let cli = Cli::parse();

    let (_header, sequence, _structure) = read_fasta_like_input(&cli.input)?;
    let initial_structure = DotBracketVec::try_from(".")?; //always start at transcript length 1
    let pairings = PairTable::try_from(&initial_structure)?;
    let t_ext = cli.t_ext; //extension time (simulation time per nucleotide)

    // build times vector
    let mut times = Vec::new(); 
    for i in 0..(sequence.len()) {
        times.push((i+1) as f64 * t_ext);
    }

    // progress bar
    println!("Simulation progress:");
    let pb = ProgressBar::new(cli.num_sims as u64);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    let sequence_len = sequence.len();
    let num_sims = cli.num_sims;
    let mut matrix = AccessibilityMatrix::new(sequence.len()); //new matrix
     


    //run n simulations
    let simulations_results: Vec<_> = (0..cli.num_sims)
        .into_par_iter()
        .map(|_| {
            let emodel = cli.energy.build_model();
            let rmodel = cli.kinetics.build_model(emodel.temperature());
            let result = run_simulation(&sequence, &pairings, &times, &emodel, &rmodel);
            pb.inc(1);
            result
        })
        .collect();
    
    pb.finish_with_message("All simulations complete!");

    // add simulations to matrix
    for i in simulations_results.iter() {
        matrix.add(i.as_ref().unwrap());
    }

    // normalize matrix
    matrix.normalize(num_sims);

    //Plot
    let title = format!("RNA Accessibility - {} simulations", num_sims);
    plot_accessibility(&matrix, &cli.output, &title)?;

    Ok(())
    }

