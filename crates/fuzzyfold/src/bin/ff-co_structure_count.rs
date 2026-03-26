//! A command line tool to run multiple co-transcriptional folding simulations and 
//! identify the structures that appear most often throughout the simulation at defined 
//! time points. Writes all recorded structures and their count to an output file.
//! 
//! # Input: 
//! - '-- file_name/file_path'
//! - File containg sequence and initial structure
//! => initial structure defines start length
//! => if initial structure is not given or has full sequence length start length is 1
//!
//! # Output:
//! - CSV file with columns: 'timepoint_idx', 'time', 'structure', 'count'
//! - output file name can be defined in command line via '--output file_name'
//! - structures are written in the form of pair lists and are sorted by count per time point
//! 
//! 
//! # Parameters
//! -'num-sims': number of simulations
//! -'t-ext': extensions time (simulation time per transcript length)
//! -'t-end': simulation time at full sequence length
//! -'t-lin': equally distributed time points per transcript length where the structure is recorded
//! -'t-log': equally distributed time points for full sequence length where the structure is recorded
//! -'k0': rate parameter, set to 1e6 to give time parameters in seconds 
//! 
//! # Notes:
//! The structures are stored and compared as pair lists, which represent the base pairs of a structure,
//! but don't take the unpaired exterior loop into account. This allows to assign structures at different 
//! transcript lengths that only differ by the exterior loop to the same macrostate.

use rayon::prelude::*;
use rand::rng;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use std::collections::HashMap;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufWriter, Write};

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
use fuzzyfold::kinetics_parsers::TimelineParameters;


#[derive(Debug, Parser)]
#[command(version, about = "Record structure frequencies during cotranscriptional folding")]

pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Number of simulations
    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimelineParameters,

    /// Output file name 
    #[arg(short, long, default_value = "structure_frequencies.csv")]
    output: PathBuf,

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

    /// new simulation run
    pub fn new() -> Self {
        Self {
            structures: Vec::new()
        }
    }

    /// Insert a new structure
    pub fn insert(&mut self, pairlist: PairList) {
        self.structures.push(pairlist);
    }

}
/// Vector that combines multiple simulation runs, compares the structures and counts
/// how often each structure occurs across the simulations per time point 
struct CombinedSimulations {
    pub counts: Vec<HashMap<PairList, usize>>,
}

impl CombinedSimulations {

    /// generate new map
    pub fn new(num_timepoints: usize) -> Self {
        Self {
            counts: vec![HashMap::new(); num_timepoints],  
        }
    }

    /// add a simulation
    pub fn add(&mut self, run: SimulationRun) {

        // check that length of the structures selected matches number of timepoints 
        assert_eq!(run.structures.len(), self.counts.len(), "Number of timepoints doesn't match");
        
        for t_idx in 0..run.structures.len() {
            let structure = run.structures[t_idx].clone();
            let timepoint = &mut self.counts[t_idx]; //timepoint that the structure was collected
            
            if let Some(count) = timepoint.get_mut(&structure) { //if structure already exists, add 1 to count
                *count += 1;
            } else { // new structure 
                timepoint.insert(structure, 1);
            }
        }
    }

    /// sort by count for each timepoint 
    pub fn sort(&self, t_idx:usize) -> Vec<(&PairList, usize)> {
        
        let mut entries: Vec<(&PairList, usize)> = self.counts[t_idx]  //convert into vector
            .iter()
            .map(|(structure, count) | (structure, *count))
            .collect();

        entries.sort_by(|a, b| b.1.cmp(&a.1)); // sort in descending order
        return entries;
    }

 }


/// vector of time points at which structures are recorded during simulation
fn build_timeline (
    sequence_len: usize,
    start_len: usize,
    t_lin: usize, // timepoints per nucleotide
    t_ext: f64, // extension time (simulation time per nucleotide)
    t_log: usize, // timepoints for posttranscriptional folding 
    t_end: f64, // posttranscriptional simulation time 
) -> Vec<f64> {

    
    let mut times_tl = vec![0.0];

    // Build times_tl vector for timeline
    // Linear segments: append `t_lin` evenly spaced points
    let mut len = start_len;

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

    return times_tl;
}



///  Build vector of accumulative simulation times as input for the simulation
fn simulation_times(
    start_len: usize,
    sequence_len: usize,
    t_ext: f64,
    t_end: f64,
) -> Vec<f64>  {
    
    let mut times = Vec::new(); 

    times.push(t_ext);
    for _ in (start_len + 1)..(sequence_len) {
        times.push(times.last().unwrap() + t_ext);
    }
    times.push(times.last().unwrap() + t_end);

    return times;
}


/// Run a single cotranscriptional folding simulation, recording structure at 
/// each time point in 'timeline'
fn run_simulation(
    sequence: &[Base],
    pairings: &PairTable,
    times: &Vec<f64>,
    timeline: &[f64],
    emodel: &impl EnergyModel,
    rmodel: &impl ff_kinetics::RateModel,
) -> Result<SimulationRun> {

    let loops = LoopStructure::try_from((&sequence[..], pairings, emodel)).unwrap();
    let mut simulator = LoopStructureSSA::from((loops, rmodel));

    let num_timepoints  = timeline.len();
    let mut run = SimulationRun::new();
    let mut t_idx = 0; 

    simulator.co_simulate(
        &mut rng(),
        times.clone(),
        |t, tinc, _, ls| {
            
            while t_idx < timeline.len() && t + tinc >= timeline[t_idx] {
                
                let structure = DotBracketVec::from(ls);
                let pt = PairTable::try_from(&structure).unwrap();
                run.insert(PairList::from(&pt));
            
                t_idx += 1;
                }
        true
        }
    );

    assert_eq!(run.structures.len(), num_timepoints, "Simulation produced {} structures but expected {}",
        run.structures.len(), num_timepoints);
   
    Ok(run)
}


/// writes the by count sorted results into an csv output file
fn write_output_file(
    path: &PathBuf,
    combined: &CombinedSimulations,
    timeline: &[f64],
    ) -> Result<()> {
    
    // open a new file 
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    // column names
    writeln!(w, "timepoint_idx, time, structure, count")?;
    
    for (t_idx, &time) in timeline.iter().enumerate() {
        for (pairlist, count) in combined.sort(t_idx) {
            writeln!(w, "{}, {:.6}, {}, {:.6}", t_idx, time, pairlist, count)?;
        }
    }
    Ok(()) 
}


fn main() -> Result<()> {
    let cli = Cli::parse();

    let (_header, sequence, structure) = read_fasta_like_input(&cli.input)?;
    let mut initial_structure = structure; 
    let mut start_len = initial_structure.len(); //transcript length from which the simulation starts
    let sequence_len = sequence.len();

    // if given structure is full sequence length, start from length one 
    if start_len == sequence_len {
        start_len = 1;
        initial_structure = DotBracketVec::try_from(".")?;
    }

    let pairings = PairTable::try_from(&initial_structure)?;
    

    // Simulation parameters
    let t_ext = cli.simulation.t_ext; //extension time (simulation time per nucleotide)
    let t_end = cli.simulation.t_end; //simulation time at full sequence length
    let t_lin = cli.simulation.t_lin; //time points recorded per transcript length
    let t_log = cli.simulation.t_log; //time points recorded at full sequence length 

    
    let timeline = build_timeline(sequence_len, start_len, t_lin, t_ext, t_log, t_end); //build timeline of timepoints
    let times = simulation_times(start_len, sequence_len, t_ext, t_end); // build vector of simulation times 

    
    // progress bar
    println!("Simulation progress:");
    let pb = ProgressBar::new(cli.num_sims as u64);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );


    //run n simulations
    let simulations_results: Vec<_> = (0..cli.num_sims)
        .into_par_iter()
        .map(|_| {
            let emodel = cli.energy.build_model();
            let rmodel = cli.kinetics.build_model(emodel.temperature());
            let result = run_simulation(&sequence, &pairings, &times, &timeline, &emodel, &rmodel);
            pb.inc(1);
            result
        })
        .collect();
    
    pb.finish_with_message("All simulations complete!");


    // Combine simulation results
    let mut combined = CombinedSimulations::new(timeline.len());

    for result in simulations_results {
        match result {
            Ok(run) => combined.add(run),
            Err(_) => println!("Error: Simulation failed!"),
        }
    }

    // sort results and write them to csv file
    let _ = write_output_file(&cli.output, &combined, &timeline);


    Ok(())
    }

