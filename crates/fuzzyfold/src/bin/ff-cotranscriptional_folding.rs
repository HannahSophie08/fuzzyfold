use clap::Args;
use clap::Parser;
use colored::*;
use anyhow::Result;

use rand::rng; //random number generator 
use ff_structure::PairTable; //base pairings
use ff_energy::EnergyModel; 
use ff_kinetics::LoopStructure;
use ff_kinetics::LoopStructureSSA;
use ff_kinetics::Metropolis;
use ff_structure::DotBracketVec;
use ff_structure::DotBracket;

use fuzzyfold::input_parsers::read_fasta_like_input; //reads input
use fuzzyfold::energy_parsers::EnergyModelArguments;
//TODO: support seeded rng.

#[derive(Debug, Args)]
pub struct RateModelParams {
    /// Metropolis rate constant (must be > 0).
    #[arg(long, default_value_t = 1.0)] //default: 1.0
    pub k0: f64, 
}

#[derive(Debug, Parser)]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]

pub struct Cli { //Command line interface 
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Transcription rate
    #[arg(long, default_value_t = 0.02)]
    trans_rate: f64,

    #[command(flatten, next_help_heading = "Kinetic model parameters")] //include fields of RateModelParams directly as CLI flags
    kinetics: RateModelParams,

    #[command(flatten, next_help_heading = "Energy model parameters")] //include fields of EnergyModelArguments directly as CLI flags
    energy: EnergyModelArguments,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // --- Build simulator ---
    let emodel = cli.energy.build_model(); 
    let rmodel = Metropolis::new(emodel.temperature(), cli.kinetics.k0); //Create Metropolis kinetic rate model

    let (header, sequence, _structure) = read_fasta_like_input(&cli.input)?; //read sequence and structure

    // start with len=1
    let mut current_len = 1;

    let mut t = 0.0; //simulation time 

    if let Some(h) = header { //print header 
         println!("{}", h.yellow())
        }

    println!("{} {:>8} {:>14} {:>14} {:>15}",
            sequence,
            "energy".green(),
            "arrival time".cyan(),
            "waiting time".cyan(),
            "mean-waiting".cyan(),
        );
    
    let mut current_struct = DotBracketVec(vec![DotBracket::Unpaired]); //initialize current structure with length 1

    while current_len < sequence.len() { //repeat until end of transcription

        let mut current_seq = &sequence[..current_len];


        let pairings = PairTable::try_from(&current_struct)?; //build PairTable
        let loops = LoopStructure::try_from((&current_seq[..], &pairings, &emodel)).unwrap(); //build loop structure 

        // --- Check if code panics, because no possible transitions for current seequence, skip simulation in this case ---
         
        let add_pair = loops // can a base pair be added?
            .get_add_neighbors_per_loop()
            .iter()
            .any(|(_, add_neighbors) | !add_neighbors.is_empty());

        let del_pair = !loops.get_del_neighbors().is_empty(); //can a base pair be deleted?

        if !add_pair && !del_pair { // no possible transitions => continue without simulation and extend structure and sequence
            t = t + cli.trans_rate;
            current_len += 1;
            current_struct.0.push(DotBracket::Unpaired);
            continue;
        }

        let mut simulator = LoopStructureSSA::from((loops, &rmodel)); //build simulator from loop structure and rate model
        let t_next = t + cli.trans_rate; //current t_max

        simulator.simulate(
            &mut rng(), //random number 
            t_next,
            |t, tinc, flux, ls| {
                println!("{} {:8.2} {:14.8e} {:14.8e} {:15.8e}", 
                    ls,
                    ls.energy() as f64 / 100.,
                    t,
                    tinc,
                    1.0 / flux,
                );
                true
            },
        );

    t = t_next; //update t

    // ---Append sequence with next nucleotide---
    current_len += 1; 
    current_seq = &sequence[..current_len];
    current_struct.0.push(DotBracket::Unpaired);

    }
    
    Ok(())
    
}