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

use ff_structure::{DotBracketVec, MultiPairTable, MultiStruct};
use ff_energy::Base;
use ff_energy::EnergyModel;
use ff_energy::NucleotideVec;
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

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimeArguments,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

pub struct Pattern {
    name: &'static str,
    paired: bool,
    opposite: Option<Base>,
    duplex_3: usize,
    duplex_5: usize,
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
    let num_ext = sequence.len() - structure.len();
    let pairings = MultiPairTable::try_from(&structure)?;

    let output = cli.simulation.output_times(num_ext)?;
    let schedule = cli.simulation.schedule(num_ext, &output)?;
    let _times = output.times();

    if let Some(h) = header {
        println!("{}", h.yellow());
    }
    println!("{}", sequence);
 
    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);


    let ac_3: Pattern = Pattern {
        name: "ac_3",
        paired: false,
        opposite: Some(Base::C),
        duplex_3: 1,
        duplex_5: 1,
    };

    let ac_5: Pattern = Pattern {
        name: "ac_5",
        paired: false,
        opposite: Some(Base::C),
        duplex_3: 2,
        duplex_5: 2,
    };

    let ac_7: Pattern = Pattern {
        name: "ac_7",
        paired: false,
        opposite: Some(Base::C),
        duplex_3: 3,
        duplex_5: 3,
    };
    
    let au_3: Pattern = Pattern {
        name: "au_3",
        paired: true,
        opposite: Some(Base::U),
        duplex_3: 1,
        duplex_5: 1,
    };

    let au_5: Pattern = Pattern {
        name: "au_5",
        paired: true,
        opposite: Some(Base::U),
        duplex_3: 2,
        duplex_5: 2,
    };

    let au_7: Pattern = Pattern {
        name: "au_7",
        paired: true,
        opposite: Some(Base::U),
        duplex_3: 3,
        duplex_5: 3,
    };

    let ag_3: Pattern = Pattern {
        name: "ag_3",
        paired: false,
        opposite: Some(Base::G),
        duplex_3: 1,
        duplex_5: 1,
    };

    let ag_5: Pattern = Pattern {
        name: "ag_5",
        paired: false,
        opposite: Some(Base::G),
        duplex_3: 2,
        duplex_5: 2,
    };

    let ag_7: Pattern = Pattern {
        name: "ag_7",
        paired: false,
        opposite: Some(Base::G),
        duplex_3: 3,
        duplex_5: 3,
    };

    let bulge_3: Pattern = Pattern {
        name: "bulge_3",
        paired: false,
        opposite: None,
        duplex_3: 1,
        duplex_5: 1,
    };

    let bulge_5: Pattern = Pattern {
        name: "bulge_5",
        paired: false,
        opposite: None,
        duplex_3: 2,
        duplex_5: 2,
    };

    let bulge_7: Pattern = Pattern {
        name: "bulge_7",
        paired: false,
        opposite: None,
        duplex_3: 3,
        duplex_5: 3,
    };

    let patterns: Vec<Pattern> = vec![ac_3, ac_5, ac_7, au_3, au_5, au_7, ag_3, ag_5, ag_7, bulge_3, bulge_5, bulge_7];

    let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, rmodel.move_types()))
        .map_err(|e| anyhow::anyhow!("failed to build the move set: {e}"))?;
    let editable_times: Vec<Vec<Vec<f64>>> = run_timecourse(moves, rmodel, &schedule, cli.num_sims as u64, &sequence, &patterns).collect();

    println!("{}", "Finished simulations!".red());

    let csv_path = cli.output.with_extension("csv");
    let seq_str = sequence.to_string();

    for (k, pattern) in patterns.iter().enumerate() {
        let path = cli.output.with_file_name(format!("{}_{}.csv", cli.output.file_stem().unwrap().to_string_lossy(), pattern.name));
        let mut writer = BufWriter::new(File::create(path)?);
        writeln!(writer, "# duplex_5={}", pattern.duplex_5)?;
        writeln!(writer, "# duplex_3={}", pattern.duplex_3)?;
        writeln!(writer, "# opposite={:?}", pattern.opposite)?;
        writeln!(writer, "# paired={}", pattern.paired)?;
        writeln!(writer, "# t-ext={:?}", cli.simulation.t_ext)?;
        writeln!(writer, "# t_end={}", cli.simulation.t_end)?;
        writeln!(writer, "# num_sims={}", cli.num_sims)?;
        
        write!(writer, "position\tnucleotide")?;
        for t in 1..=editable_times.len() {
            write!(writer, "\tT{}", t)?;
        }
        writeln!(writer)?;

        for (i, nt) in seq_str.chars().enumerate() {
            write!(writer, "{}\t{}", i + 1, nt)?;
            let is_a = sequence[i] == Base::A;
            for times in &editable_times {
                if is_a {
                    write!(writer, "\t{}", times[k][i])?;
                } else {
                    write!(writer, "\tNA")?;
                }
            }
            writeln!(writer)?;
        }
        writer.flush()?;
    }

    Ok(())
}

fn run_timecourse<W, K>(
    moves: W,
    rmodel: K,
    sim_times: &[f64],
    num_sims: u64,
    sequence: &NucleotideVec,
    patterns: &Vec<Pattern>,
) -> impl ParallelIterator<Item = Vec<Vec<f64>>>
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

    let adenosines: Vec<usize> = sequence.iter().enumerate()
            .filter(|(_, base)| **base == Base::A)
            .map(|(i, _)| i)
            .collect();

    (0..num_sims)
        .into_par_iter()
        .map_init(
            move || pb.clone(), // each thread gets a clone
            move |pb, _
            | { 
                let mut editable_time: Vec<Vec<f64>> = vec![vec![0.0; sequence.len()]; patterns.len()];
                let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
                simulator.co_simulate(
                    &mut rng(),
                    sim_times,
                    |t, tinc, _, w| {
                        let current_struct = w.current_structure();
                        let mpt = MultiPairTable::try_from(&current_struct).expect("Failed to build pair table from current structure");
                        for (n, pattern) in patterns.iter().enumerate() {
                            let editable_positions = position_check(sequence, &current_struct, &mpt, &adenosines, &pattern);
                            for (i, pos) in editable_positions.iter().enumerate() {
                                if *pos {
                                    editable_time[n][adenosines[i]] += tinc;
                                }
                            }
                        } 
                        true
                    },
                );
                pb.inc(1);
                editable_time
            },
        )
}


// checks whether position is unpaired and across from a 'C', and whether the 
// 5' and 3' neighbors are paired
fn position_check (sequence: &NucleotideVec, structure: &DotBracketVec, mpt: &MultiPairTable, positions: &Vec<usize>, pattern: &Pattern) -> Vec<bool> {

    let mut result = Vec::new();

    'outer_pos: for position in positions.iter() {

        if *position < pattern.duplex_5 || *position + pattern.duplex_3 >= structure.len() {
            result.push(false);
            continue;
        }
        
        if !pattern.paired {
            if mpt.get(*position) != MultiStruct::Unpaired {
                result.push(false);
                continue;
            }
        } else {
            if mpt.get(*position) == MultiStruct::Unpaired {
                result.push(false);
                continue;
            }
        }

        let mut outer = Vec::new();
        for i in 1..=pattern.duplex_5 {
            outer.push((*position - i) as usize);
        }

        let mut inner = Vec::new();
        for i in 1..=pattern.duplex_3 {
            inner.push((*position + i) as usize);
        }

        let mut outer_partner = Vec::new();
        for o in outer {
             match mpt.get(o) {
                MultiStruct::Paired(p) => outer_partner.push(p),
                _ => {
                    result.push(false);
                    continue 'outer_pos;
                }
            }
        } 

        let mut inner_partner = Vec::new();
        for i in inner {
            match mpt.get(i) {
                MultiStruct::Paired(p) => inner_partner.push(p),
                _ => {
                    result.push(false);
                    continue 'outer_pos;
                }
            }
        }
        
        for idx in 0..outer_partner.len() {
            if idx + 1 < outer_partner.len() {
                if outer_partner[idx] +  1 != outer_partner[idx + 1] {
                    result.push(false );
                    continue 'outer_pos; 
                }
            }
        }

        for idx in 0..inner_partner.len() {
            if idx + 1 < inner_partner.len() {
                if inner_partner[idx+1] + 1 != inner_partner[idx] {
                    result.push(false );
                    continue 'outer_pos; 
                }
            }
        }

        let duplex_len = pattern.duplex_5 + pattern.duplex_3;

        if !pattern.opposite.is_none() {
            if outer_partner[pattern.duplex_5 -1] > inner_partner[pattern.duplex_3 - 1] {
            
                if outer_partner[pattern.duplex_5 -1] - inner_partner[pattern.duplex_3 - 1] == duplex_len && sequence[outer_partner[0] - 1] == pattern.opposite.unwrap() {
                    result.push(true);
                    continue;
                }
            }
        } else {
            if outer_partner[pattern.duplex_5 -1] > inner_partner[pattern.duplex_3 - 1] {
                if outer_partner[pattern.duplex_5 -1] - inner_partner[pattern.duplex_3 - 1] == duplex_len - 1{
                    result.push(true);
                    continue;
                }
            } 
        }
     
        result.push(false);
    }
        
    return result  
}