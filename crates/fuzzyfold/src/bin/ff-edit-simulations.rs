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
use ff_energy::EnergyModel;
use ff_energy::NucleotideVec;
use ff_energy::Base;
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
    let t_sep = cli.simulation.t_sep.expect("t-sep must exist after validation!");
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
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel.clone(), NoShift))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                    eprintln!("[moveset] selected=NoShift");
                run_timecourse(moves, rmodel, emodel, &sequence, &sim_times, cli.num_sims as u64, &times, NoShift).collect()
            },
            (true, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel.clone(), ThreeWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                     eprintln!("[moveset] selected=ThreeWayOnly");
                run_timecourse(moves, rmodel, emodel, &sequence, &sim_times, cli.num_sims as u64, &times, ThreeWayOnly).collect()
            },
            (false, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel.clone(), FourWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                     eprintln!("[moveset] selected=FourWayOnly");
                run_timecourse(moves, rmodel, emodel, &sequence, &sim_times, cli.num_sims as u64, &times, FourWayOnly).collect()
            },
            (true, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel.clone(), ThreeAndFour))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                     eprintln!("[moveset] selected=ThreeAndFour");
                run_timecourse(moves, rmodel, emodel, &sequence, &sim_times, cli.num_sims as u64, &times, ThreeAndFour).collect()
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
    writeln!(writer, "# sequence={}", sequence)?;
    writeln!(writer, "# num_sims={}", cli.num_sims)?;
    writeln!(writer, "# t_sep={:?}", t_sep)?;
    writeln!(writer, "# t_ext={:?}",  cli.simulation.t_ext)?;
    writeln!(writer, "# t_end={}", cli.simulation.t_end)?;
    writeln!(writer, "time,structure,count")?;

    for (t_idx, structures) in master.iter().enumerate() {
        let t = times[t_idx];
        for (s, count) in structures.iter() {
            writeln!(writer, "{},{},{}", t, s, count)?;
        }
    }
    println!("Wrote csv file: {}", format!("{}",csv_path.display()).green()); 

    Ok(())
}


fn run_timecourse<K, E, S>(
    moves: LoopNeighbors<E, S>,
    rmodel: K,
    emodel: Arc<E>,
    sequence: &NucleotideVec,
    sim_times: &[f64],
    num_sims: u64,
    times: &[f64],
    shift_policy: S,
) -> impl ParallelIterator<Item = Vec<DotBracketVec>>
where
    K: RateModel + Clone + Send + Sync,
    E: EnergyModel + Send + Sync + 'static,
    S: ShiftPolicy + Copy + Send + Sync,
    SSA<LoopNeighbors<E, S>, K>: From<(LoopNeighbors<E, S>, K)>,
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
            move |pb, _| {
                let mut timeline: Vec<DotBracketVec> = Vec::with_capacity(times.len());
                let mut current_sequence = sequence.clone();
                let mut current_moves = moves.clone();
                let mut current_sim_times = sim_times.to_vec();
                let mut t_idx = 0usize;
                let mut t_offset = 0.0;
                let l0 = sequence.len() - sim_times.len() + 1;

                loop {
                    let mut simulator = SSA::from((current_moves.clone(), rmodel.clone()));
                    let mut edit_info: Option<(DotBracketVec, NucleotideVec, f64)> = None;
                    simulator.co_simulate(
                        &mut rng(),
                        &current_sim_times,
                        |t, tinc, _, w| {
                            let total_t = t + t_offset;
                            while t_idx < times.len() && total_t + tinc >= times[t_idx] {
                                timeline.push(w.current_structure());
                                let (edited, edited_sequence) = edit_or_move_on(&current_sequence, &w.current_structure());
                                if edited {
                                    edit_info = Some((w.current_structure(), edited_sequence, t));
                                    return false;
                                }
                                t_idx += 1;

                            }
                            true
                        },
                    );
                    match edit_info {
                        Some((structure, new_sequence, t)) => {
                            let pairings = PairTable::try_from(&structure).expect("failed to build PairTable from edited structure");
                            let transcript_len = structure.len();
                            current_sequence = new_sequence;
                            t_offset = t;
                            current_moves = LoopNeighbors::try_from((current_sequence.clone(), &pairings, emodel.clone(), shift_policy)).expect("failed to construct AddDelMoves");
                            let used = transcript_len - l0;
                            if sim_times[used] < t {
                                let remaining_time = sim_times[used] - t;
                                let mut new_sim_times = Vec::with_capacity(1 + sim_times.len() - used - 1);
                                new_sim_times.push(remaining_time);
                                new_sim_times.extend_from_slice(&sim_times[used + 1..]);
                                current_sim_times = new_sim_times;
                            } else {
                                current_sim_times = sim_times[used..].to_vec();
                            }
                            continue;
                        },
                        None => break,
                    };
                }
                
                let n = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                if n % 10 == 0 || n == num_sims {
                    eprintln!("progress: {}/{}", n, num_sims);
                }

                pb.inc(1);
                timeline
            },
        )
}

fn random_choice (probability: f64) -> bool {
    if rand::random::<f64>() < probability {
        return true;
    } 
    return false;
}

fn edit_or_move_on (sequence: &NucleotideVec, structure: &DotBracketVec) -> (bool, NucleotideVec) {

    let (adenosines, editing) = editing(sequence, structure);
    let mut edited_sequence = sequence.clone();
    let mut bases: Vec<Base> = sequence.iter().copied().collect();
    let mut edited = false;

    for (i, a) in adenosines.iter().enumerate() {
        if editing[i] {
            bases[*a] = Base::I;
            edited_sequence = NucleotideVec::from(ff_energy::NucleotideVec(bases.clone())); 
            println!("Edited position {} at transcript length {}", a, structure.len());
            edited = true;
        }
    }

    return (edited, edited_sequence)
}

fn editing (sequence: &NucleotideVec, structure: &DotBracketVec) -> (Vec<usize>, Vec<bool>) {

    let adenosines = sequence.iter().enumerate()
            .filter(|(_, base)| **base == Base::A)
            .map(|(i, _)| i)
            .collect();

    let edit_11 = position_check(sequence, structure, &adenosines, 1, 1);
    let edit_36 = position_check(sequence, structure, &adenosines, 3, 6);
    let edit_58 = position_check(sequence, structure, &adenosines, 5, 8);

    let mut edited = Vec::new();
    for (i, _a) in adenosines.iter().enumerate() {
        let mut choice = false;
        if edit_58[i] {
            choice = random_choice(0.3);
        } else if edit_36[i] {
            choice = random_choice(0.2);
        }  else if edit_11[i] {
            choice = random_choice(0.1);
        }
        edited.push(choice);
    }
    return (adenosines, edited)
}
// checks whether position is unpaired and across from a 'C', and whether the 
// 5' and 3' neighbors are paired
fn position_check (sequence: &NucleotideVec, structure: &DotBracketVec, positions: &Vec<usize>, duplex_5: usize, duplex_3: usize) -> Vec<bool> {

    let mut result = Vec::new();
    let pt = PairTable::try_from(structure).unwrap();

    'outer_pos: for position in positions.iter() {

        if *position < duplex_5 || *position + duplex_3 >= structure.len() {
            result.push(false);
            continue;
        }
            
        if structure[*position] != DotBracket::Unpaired {
            result.push(false);
            continue;
        }

        let mut outer = Vec::new();
        for i in 1..=duplex_5 {
            outer.push((*position - i) as usize);
        }

        let mut inner = Vec::new();
        for i in 1..=duplex_3 {
            inner.push((*position + i) as usize);
        }

        let mut outer_partner = Vec::new();
        for o in outer {
            if !pt[o].is_some() {
                result.push(false );
                continue 'outer_pos; 
            } else {
                outer_partner.push(pt[o].unwrap() as usize);
            }
        } 

        let mut inner_partner = Vec::new();
        for i in inner {
            if !pt[i].is_some() {
                result.push(false );
                continue 'outer_pos; 
            } else {
                inner_partner.push(pt[i].unwrap() as usize);
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

        let duplex_len = duplex_5 + duplex_3;

        if outer_partner[duplex_5 -1] > inner_partner[duplex_3 - 1] {
            if outer_partner[duplex_5 -1] - inner_partner[duplex_3 - 1] == duplex_len && sequence[outer_partner[0] - 1] == Base::C {
                result.push(true);
                continue;
            }
        }
     
        result.push(false);
    }
        
    return result  
}