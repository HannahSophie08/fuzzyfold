use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::sync::Arc;
use std::path::Path;
use std::path::PathBuf;

use rand::seq;
use rayon::prelude::*;
use rand::rng;
use colored::*;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use serde_json::to_string_pretty;

use ff_structure::PairTable;
use ff_energy::EnergyModel;
use ff_kinetics::RateModel;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::shift_policy::*;
use ff_kinetics::SSA;
use ff_kinetics::timeline::Timeline;
use ff_kinetics::timeline_plotting::plot_occupancy_over_time;
use ff_kinetics::MacrostateRegistry;

use ff_alu_analysis::input_parsers::read_eval_input;
use ff_alu_analysis::energy_parsers::EnergyModelArguments;
use ff_alu_analysis::kinetics_parsers::RateModelArguments;
use ff_alu_analysis::kinetics_parsers::TimelineParameters;

#[derive(Debug, Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(long, value_name = "FILE", required = false)]
    annotations: Vec<PathBuf>,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimelineParameters,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

pub struct Alu {
    name: String,
    start: usize,
    end: usize,
    orientation: bool, // inverted (?)
}

impl Alu {

    fn new(name: String, start: usize, end: usize, orientation: bool) -> Self{
        Self {
            name,
            start,
            end,
            orientation,
        }
    }
}


pub struct pair_record {
    alus: Vec<Alu>,
    matrix: 
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.simulation.validate()?;

    // --- Build simulator ---
    let emodel = Arc::new(cli.energy.build_model());
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) = read_eval_input(&cli.input, is_rna)?;
    let pairings = PairTable::try_from(&structure)?;

    if let Some(h) = header {
        println!("{}", h.yellow());
    } 
    println!("{}", sequence);

    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    let times = cli.simulation.get_output_times(&sequence, &structure);
    

    // Output paths 
    let tln_path = cli.output.with_extension("tln");
    let svg_path = cli.output.with_extension("svg");
    let nxy_path = cli.output.with_extension("nxy");

    
    println!("{}", "Finished simulations!".red());

    // save / print / plot.
    let mut writer = BufWriter::new(File::create(nxy_path.clone())?);
    write!(writer, "{}", master)?;
    println!("Wrote nxy file: {}", format!("{}",nxy_path.display()).green());
    plot_occupancy_over_time(&master, svg_path.clone(), cli.simulation.t_ext, cli.simulation.t_end);
    println!("Plotted svg file: {}", svg_path.display());
    let serial = master.to_serializable();
    let json = to_string_pretty(&serial).unwrap();
    fs::write(tln_path.clone(), json).unwrap();
    println!("Wrote tln file: {}", tln_path.display());

    Ok(())
}


fn run_timecourse<W, K, E>(
    moves: W,
    rmodel: K,
    simulation_times: &[f64],
    num_sims: u64,
    registry: Arc<MacrostateRegistry<E>>,
    times: &[f64],
) -> impl ParallelIterator<Item = Timeline<E>>
where
    W: Walker + Clone + Send + Sync,
    K: RateModel + Clone + Send + Sync,
    E: EnergyModel + Send + Sync,
    SSA<W, K>: From<(W, K)>,
{   
    let pb = ProgressBar::new(num_sims);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    if simulation_times.len() > 1 { // Cotranscriptional simulation 

        (0..num_sims)
            .into_par_iter()
            .map_init(
                move || pb.clone(), // each thread gets a clone
                move |pb, _| {

                    //let registry = Arc::clone(&registry);
                    //let mut timeline = Timeline::new(times, registry);

                    let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
                    let mut t_idx = 0;
                    simulator.co_simulate(
                        &mut rng(),
                        simulation_times,
                        |t, tinc, _, w| {
                            while t_idx < times.len() && t + tinc >= times[t_idx] {
                                let structure = w.current_structure();
                                //timeline.assign_structure(t_idx, &structure);
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
    

    (0..num_sims)
        .into_par_iter()
        .map_init(
            move || pb.clone(), // each thread gets a clone
            move |pb, _| {
                let registry = Arc::clone(&registry);
                let mut timeline = Timeline::new(times, registry);

                let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
                let mut t_idx = 0;
                simulator.simulate(
                    &mut rng(),
                    t_end,
                    |t, tinc, _, w| {
                        while t_idx < times.len() && t + tinc >= times[t_idx] {
                            let structure = w.current_structure();
                            timeline.assign_structure(t_idx, &structure);
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
 