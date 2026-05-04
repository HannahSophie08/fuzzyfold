//! Stochastic folding simulator 
//! 
//! This binary performs stochastic folding simulations both at full sequence length and 
//! cotrnacriptional, running multiple simulations in parallel and merging them into a master timeline. 
//! This timeline is used to produce an occupancy plot. The plot shows the occupancy of predefined 
//! macrostates in average over time. 
//! 
//! Mode: 
//! Full length simulation: give structure of full sequence length 
//! Cotranscriptonal simulation: give no structure or structure shorter than sequence length
//! 
//! --- Full length simulation --- 
//! 
//! Parameters:
//!  
//! - t-end: total simulationt time 
//! - t-sep: last timepoint on linear scale (timepoint where linear scale in occupany plot switches
//! to logarithmic scale and seperator is placed)
//! - t-lin: timepoints on linear time scale [0...t-sep]
//! - t-log: time points on logarithmic timescale [t-sep...t-log]
//! - num-sims: number of simulations performed
//! - k0: kinetic rate constant 
//! 
//! Input
//! - user has to give structure and t-sep
//! - t-lin and t-log are optional (default: 100)
//! - num-sim optional (default: 1)
//! 
//! --- Cotranscriptional simulation ---
//! 
//! Parameters
//!  
//! - t-ext: extension time (simulation time per transcript length)
//! - t-end: time of postranscriptional simulation
//! - t-sep: last timepoint on linear scale (timepoint where linear scale in occupany plot switches
//! to logarithmic scale and seperator is placed)
//! - t-lin: recorded time steps per transcript length [without t-sep]/ timepoints recorded on linear timescale [with t-sep]
//! - t-log: time points recorded for the posttranscriptional simulation on a logarithmic timescale [without t-sep]/
//! timepoints recorded on logarithmic timescale [with t-sep]
//! - num-sims: number of simulations performed
//! - k0: kinetic rate constant 
//! 
//! Input & mode 
//! - Cotranscriptional simulation 
//! -> user has to give either no structure or a structure, that is shorter then full length, that 
//!     will be used as the start structure and t-ext (extension time)
//! - optional: t-sep 
//! -> no t-sep: linear scale ends at end of transcription (t-lin: timepoints recorded per extension step;
//! t-log: timepoints recorded during posttranscriptional folding)
//! -> t-sep: linear timescale ends at t-sep (t-lin: timepoints recorded on linear timescale; t-log: time points recorded on
//! logarithmic timescale)
//! - optional: t-lin and t-log are optional (default: 100)
//! - optional: num-sim (default: 1)
//! 
//! --- General ---
//! 
//! Output: 
//! - user has to give output file
//! - output formats: 
//! -> svg (occupancy plot)
//! -> nxy (occupancies)
//! -. tln (timeline)
//! 
//! Load timeline:
//! The 'tln' file can be used to accumulate results incrementially. Existing timlines can be reloaded and when a simulation is
//! run with it, it gets updated. (Parameters have to be set to the same values!)
//! 



use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::sync::Arc;
use std::path::Path;
use std::path::PathBuf;

use ff_structure::DotBracketVec;
use ff_structure::PairList;
use rayon::prelude::*;
use rand::rng;
use colored::*;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use rustc_hash::FxHashMap;
use plotters::prelude::*;
use plotters::style::Palette;
use plotters::style::IntoFont;
use plotters::prelude::IntoLogRange;
use plotters::style::Color;

use ff_structure::PairTable;
use ff_energy::EnergyModel;
use ff_kinetics::RateModel;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::shift_policy::*;
use ff_kinetics::SSA;


use ff_alu_analysis::input_parsers::read_cotr_input;
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

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[arg(long, value_name = "START-END", num_args = 1..)]
    regions: Vec<String>,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimelineParameters,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

impl Category {
    fn to_key(&self) -> String {
        match self {
            Category::Within(a)      => format!("Within_{}", a + 1),
            Category::Between(a, b)  => format!("Between_{}_{}", a + 1, b + 1),
            Category::WithRest(a)    => format!("WithRest_{}", a + 1),
        }
    }
}

#[derive(Clone, Eq, Hash, PartialEq)]
enum Category {
    Within(usize),
    Between(usize, usize),
    WithRest(usize),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.simulation.validate()?;

    // --- Build simulator ---
    let emodel = Arc::new(cli.energy.build_model());
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) = read_cotr_input(&cli.input, is_rna)?;
    let pairings = PairTable::try_from(&structure)?;

    let regions: Vec<(usize, usize)> = cli.regions.iter()
    .map(|r| {
        let (a, b) = r.split_once('-').expect("region must be in format START-END");
        (a.parse::<usize>().expect("invalid start"), b.parse::<usize>().expect("invalid end"))
    })
    .collect();

    let svg_path = cli.output.with_extension("svg");
    let csv_path = cli.output.with_extension("csv");

    // --- Check input ---
    if structure.len() < sequence.len() && cli.simulation.t_ext.is_none() {
        panic!("Error:'t-ext' (extension time) missing for cotranscriptional simulation!"); 
    }

    if structure.len() == sequence.len() && cli.simulation.t_ext.is_some() {
        panic!("Error: 't-ext' (extension time) given in combination full length structure. For cotranscriptional simulation give either 
        shorter start structure or no structure!");
    }


    if let Some(h) = header {
        println!("{}", h.yellow());
    }

    println!("{}", sequence);

    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    let times = cli.simulation.get_output_times(&sequence, Some(&structure)); // build times vector for output times 

    let mut sim_times = vec![cli.simulation.t_ext.unwrap(); sequence.len() - structure.len()];
    sim_times.push(cli.simulation.t_end);

    let all_structures: Vec<Vec<DotBracketVec>> =
        match (rmodel.k3ws().is_some(), rmodel.k4ws().is_some()) {
            (false, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, NoShift))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
            (true, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
            (false, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, FourWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
            (true, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeAndFour))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64, &times).collect()
            },
        };

    
    let sorted_structures = sort_by_timepoint(&times, &all_structures, cli.num_sims);

    let categories =  count_categories(&sorted_structures, &regions)?;

    let co_time = cli.simulation.t_ext.unwrap() * (sequence.len() - structure.len()) as f64;
    let post_time = co_time + cli.simulation.t_end;

    plot_all_categorie_over_time(&times, &categories, co_time, post_time, svg_path.clone());

    println!("{}", "Finished simulations!".red());

    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
    write_categories(&mut writer, &times, &categories)?;
   
    
        
    println!("Plotted svg file: {}", svg_path.display());


    Ok(())
}




/// Cotranscriptional folding simulation 
fn run_cotimecourse<W, K>(
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
           
            let mut structures: Vec<DotBracketVec> = Vec::with_capacity(times.len());

            let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
            let mut t_idx = 0;
            simulator.co_simulate(
                &mut rng(),
                sim_times,
                |t, tinc, _, w| {
                    while t_idx < times.len() && t + tinc >= times[t_idx] {
                        structures.push(w.current_structure());
                        t_idx += 1;
                    }
                    true
                },
            );

            pb.inc(1);
            structures
        },
    )

}

fn sort_by_timepoint(times: &Vec<f64>, structures: &Vec<Vec<DotBracketVec>>, num_sims: usize) -> Vec<Vec<DotBracketVec>> {

    let num_times = times.len();
    let mut sorted: Vec<Vec<DotBracketVec>> = vec![Vec::with_capacity(num_sims); num_times];

    for s in structures.iter() {
        for (t_idx, structure) in s.iter().enumerate() {
            sorted[t_idx].push(structure.clone());
        }
    }
    
    sorted
}

fn count_categories(structures: &[Vec<DotBracketVec>], regions: &Vec<(usize, usize)>,) -> Result<Vec<FxHashMap<Category, f64>>>{

    let mut result = Vec::with_capacity(structures.len());
   
    for timepoint_structures in structures.iter() {

        let mut counts: FxHashMap<Category, usize> = FxHashMap::default();
        let mut total_pairs = 0; // total pairs in region 


        for structure in timepoint_structures.iter() {

            let pl = PairList::try_from(structure)?;

            for (i, j) in pl.iter() {
                let region_i = regions.iter().position(|&(start, end) | (*i as usize) >= start && (*i as usize) <= end);
                let region_j = regions.iter().position(|&(start, end) | (*j as usize) >= start && (*j as usize) <= end);

                let category = match (region_i, region_j) {
                    (Some(a), Some(b)) if a == b => Category::Within(a),
                    (Some(a), Some(b)) => Category::Between(a.min(b), a.max(b)),
                    (Some(a), None) => Category::WithRest(a),
                    (None, Some(b)) => Category::WithRest(b),
                    (None, None) => continue,
                };

                *counts.entry(category).or_insert(0) += 1;
                total_pairs += 1;
            
            }
        }

        let percentages = counts.into_iter()
            .map(|(cat, count)| (cat, count as f64 / total_pairs as f64 * 100.0))
            .collect();

        result.push(percentages); 
    }
    Ok(result)
}


fn plot_all_categorie_over_time(
    times: &[f64],
    categories: &[FxHashMap<Category, f64>],
    t_split: f64,
    t_end: f64,
    filename: impl AsRef<Path>,
    ) {
    
    assert!(t_split > 0.0 && t_end > t_split, "Require 0 < t_split < t_end");

    // Image size; tweak as you like
    //let root = BitMapBackend::new(filename, (1024, 480)).into_drawing_area();
    let root = SVGBackend::new(filename.as_ref(), (1024, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    root.titled("Category occupancy over time", ("sans-serif", 28)).unwrap();
    root.draw_text(
        "time",
        &("sans-serif", 22).into_font().into_text_style(&root),
        (496, 450),   // roughly centered at bottom
    ).unwrap();


    let eps = 1e-9; // epsilon for plot labels
    // Split into two panels: 50% for linear (left), 50% for log (right)
    let (left, right) = root.split_horizontally(512);

    // ---- Left: linear panel ----
    let mut chart_left = ChartBuilder::on(&left)
        .caption("Linear plot", ("sans-serif", 18))
        .margin(20)
        .margin_top(40)
        .margin_right(0)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..(t_split+eps), 0.0..1.0).unwrap();
    chart_left
        .configure_mesh()
        //.x_desc("liner scale")
        .y_desc("occupancy")
        .x_labels(6)
        .y_labels(10)
        .light_line_style(RGBColor(220, 220, 220))
        .axis_desc_style(("sans-serif", 22))
        .label_style(("sans-serif", 18))
        .draw()
        .unwrap();

    // draw separator at x = t_lin (right edge of this panel)
    chart_left.draw_series(std::iter::once(PathElement::new(
        vec![(t_split, 0.0), (t_split, 1.0)],
        BLACK.mix(0.7),
    ))).unwrap();

    // ---- Right: log panel ----
    let mut chart_right = ChartBuilder::on(&right)
        .caption("Logarithmic plot", ("sans-serif", 18))
        .margin(20)
        .margin_top(40)
        .margin_left(0)
        .margin_right(40)
        .x_label_area_size(40)
        .y_label_area_size(0) // hide y labels on right
        .build_cartesian_2d(((t_split - eps)..(t_end + eps)).log_scale(), 0.0..1.0)
        .unwrap();

    chart_right
        .configure_mesh()
        //.x_desc("log scale")
        .x_labels(6)
        .x_label_formatter(&|x| if *x < 0.01 {format!("{:.1e}", x)} else {format!("{}", x)})  // scientific notation
        .y_labels(10) // hide y ticks on right
        .light_line_style(RGBColor(220, 220, 220))
        .label_style(("sans-serif", 18))
        .draw().unwrap();

    // repeat separator at x = t_lin (left edge of this panel)
    chart_right.draw_series(std::iter::once(PathElement::new(
        vec![(t_split, 0.0), (t_split, 1.0)],
        BLACK.mix(0.7),
    ))).unwrap();


    
    let mut all_categories: Vec<&Category> = categories.iter()
        .flat_map(|m| m.keys())
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();

    all_categories.sort_by_key(|c| match c {
        Category::Within(a) => (0, *a, 0),
        Category::Between(a, b) => (1, *a, *b),
        Category::WithRest(a) => (2, *a, 0),
    });


    // Find global Y max for normalization
    for (i, category) in all_categories.iter().enumerate() {
        let color = Palette99::pick(i).mix(0.9); // pick a distinct color

        let series: Vec<(f64, f64)> = times.iter().zip(categories.iter())
            .map(|(&t, map) | (t, map.get(category).copied().unwrap_or(0.0) / 100.0))
            .collect();

        let label = match category {
            Category::Within(a) => format!("Within region {}", a + 1),
            Category::Between(a, b) => format!("Between region {} and {}", a + 1, b + 1),
            Category::WithRest(a) => format!("Region {} with rest", a + 1),
        };


        chart_left.draw_series(LineSeries::new(
                series.iter().cloned().filter(|(t, _)| *t <= t_split + eps),
                color.stroke_width(2),
        )).unwrap();


        chart_right.draw_series(LineSeries::new(
            series.iter().cloned().filter(|(t, _)| *t >= t_split - eps),
            color.stroke_width(2),
        )).unwrap()
            .label(label)   // <-- label for legend
            .legend(move |(x, y)| 
                PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
            );
    }
    
    // after loop:
    chart_right
        .configure_series_labels()
        .border_style(BLACK)
        .background_style(WHITE.mix(0.8))
        .position(SeriesLabelPosition::UpperRight)
            .label_font(("sans-serif", 16).into_font())   // <-- legend font size
        .draw().unwrap();
    
    root.present().unwrap(); // write the PNG
}


fn write_categories( 
    writer: &mut impl Write,
    times: & [f64],
    categories: &[FxHashMap<Category, f64>],
) -> Result<()> {
    writeln!(writer, "time,category,value")?;
    for (t, map) in times.iter().zip(categories.iter()) {
        for (cat, val) in map {
            writeln!(writer, "{},{},{}", t, cat.to_key(), val)?;
        }
    }
    Ok(())
}