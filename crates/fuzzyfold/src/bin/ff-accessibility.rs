//! A command line tool to visualize nucleotide accessibility durin cotranscriptional folding.
//!
//! This tool runs stochastic cotranscriptional folding simulations, and tracks the accessibility 
//! at each position for each transcript length. The accessibility is the fraction of simulations
//! where a nucleotide is unpaired.
//! The accessibility is visualized as a heat map showing position vs. transcript length.
//! The simulation begins at transcript length 1 and one structure per transcript length is recorded. 
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
use ff_structure::DotBracketVec;
use ff_energy::Base;
use ff_energy::EnergyModel;
use ff_kinetics::LoopStructure;
use ff_kinetics::LoopStructureSSA;
use plotters::prelude::*;
use plotters::style::{RGBColor, Color};

use fuzzyfold::input_parsers::read_fasta_like_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;


#[derive(Debug, Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time and visualization of the accessibillity during transcription")]

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

    /// Output file name 
    #[arg(short, long, default_value = "accessibility.png")]
    output: String,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}


/// Matrix to store accessibility for all transcript lengths
struct AccessibilityMatrix {
    matrix: Vec<Vec<f64>>, 
    sequence_length: usize, 
}

impl AccessibilityMatrix {

    // new matrix
    fn new(sequence_length: usize) -> Self {
        Self {
            matrix: vec![vec![0.0; sequence_length]; sequence_length],
            sequence_length,
        }
    }

    // Add accessibility results of a single simulation
    fn add(&mut self, simulation: &[Vec<bool>]) {
        
        for transcript_len in 0..simulation.len() {
            for pos in 0..transcript_len {
                if simulation[transcript_len][pos] {
                    self.matrix[transcript_len][pos] += 1.0;
                }
            }
        }
    }

    // normalize matrix by number of simulations 
    fn normalize(&mut self, num_sims: usize) {
        for transcript_len in 0..self.sequence_length {
            for pos in 0..transcript_len {
                self.matrix[transcript_len][pos] /= num_sims as f64;
            }
        }
    }

}

// Runs a single cotranscriptional folding simulation
fn run_simulation(
    sequence: &[Base],
    pairings: &PairTable,
    times: &Vec<f64>,
    emodel: &impl EnergyModel,
    rmodel: &impl ff_kinetics::RateModel,
) -> Result<Vec<Vec<bool>>> {
    let loops = LoopStructure::try_from((&sequence[..], pairings, emodel)).unwrap();
    let mut simulator = LoopStructureSSA::from((loops, rmodel));

    let mut simulation_data = vec![vec![false; sequence.len()]; sequence.len()];
    simulation_data[0][0] = true;
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


// Map accessibility value to RGB color 
fn accessibility_color(accessibility: f64) -> RGBColor {
    if accessibility < 0.5 {
        let t = accessibility * 2.0;
        RGBColor(
            (30.0 + t * 170.0) as u8,
            (0.0 + t * 50.0) as u8,
            (50.0 - t * 30.0) as u8,
        )
    } else {
        let t = (accessibility - 0.5) * 2.0;
        RGBColor(
            (200.0 + t * 55.0) as u8,
            (50.0 + t * 180.0) as u8,
            (20.0 - t * 20.0) as u8,
        )
    }
}


// generate accessibility profile: heatmap visualizing nucleotide accessibility during contransriptional folding 
fn plot_accessibility (access_matrix: &AccessibilityMatrix, file_name: &str, title: &str) -> Result<()>{
     
    let seq_len = access_matrix.sequence_length;
    
    let width = (seq_len as u32) + 200;
    let height = (seq_len as u32) + 260; //extra space for legend
    
    let root = BitMapBackend::new(file_name, (width, height)).into_drawing_area();
    let (main_area, legend_area) = root.split_vertically(seq_len as u32 + 200);

    //heatmap chart
    let mut chart = ChartBuilder::on(&main_area)
        .caption(title, ("sans-serif", 30).into_font().color(&BLACK))
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0..seq_len, seq_len..0)?; //transcript length increases downward

    chart.configure_mesh()
        .x_desc("Nucleotide Position")
        .y_desc("Transcript Length")
        .x_labels(10)
        .y_labels(10)
        .axis_desc_style(("sans-serif", 20))
        .light_line_style(&WHITE)
        .draw()?;

    for transcript_len in 1..=seq_len {
        for pos in 0..transcript_len {
            let accessibility = access_matrix.matrix[transcript_len - 1][pos];
            let color = accessibility_color(accessibility);

            chart.draw_series(std::iter::once(Rectangle::new(
                [(pos, transcript_len - 1), (pos + 1, transcript_len)],
                color.filled(),
            )))?;
        }
    }

    //legend
    let mut legend_chart = ChartBuilder::on(&legend_area)
        .margin(10)
        .margin_left(200)
        .margin_right(200)
        .x_label_area_size(30)
        .build_cartesian_2d(0.0..1.0, 0..1)?;
    
    legend_chart.configure_mesh()
        .disable_y_mesh()
        .disable_y_axis()
        .x_labels(3)
        .x_label_formatter(&|x| format!("{:.1}", x))
        .x_desc("Accessibility")
        .axis_desc_style(("sans-serif", 18))
        .label_style(("sans-serif", 14))
        .draw()?;

    //color gradient legend
    let steps = 200;
    for i in 0..steps {
        let x_start = i as f64 / steps as f64;
        let x_end = (i + 1) as f64 / steps as f64;
        let color = accessibility_color(x_start);
        
        legend_chart.draw_series(std::iter::once(Rectangle::new(
            [(x_start, 0), (x_end, 1)],
            color.filled(),
        )))?;
    }

    root.present()?;
    Ok(())

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





  




    

    

