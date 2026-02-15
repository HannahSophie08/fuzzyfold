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
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time and visualitualize accessibillity during transcription")]

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

    /// Timepoints per sequence length 
    #[arg(long, default_value_t = 1)] 
    t_lin: usize,

    /// Output file name 
    #[arg(short, long, default_value = "accessibility.png")]
    output: String,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}


//Matrix to store accessibility
struct AccessibilityMatrix {
    data: Vec<Vec<f64>>,
    sequence_length: usize,
}

impl AccessibilityMatrix {

    fn new(sequence_length: usize) -> Self {
        Self {
            data: vec![vec![0.0; sequence_length]; sequence_length],
            sequence_length,
        }
    }

    fn add(&mut self, simulation: &[Vec<bool>]) {
        
        for transcript_len in 0..simulation.len() {
            for pos in 0..simulation[transcript_len].len() {
                if simulation[transcript_len][pos] {
                    self.data[transcript_len][pos] += 1.0;
                }
            }
        }
    }

    fn normalize(&mut self, num_sims: usize) {
        for transcript_len in 0..self.sequence_length {
            for pos in 0..self.sequence_length {
                self.data[transcript_len][pos] /= num_sims as f64;
            }
        }
    }

}

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
    let mut last_recorded_length = 0;

    simulator.co_simulate(
        &mut rng(),
        times.clone(),
        |_t, tinc, _, ls| {
            if tinc == 0.0 {
                let structure = DotBracketVec::from(ls);
                let transcript_len = structure.len();
                let pt = PairTable::try_from(&structure).unwrap();
                
                if transcript_len != last_recorded_length && transcript_len <= sequence.len() {
                    
                    for pos in 0..transcript_len {
                        simulation_data[transcript_len - 1][pos] = match pt.get(pos) {
                            Some(None) => true,      // Unpaired (accessible)
                            Some(Some(_)) => false,  // Paired (unaccessible)
                            None => false,           
                        };
                    }
                    
                    last_recorded_length = transcript_len;
                }
            }
            true
        }
    );
    for pos in 0..1 {
        simulation_data[0][pos] = true;
    }

    Ok(simulation_data)
}



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

fn plot_accessibility (matrix: &AccessibilityMatrix, file_name: &str, title: &str) -> Result<()>{

    let root = BitMapBackend::new(file_name, (820, 860)).into_drawing_area();
    root.fill(&WHITE)?;

    let seq_len = matrix.sequence_length;

    let (main_area, legend_area) = root.split_vertically(800); 

    let mut chart = ChartBuilder::on(&main_area)
        .caption(title, ("sans-serif", 30).into_font().color(&BLACK))
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0..seq_len, seq_len..0)?;

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
            let accessibility = matrix.data[transcript_len - 1][pos];
            let color = accessibility_color(accessibility);

            chart.draw_series(std::iter::once(Rectangle::new(
                [(pos, transcript_len - 1), (pos + 1, transcript_len)],
                color.filled(),
            )))?;
        }
    }
    
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

    // Build models
    let emodel = cli.energy.build_model();
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let (_header, sequence, _structure) = read_fasta_like_input(&cli.input)?;
    let initial_structure = DotBracketVec::try_from(".")?;
    let pairings = PairTable::try_from(&initial_structure)?;
    let t_ext = cli.t_ext; //extension time (simulation time per nucleotide)
    let t_lin = cli.t_lin;
    let mut times = Vec::new();
    for i in 0..(sequence.len()) {
        times.push((i+1) as f64 * t_ext);
    }

    println!("Simulation progress:");
    let pb = ProgressBar::new(cli.num_sims as u64);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    let num_sims = cli.num_sims;
    let mut matrix = AccessibilityMatrix::new(sequence.len());

    let simulations_results: Vec<_> = (0..cli.num_sims)
        .into_par_iter()
        .map(|_| {
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





  




    

    

