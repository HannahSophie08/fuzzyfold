use std::process::{Command, Stdio};
use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::path::Path;
use std::path::PathBuf;
use clap::Parser;
use anyhow::Result;
use plotters::prelude::*;
use plotters::style::Palette;
use plotters::style::IntoFont;
use plotters::style::Color;

use fuzzyfold::input_parsers::read_fasta_like_input;
use ff_structure::DotBracketVec;
use ff_structure::PairList;
use ff_structure::PairTable;
use rustc_hash::FxHashMap;
use ff_alu_analysis::category::Category;


#[derive(Debug, Parser)]

pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: PathBuf,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[arg(long, value_name = "START-END", num_args = 1..)]
    regions: Vec<String>,

    #[arg(long, value_name = "NAME", num_args = 1..)]
    region_names: Vec<String>,

    #[arg(short, long, default_value_t = 100)]
    num_samples: usize,
}

fn main() -> Result<()> {
    
    let cli = Cli::parse();

    let (header, sequence, mut structure) = read_fasta_like_input(&cli.input)?;

    let regions: Vec<(usize, usize)> = cli.regions.iter()
    .map(|r| {
        let (a, b) = r.split_once('-').expect("region must be in format START-END");
        (a.parse::<usize>().expect("invalid start"), b.parse::<usize>().expect("invalid end"))
    })
    .collect();

    let csv_path = cli.output.with_extension("csv");
    let svg_path = cli.output.with_extension("svg");
    let mut writer = BufWriter::new(File::create(csv_path.clone())?);
    writeln!(writer, "length,category,value")?;
    
    let mut all_categories: Vec<FxHashMap<Category, f64>> = Vec::new();
    let mut lengths: Vec<usize> = Vec::new();

    for l in 1..sequence.len() {
        let subseq = sequence[..l];

        let structures = run_rnasubopt(subseq, cli.num_samples)?;
        if structures.is_empty() {
            continue
        }

        let categories = count_categories(&structures, &regions)?;

        write_row(&mut writer, l, &categories)?;
        lengths.push(l);
        all_categories.push(categories);

    }

    let region_ends: Vec<usize> = regions.iter().map(|&(_, end)| end).collect();

    plot_categories_over_length(&lengths, &all_categories, &cli.region_names, &region_ends, &svg_path);

    Ok(())
}


fn run_rnasubopt(sequence: &str, num_samples: usize) -> Result<Vec<DotBracketVec>> {

    let mut child = Command::new("RNAsubopt")
        .arg("-p")
        .arg(num_samples.to_string())
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;

    child.stdin
        .as_mut()
        .unwrap()
        .write_all(sequence.as_bytes())?;

    let output = child.wait_with_output()?;
    let stdout = String::from_utf8(output.stdout)?;

    let structures = stdout
        .lines()
        .skip(2)
        .map(|l| l.split_whitespace().next().unwrap().to_string())
        .collect();

    Ok(structures)
}


fn count_categories(structures: &[DotBracketVec], regions: &Vec<(usize, usize)>) -> Result<FxHashMap<Category, f64>>{

    let mut counts: FxHashMap<Category, usize> = FxHashMap::default();
    let mut total_pairs = 0;

    let mut skipped = 0;
    for structure in structures.iter() {

        let pl = PairList::try_from(structure)?;

        for (i, j) in pl.iter() {
            let region_i = regions.iter().position(|&(start, end) | (*i as usize) >= start && (*i as usize) <= end);
            let region_j = regions.iter().position(|&(start, end) | (*j as usize) >= start && (*j as usize) <= end);

            let category = match (region_i, region_j) {
                (Some(a), Some(b)) if a == b => Category::Within(a),
                (Some(a), Some(b)) => Category::Between(a.min(b), a.max(b)),
                (Some(a), None) => Category::WithRest(a),
                (None, Some(b)) => Category::WithRest(b),
                (None, None) => { skipped += 1; continue },
            };

            *counts.entry(category).or_insert(0) += 1;
            total_pairs += 1;
        
        }
    }
    eprintln!("Skipped {} pairs (both outside regions)", skipped);
    eprintln!("Counted {} pairs", total_pairs);


    let percentages = counts.into_iter()
        .map(|(cat, count)| (cat, count as f64 / total_pairs as f64 * 100.0))
        .collect();
    
    Ok(percentages)
}

fn write_row(
    writer: &mut impl Write,
    length: usize,
    categories: &FxHashMap<Category, f64>,
) -> Result<()> {
    for (cat, val) in categories {
        writeln!(writer, "{},{},{}", length, cat.to_key(), val)?;
    }
    Ok(())
}

fn region_name<'a>(names: &'a [String], idx: usize) -> String {
    names.get(idx)
        .cloned()
        .unwrap_or_else(|| format!("region {}", idx + 1))
}


fn plot_categories_over_length(
    lengths: &[usize],
    categories: &[FxHashMap<Category, f64>],
    region_names: &[String],
    region_ends: &[usize],
    filename: impl AsRef<Path>,
) {
    let l_max = *lengths.last().unwrap_or(&1);

    let root = SVGBackend::new(filename.as_ref(), (1024, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    root.titled("Category occupancy over transcript length", ("sans-serif", 28)).unwrap();
    root.draw_text(
        "length (nt)",
        &("sans-serif", 22).into_font().into_text_style(&root),
        (496, 450),
    ).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .margin_top(40)
        .margin_right(40)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(0f64..(l_max as f64), 0.0..1.0)
        .unwrap();

    chart
        .configure_mesh()
        .y_desc("occupancy")
        .x_labels(10)
        .y_labels(10)
        .light_line_style(RGBColor(220, 220, 220))
        .axis_desc_style(("sans-serif", 22))
        .label_style(("sans-serif", 18))
        .draw()
        .unwrap();

    // draw vertical dashed lines at region ends
    let dash_len = 0.04_f64;
    let gap_len  = 0.02_f64;
    for &re in region_ends {
        let x = re as f64;
        let mut y = 0.0_f64;
        while y < 1.0 {
            let y_end = (y + dash_len).min(1.0);
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(x, y), (x, y_end)],
                BLACK.mix(0.5).stroke_width(1),
            ))).unwrap();
            y = y_end + gap_len;
        }
    }

    // collect and sort categories
    let mut all_categories: Vec<&Category> = categories.iter()
        .flat_map(|m| m.keys())
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();

    all_categories.sort_by_key(|c| match c {
        Category::Within(a)     => (0, *a, 0),
        Category::Between(a, b) => (1, *a, *b),
        Category::WithRest(a)   => (2, *a, 0),
    });

    for (i, category) in all_categories.iter().enumerate() {
        let color = Palette99::pick(i).mix(0.9);

        let series: Vec<(f64, f64)> = lengths.iter().zip(categories.iter())
            .map(|(&l, map)| (l as f64, map.get(category).copied().unwrap_or(0.0) / 100.0))
            .collect();

        let label = match category {
            Category::Within(a)     => format!("Within {}", region_name(region_names, *a)),
            Category::Between(a, b) => format!("Between {} and {}", region_name(region_names, *a), region_name(region_names, *b)),
            Category::WithRest(a)   => format!("{} with rest", region_name(region_names, *a)),
        };

        chart.draw_series(LineSeries::new(
            series,
            color.stroke_width(2),
        )).unwrap()
            .label(label)
            .legend(move |(x, y)|
                PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
            );
    }

    chart
        .configure_series_labels()
        .border_style(BLACK)
        .background_style(WHITE.mix(0.8))
        .position(SeriesLabelPosition::UpperRight)
        .label_font(("sans-serif", 16).into_font())
        .draw()
        .unwrap();

    root.present().unwrap();
}