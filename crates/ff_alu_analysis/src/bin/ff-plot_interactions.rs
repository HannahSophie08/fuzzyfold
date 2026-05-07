
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::path::PathBuf;


use clap::Parser;
use anyhow::Result;
use rustc_hash::FxHashMap;
use plotters::prelude::*;
use plotters::style::Palette;
use plotters::style::IntoFont;
use plotters::prelude::IntoLogRange;
use plotters::style::Color;


use ff_alu_analysis::category::Category;


#[derive(Debug, Parser)]
#[command(version, about = "Plot category occupancy from ff-interactions CSV output.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: PathBuf,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[arg(long)]
    t_split: f64,

    #[arg(long, value_name = "NAME", num_args = 0..)]
    region_names: Vec<String>,

}



fn main() -> Result<()> {

    let cli = Cli::parse();

    let (t_split, t_end, times, categories) = read_csv(&cli.input)?; 


    let all_categories_path = cli.output.with_extension("_all_categories.svg");
    let between_categories_path =  cli.output.with_extension("between_categories.svg");

    plot_all_categorie_over_time(&times, &categories, t_split, t_end, all_categories_path.clone());

    let between_categories = renormalize_between(&categories);
    plot_betweens(&times, &between_categories, t_split, t_end, between_categories_path.clone());

    Ok(())
}


fn read_csv(path: &PathBuf) ->  Result<(f64, f64, Vec<f64>, Vec<FxHashMap<Category, f64>>)> {

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

     // Read metadata lines
    let t_split: f64 = lines.next().unwrap()?
        .trim_start_matches("# t_split=")
        .parse()?;
    let t_end: f64 = lines.next().unwrap()?
        .trim_start_matches("# t_end=")
        .parse()?;

    lines.next(); //skip header

    let mut rows: Vec<(f64, Category, f64)> = Vec::new();
    for line in lines {
        let line = line?;
       let mut parts = line.splitn(3, ',');
        let time: f64 = parts.next().unwrap().parse()?;
        let category = Category::from_key(parts.next().unwrap())?;
        let value: f64 = parts.next().unwrap().parse()?;
        rows.push((time, category, value)); 
    }

    let mut times: Vec<f64> = Vec::new();
    let mut categories: Vec<FxHashMap<Category, f64>> = Vec::new();

    for (time, category, value) in rows {
        if times.last().copied() != Some(time) {
            times.push(time);
            categories.push(FxHashMap::default());
        }
        categories.last_mut().unwrap().insert(category, value);
    }

    Ok((t_split, t_end, times, categories))
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




fn renormalize_between(categories: &[FxHashMap<Category, f64>]) -> Vec<FxHashMap<Category, f64>> {
    categories.iter().map(|map| {
        let between_sum: f64 = map.iter()
            .filter(|(cat, _)| matches!(cat, Category::Between(_, _)))
            .map(|(_, v)| v)
            .sum();

        map.iter()
            .filter(|(cat, _)| matches!(cat, Category::Between(_, _)))
            .map(|(cat, v)| (
                cat.clone(),
                if between_sum > 0.0 { v / between_sum * 100.0 } else { 0.0 }
            ))
            .collect()
    }).collect()
}



fn plot_betweens(
    times: &[f64],
    categories: &[FxHashMap<Category, f64>],
    t_split: f64,
    t_end: f64,
    filename: impl AsRef<Path>,
) {
    assert!(t_split > 0.0 && t_end > t_split, "Require 0 < t_split < t_end");

    let root = SVGBackend::new(filename.as_ref(), (1024, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    root.titled("Between-region occupancy over time", ("sans-serif", 28)).unwrap();
    root.draw_text(
        "time",
        &("sans-serif", 22).into_font().into_text_style(&root),
        (496, 450),
    ).unwrap();

    let eps = 1e-9;
    let (left, right) = root.split_horizontally(512);

    let mut chart_left = ChartBuilder::on(&left)
        .caption("Linear plot", ("sans-serif", 18))
        .margin(20).margin_top(40).margin_right(0)
        .x_label_area_size(40).y_label_area_size(50)
        .build_cartesian_2d(0.0..(t_split + eps), 0.0..1.0).unwrap();
    chart_left.configure_mesh()
        .y_desc("occupancy").x_labels(6).y_labels(10)
        .light_line_style(RGBColor(220, 220, 220))
        .axis_desc_style(("sans-serif", 22))
        .label_style(("sans-serif", 18))
        .draw().unwrap();
    chart_left.draw_series(std::iter::once(PathElement::new(
        vec![(t_split, 0.0), (t_split, 1.0)], BLACK.mix(0.7),
    ))).unwrap();

    let mut chart_right = ChartBuilder::on(&right)
        .caption("Logarithmic plot", ("sans-serif", 18))
        .margin(20).margin_top(40).margin_left(0).margin_right(40)
        .x_label_area_size(40).y_label_area_size(0)
        .build_cartesian_2d(((t_split - eps)..(t_end + eps)).log_scale(), 0.0..1.0)
        .unwrap();
    chart_right.configure_mesh()
        .x_labels(6)
        .x_label_formatter(&|x| if *x < 0.01 { format!("{:.1e}", x) } else { format!("{}", x) })
        .y_labels(10)
        .light_line_style(RGBColor(220, 220, 220))
        .label_style(("sans-serif ", 18))
        .draw().unwrap();
    chart_right.draw_series(std::iter::once(PathElement::new(
        vec![(t_split, 0.0), (t_split, 1.0)], BLACK.mix(0.7),
    ))).unwrap();

    // Collect and sort only Between categories
    let mut between_categories: Vec<&Category> = categories.iter()
        .flat_map(|m| m.keys())
        .filter(|c| matches!(c, Category::Between(_, _)))  // <-- only Betweens
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();

    between_categories.sort_by_key(|c| match c {
        Category::Between(a, b) => (*a, *b),
        _ => unreachable!(),
    });

    for (i, category) in between_categories.iter().enumerate() {
        let color = Palette99::pick(i).mix(0.9);

        let series: Vec<(f64, f64)> = times.iter().zip(categories.iter())
            .map(|(&t, map)| (t, map.get(category).copied().unwrap_or(0.0) / 100.0))
            .collect();

        let Category::Between(a, b) = category else { unreachable!() };
        let label = format!("Between region {} and {}", a + 1, b + 1);

        chart_left.draw_series(LineSeries::new(
            series.iter().cloned().filter(|(t, _)| *t <= t_split + eps),
            color.stroke_width(2),
        )).unwrap();

        chart_right.draw_series(LineSeries::new(
            series.iter().cloned().filter(|(t, _)| *t >= t_split - eps),
            color.stroke_width(2),
        )).unwrap()
            .label(label)
            .legend(move |(x, y)|
                PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
            );
    }

    chart_right
        .configure_series_labels()
        .border_style(BLACK)
        .background_style(WHITE.mix(0.8))
        .position(SeriesLabelPosition::UpperRight)
        .label_font(("sans-serif", 16).into_font())
        .draw().unwrap();

    root.present().unwrap();
}
