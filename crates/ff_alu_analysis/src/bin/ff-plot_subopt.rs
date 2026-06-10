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

    #[arg(long, value_name = "NAME", num_args = 1..)]
    region_names: Vec<String>,

    #[arg(long, value_name = "TIME", num_args = 1..)]
    region_ends: Vec<f64>,
}



fn main() -> Result<()> {
    let cli = Cli::parse();

    let (lengths, categories) = read_csv(&cli.input)?;

    let path = cli.output.with_extension(".svg");

    plot_categories_over_length(
        &lengths, &categories, &cli.region_names, &cli.region_ends,path,
    );


    Ok(())
}

fn region_name<'a>(names: &'a [String], idx: usize) -> String {
    names.get(idx)
        .cloned()
        .unwrap_or_else(|| format!("region {}", idx + 1))
}


fn read_csv(path: &PathBuf) -> Result<(Vec<f64>, Vec<FxHashMap<Category, f64>>)> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    lines.next(); // skip header

    let mut rows: Vec<(f64, Category, f64)> = Vec::new();
    for line in lines {
        let line = line?;
        let mut parts = line.splitn(3, ',');
        let length: f64  = parts.next().unwrap().parse()?;
        let category   = Category::from_key(parts.next().unwrap())?;
        let value: f64 = parts.next().unwrap().parse()?;
        rows.push((length, category, value));
    }

    let mut lengths: Vec<f64> = Vec::new();
    let mut categories: Vec<FxHashMap<Category, f64>> = Vec::new();

    for (length, category, value) in rows {
        if lengths.last().copied() != Some(length) {
            lengths.push(length);
            categories.push(FxHashMap::default());
        }
        categories.last_mut().unwrap().insert(category, value);
    }

    Ok((lengths, categories))
}


fn plot_categories_over_length(
    lengths: &[f64],
    categories: &[FxHashMap<Category, f64>],
    region_names: &[String],
    region_ends: &[f64],
    filename: impl AsRef<Path>,
) {
    let l_max = *lengths.last().unwrap_or(&1.0);

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
        .margin_right(200)
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

        let between = matches!(category, Category::Between(_, _));

        if between {
            chart.draw_series(LineSeries::new(
            series,
            color.stroke_width(2),
            )).unwrap()
        } else {
            chart.draw_series(LineSeries::new(
                series,
                color.mix(0.8).stroke_width(1),  // ← thinner and more transparent
            )).unwrap()
        };
    } // end of for loop

    // manual legend
    let legend_x = 840_i32;
    let legend_y_start = 50_i32;
    let row_height = 22_i32;
    let font = ("sans-serif", 14).into_font();

    for (i, category) in all_categories.iter().enumerate() {
        let color = Palette99::pick(i).mix(0.9);
        let y = legend_y_start + (i as i32) * row_height;

        root.draw(&PathElement::new(
            vec![(legend_x, y + 6), (legend_x + 20, y + 6)],
            ShapeStyle { color: color.to_rgba(), filled: true, stroke_width: 2 },
        )).unwrap();

        let label = match category {
            Category::Within(a)     => format!("Within {}", region_name(region_names, *a)),
            Category::Between(a, b) => format!("Between {} and {}", region_name(region_names, *a), region_name(region_names, *b)),
            Category::WithRest(a)   => format!("{} with rest", region_name(region_names, *a)),
        };

        root.draw_text(
            &label,
            &font.clone().into_text_style(&root),
            (legend_x + 26, y),
        ).unwrap();
    }


    root.present().unwrap();
}




