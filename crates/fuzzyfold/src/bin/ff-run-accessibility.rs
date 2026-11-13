use std::env;
use std::io::Write;
use std::process::{Command, Stdio};
use std::path::PathBuf;
use std::fs::File;

use anyhow::{Result, Context, bail};
use clap::Parser;
use csv::ReaderBuilder;

use fuzzyfold::input_parsers::read_fasta_file;

#[derive(Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

     /// Coordinates file 
    #[arg(long, value_name = "COORDINATES")]
    coordinates: String,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(long)]
    t_lin: Option<usize>,

    #[arg(long)]
    t_log: usize,

    #[arg(long)]
    t_ext: f64,

    #[arg(long)]
    t_end: f64,

    #[arg(long)]
    k0: f64,

    #[arg(long)]
    k3ws: Option<f64>,

    #[arg(long)]
    row: Option<usize>,
} 

fn main() -> Result<()> {

    let cli = Cli::parse();
    let is_rna = false;
    let (header, sequence, _structure) = read_fasta_file(&cli.input, is_rna)?;
    let genome = sequence.to_string();

    let rows_to_process: Vec<(usize, usize)> = match cli.row {
        Some(idx) => vec![read_coordinate_row(&cli.coordinates, idx)?],
        None => read_coordinates(&cli.coordinates)?,
    };

    let accessibility_path = find_ff_accessibility("ff-accessibility")?;

    let csv_path = cli.output.with_extension("csv");
    let mut file = File::create(&csv_path)
        .with_context(|| format!("failed to create output file {:?}", csv_path))?;

    if let Some(h) = &header {
        writeln!(file, "# {}", h)?;
    }

    let mut writer = csv::Writer::from_writer(file);

    writer.write_record(&["start", "end", "time", "length", "accessibilities"])?;

    let total = rows_to_process.len();

    for (i, (start, end)) in rows_to_process.iter().enumerate() {

        println!("{}/{}", i+1, total);
        let seq = get_sequence(*start, *end, genome.clone())?;

        let points = run_ff_accessibility(&accessibility_path, &seq, cli.num_sims, cli.t_lin, cli.t_log, cli.t_ext, cli.t_end, cli.k0, cli.k3ws)
            .with_context(|| format!("simulation failed on coordinates {}-{}", start, end))?;

        for (time, accessibility) in points {
            let result = accessibility
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(";");

            writer.write_record(&[start.to_string(), end.to_string(), time.to_string(),
                                  accessibility.len().to_string(), result])?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn find_ff_accessibility(name: &str) ->Result<PathBuf> {

    let exe_path = env::current_exe().context("failed to determine current executable path!")?;
    let dir = exe_path.parent().context("executable has no parent directory!")?;

    Ok(dir.join(name))

}

fn read_coordinates(path: &str) -> Result<Vec<(usize, usize)>> {

    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("failed to open coordinates file"))?;

    rdr.records()
        .map(|record| {
            let record = record.context("failed to read coordinate row")?;
            let start: usize = record.get(1)
                .context("missing starting column!")?
                .parse()
                .context("failed to parse start")?;
            let end: usize = record.get(2)
                .context("missing end column")?
                .parse()
                .context("failed to parse end")?;
            Ok((start, end))
        })
        .collect()
}

fn read_coordinate_row(path: &str, idx: usize) -> Result<(usize, usize)> {
    
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("failed to open coordinates file"))?;

    let record = rdr.records()
        .nth(idx)
        .with_context(|| format!("row index {} out of range", idx))?
        .with_context(|| format!("failed to read coordinate row {}", idx))?;
    let start: usize = record.get(1)
        .context("missing starting column!")?
        .parse()
        .context("failed to parse start")?;
    let end: usize = record.get(2)
        .context("missing end column")?
        .parse()
        .context("failed to parse end")?;

    Ok((start, end))
}


fn get_sequence(start: usize, end: usize, genome: String) -> Result<String> {

    if start > end {
        bail!("start ({}) must be smaller than end ({})", start, end);
    }

    if end > genome.len() {
        bail!("end ({}) is out of bounds of the genome", end);
    }

    let seq_up = "GTGACTGTGGAGATGAGGATCACCCATCT";
    let seq_dn = "AGA"; 

    Ok(format!("{}{}{}", seq_up, &genome[start..end], seq_dn))
}

fn run_ff_accessibility(binary_path: &PathBuf, sequence: &str, num_sims: usize, t_lin: Option<usize>, t_log: usize, t_ext: f64, t_end: f64, k0: f64, k3ws: Option<f64>) -> Result<Vec<(f64, Vec<f64>)>> {

    let mut cmd = Command::new(binary_path);
    cmd.arg("-")
        .arg("--num-sims").arg(num_sims.to_string())
        .arg("--t-log").arg(t_log.to_string())
        .arg("--t-ext").arg(t_ext.to_string())
        .arg("--t-end").arg(t_end.to_string())
        .arg("--k0").arg(k0.to_string())
        .arg("--dna")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    if let Some(t_lin) = t_lin {
        cmd.arg("--t-lin").arg(t_lin.to_string());
    }
    if let Some(k3ws) = k3ws {
        cmd.arg("--k3ws").arg(k3ws.to_string());
    }
    let mut child = cmd.spawn().context("failed to spawn simulation binary")?;

    child.stdin
        .take()
        .unwrap()
        .write_all(sequence.as_bytes())?;

    let output = child.wait_with_output()?;
    if !output.status.success() {
        anyhow::bail!("Simulation failed: {}", String::from_utf8_lossy(&output.stderr))
    }
    let stdout = String::from_utf8(output.stdout)?;
    parse_accessibility_lines(&stdout)
}

fn parse_accessibility_lines(stdout: &str) -> Result<Vec<(f64, Vec<f64>)>> {
    let mut points = Vec::new();
    for line in stdout.lines() {
        let Some(open) = line.find('[') else { continue };
        // Only lines that start with a time; skips the header and parameter lines.
        let Ok(time) = line[..open].trim().parse::<f64>() else { continue };
        let close = line.find(']').context("vector not terminated!")?;
        let inner = line[open + 1..close].trim();
        let values = if inner.is_empty() {
            Vec::new()
        } else {
            inner.split(',')
                .map(|v| v.trim().parse::<f64>().context("failed to parse accessibilities"))
                .collect::<Result<Vec<_>>>()?
        };
        points.push((time, values));
    }
    if points.is_empty() {
        bail!("ff-accessibility produced no time points");
    }
    Ok(points)
}