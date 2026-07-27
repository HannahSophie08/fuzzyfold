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
} 

fn main() -> Result<()> {

    let cli = Cli::parse();
    let is_rna = false;
    let (header, sequence, _structure) = read_fasta_file(&cli.input, is_rna)?;

    let genome = sequence.to_string();

    let coordinates = read_coordinates(&cli.coordinates)?;

    let csv_path = cli.output.with_extension("csv");
    let mut file = File::create(&csv_path)
        .with_context(|| format!("failed to create output file {:?}", csv_path))?;

    if let Some(h) = &header {
        writeln!(file, "# {}", h)?;
    }

    let mut writer = csv::Writer::from_writer(file);

    writer.write_record(&["start", "end", "accessibilities"])?;

    for (start, end) in &coordinates {
        let seq = get_sequence(*start, *end, genome.clone())?;
        let accessibility: Vec<f64> = run_ff_accessibility("./target/release/ff-accessibility", &seq, cli.num_sims)
            .with_context(|| format!("simulation failed on coordinates {}-{}", start, end))?;

        let result = accessibility 
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(";");
        
        writer.write_record(&[start.to_string(), end.to_string(), result])?;
    }
    writer.flush()?;
    Ok(())
}

fn read_coordinates(path: &str) -> Result<Vec<(usize, usize)>> {

    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("failed to open coordinates read_fasta_file"))?;

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

fn run_ff_accessibility(binary_path: &str, sequence: &str, num_sims: usize) -> Result<Vec<f64>> {

    let mut child = Command::new(binary_path)
        .arg("-")
        .arg("--num-sims").arg(num_sims.to_string())
        .arg("--t-lin").arg("0")
        .arg("--t-log").arg("1")
        .arg("--t-ext").arg("0.02")
        .arg("--t-end").arg("1")
        .arg("--k0").arg("1e5")
        .arg("--dna")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .context("failed to spawn simulation binary")?;

    child.stdin
        .take()
        .unwrap()
        .write_all(sequence.as_bytes())?;

    let output = child.wait_with_output()?;
    if !output.status.success() {
        anyhow::bail!("Simulation failed: {}", String::from_utf8_lossy(&output.stderr))
    }
    let stdout = String::from_utf8(output.stdout)?;
    parse_accessibility_line(&stdout)
}

fn parse_accessibility_line(stdout: &str) -> Result<Vec<f64>> {

    let line = stdout
        .lines()
        .rev()
        .find(|l| !l.trim().is_empty())
        .context("ff-accessibility produced not output!")?;

    let start = line.find('[').context("no vector found in output")?;
    let end = line.find(']').context("vector not terminated!")?;

    line[start + 1..end]
        .split(',')
        .map(|v| v.trim().parse::<f64>().context("failed to parse accessibilities"))
        .collect()
}
