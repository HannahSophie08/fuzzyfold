use std::io::Write;
use std::process::{Command, Stdio};
use std::path::PathBuf;
use std::fs::File;

use anyhow::{Result, Context, bail};
use clap::Parser;
use csv::ReaderBuilder;

use fuzzyfold::input_parsers::read_fasta_file;
use ff_structure::DotBracket;
use ff_structure::DotBracketVec; 
use ff_structure::PairTable;

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
    num_samples: usize,

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

    let csv_path = cli.output.with_extension("csv");
    let mut file = File::create(&csv_path)
        .with_context(|| format!("failed to create output file {:?}", csv_path))?;

    if let Some(h) = &header {
        writeln!(file, "# {}", h)?;
    }

    let mut writer = csv::Writer::from_writer(file);

    writer.write_record(&["start", "end", "accessibilities"])?;

    let total = rows_to_process.len();

    for (i, (start, end)) in rows_to_process.iter().enumerate() {

        println!("{}/{}", i+1, total);
        let seq = get_sequence(*start, *end, genome.clone())?;
        let structures = run_rnasubopt(&seq, cli.num_samples)?;
        let accessibility: Vec<f64> = get_accessibility(structures);
        
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

fn run_rnasubopt(sequence: &str, num_samples: usize) -> Result<Vec<DotBracketVec>> {

    let mut child = Command::new("RNAsubopt")
        .arg("-p")
        .arg(num_samples.to_string())
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;

    child.stdin
        .take()
        .context("failed to open RNAsubopt stdin")?
        .write_all(sequence.as_bytes())?;

    let output = child.wait_with_output()?;

    if !output.status.success() {
        bail!("RNAsubopt failed: {}", String::from_utf8_lossy(&output.stderr));
    }
    
    let stdout = String::from_utf8(output.stdout)?;

    let structures = stdout
        .lines()
        .skip(2)
        .map(|l| {
            let dot_bracket = l.split_whitespace().next().unwrap();
            let pt = PairTable::try_from(dot_bracket).expect("invalid dot-bracket");
            DotBracketVec::from(&pt)
        })
        .collect();

    Ok(structures)
}

fn get_accessibility(structures: Vec<DotBracketVec>) -> Vec<f64> {

    let n = structures.len(); 

    let seq_len = structures.first().map(|s| s.0.len()).unwrap_or(0);

    let mut accessibility: Vec<usize> = vec![0; seq_len];

    for structure in structures {
        for (pos, state) in structure.iter().enumerate() {
            if matches!(state, DotBracket::Unpaired) {
                accessibility[pos] += 1;
            }
        }
    } 

    accessibility.iter().map(|&count| count as f64 / n as f64).collect()
}


