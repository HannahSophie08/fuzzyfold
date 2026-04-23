use clap::Args;
use anyhow::bail;
use anyhow::Result;
use ff_kinetics::Arrhenius;
use ff_structure::DotBracketVec;

#[derive(Debug, Args)]
pub struct RateModelArguments {
    /// Rate constant for add/delete moves.
    #[arg(long, default_value_t = 1e5)]
    pub k0: f64,

    /// Rate constant for three-way shift moves (optional, default = off).
    #[arg(long)]
    pub k3ws: Option<f64>,

    /// Rate constant for four-way shift moves (optional, default = off).
    #[arg(long)]
    pub k4ws: Option<f64>,
}

impl RateModelArguments {
    /// Validate that all parameters make sense.
    pub fn build_model(&self, celsius: f64) -> Arrhenius {
        Arrhenius::new(celsius, self.k0, self.k3ws, self.k4ws)
    }
}


#[derive(Debug, Args)]
pub struct TimelineParameters {
    /// The last time point of the linear scale.
    #[arg(long, default_value_t = 0.0)]
    pub t_ext: f64,

    /// Simulation stop time.
    #[arg(long, default_value_t = 1.0)]
    pub t_end: f64,

    /// Number of time points on the linear scale: [0..t-ext]
    #[arg(long, default_value_t = 100)]
    pub t_lin: usize,

    /// Number of time points on the logarithmic scale: [t-ext..t-end]
    #[arg(long, default_value_t = 100)]
    pub t_log: usize,

    /// Which mode? t_sep given: seperator between linear and logarithimic part at t_sep, 
    /// t-lin: timepoints recorded on a linear timescale between 0 and t-sep, t-log: timepoints 
    /// recorded on a logarithmic timescale between t-sep and total time.
    #[arg(long, default_value_t = -inf)]
    pub t_sep: usize,   
}

impl TimelineParameters {
    /// Validate that all parameters make sense.
    pub fn validate(&self) -> Result<()> {
        if self.t_end <= self.t_sep {
            bail!("t_end ({}) must be greater than t_sep ({})", self.t_end, self.t_sep);
        }
        if self.t_lin == 0 && self.t_log > 1 {
            bail!("t_lin must be > 0 if t_log > 1 (got t_lin={}, t_log={})", self.t_lin, self.t_log);
        }
        Ok(())
    }

    pub fn get_output_times(&self, sequence: String, structure: DotBracketVec) -> Vec<f64> {
        let t_end = self.t_end;
        let t_ext = self.t_ext;
        let t_lin = self.t_lin;
        let t_log = self.t_log;
        let t_sep = self.t_sep;
        let mut times = vec![0.0];

        if t_ext == 0.0 { // full length simulation
            if t_sep == -inf { // no t_sep given => t_lin and t_log are applied globally
                 // Linear time points: append `t_lin` evenly spaced points on t_end
                let start = *times.last().unwrap();
                let step = t_end / t_lin as f64;
                for i in 1..=t_lin {
                    times.push(start + i as f64 * step);
                }

                // Logarithmic time points 
                let start = 0.0;
                let log_start = start.ln();
                let log_end = t_end.ln();
                for i in 1..t_log {
                    let frac = i as f64 / t_log as f64;
                    let value = (log_start + frac * (log_end - log_start)).exp();
                    times.push(value);
                }
                times.push(t_end);
                times

            } else { 
            // t_sep given => t_lin evenly spaced timepoints on a linear timescale between 0 and t-sep and 
            // t_log timepoints on a logarithmic timescale 
                // Linear segments: append `t_lin` evenly spaced points
                let start = *times.last().unwrap();
                let step = t_ext / t_lin as f64;
                for i in 1..=t_lin {
                    times.push(start + i as f64 * step);
                }

                // Logarithmic tail
                let start = *times.last().unwrap();
                let log_start = start.ln();
                let log_end = t_end.ln();
                for i in 1..t_log {
                    let frac = i as f64 / t_log as f64;
                    let value = (log_start + frac * (log_end - log_start)).exp();
                    times.push(value);
                }
                times.push(t_end);

                times
            }

        } else {  // Co-transcriptional simulation 
            
            if t_sep == -inf { // no t_sep given => t_lin = timepoints per transcript length and t_log timepoints selected for full length simulation 
                // Co-transcriptional folding: linear timescale, posttranscriptional folding: logarithmic timescale 
                // Linear time points: append `t_lin` evenly spaced points on t_end
                let start = *times.last().unwrap();
                let step = t_ext/ t_lin as f64;
                for i in 1..=t_lin {
                    times.push(start + i as f64 * step);
                }

                // Logarithmic time points 
                let start = 0.0;
                let log_start = start.ln();
                let log_end = t_end.ln();
                for i in 1..t_log {
                    let frac = i as f64 / t_log as f64;
                    let value = (log_start + frac * (log_end - log_start)).exp();
                    times.push(value);
                }
                times.push(t_end);
                times

            } else { 
            // t_sep given => t_lin evenly spaced timepoints on a linear timescale between 0 and t-sep and 
            // t_log timepoints on a logarithmic timescale 
                // Linear segments: append `t_lin` evenly spaced points
                let start = *times.last().unwrap();
                let step = t_ext / t_lin as f64;
                for i in 1..= t_lin {
                    times.push(start + i as f64 * step);
                }

                // Logarithmic tail
                let start = *times.last().unwrap();
                let log_start = start.ln();
                let log_end = t_end.ln();
                for i in 1..t_log {
                    let frac = i as f64 / t_log as f64;
                    let value = (log_start + frac * (log_end - log_start)).exp();
                    times.push(value);
                }
                times.push(t_end);

                times
            }
            }

        // Linear segments: append `t_lin` evenly spaced points
        let start = *times.last().unwrap();
        let step = t_ext / t_lin as f64;
        for i in 1..=t_lin {
            times.push(start + i as f64 * step);
        }

        // Logarithmic tail
        let start = *times.last().unwrap();
        let log_start = start.ln();
        let log_end = t_end.ln();
        for i in 1..t_log {
            let frac = i as f64 / t_log as f64;
            let value = (log_start + frac * (log_end - log_start)).exp();
            times.push(value);
        }
        times.push(t_end);

        times
    }
}


