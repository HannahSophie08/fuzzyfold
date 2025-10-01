use clap::Args;
use anyhow::Result;
use anyhow::bail;

#[derive(Debug, Args)]
pub struct RateModelParams {
    /// Metropolis rate constant (must be > 0).
    #[arg(long, default_value_t = 1e6)]
    pub k0: f64,
}

#[derive(Debug, Args)]
pub struct TimelineParameters {
    /// The last time point of the linear scale.
    #[arg(long, default_value_t = 1e-5)]
    pub t_ext: f64,

    /// Simulation stop time.
    #[arg(long, default_value_t = 1.0)]
    pub t_end: f64,

    /// Number of time points on the linear scale: [0..t-ext]
    #[arg(long, default_value_t = 1)]
    pub t_lin: usize,

    /// Number of time points on the logarithmic scale: [t-ext..t-end]
    #[arg(long, default_value_t = 20)]
    pub t_log: usize,
}

impl TimelineParameters {
    /// Validate that all parameters make sense.
    pub fn validate(&self) -> Result<()> {
        if self.t_end <= self.t_ext {
            bail!("t_end ({}) must be greater than t_ext ({})", self.t_end, self.t_ext);
        }
        if self.t_lin == 0 && self.t_log > 1 {
            bail!("t_lin must be > 0 if t_log > 1 (got t_lin={}, t_log={})", self.t_lin, self.t_log);
        }
        Ok(())
    }

    pub fn get_output_times(&self) -> Vec<f64> {
        let t_end = self.t_end;
        let t_ext = self.t_ext;
        let t_lin = self.t_lin;
        let t_log = self.t_log;
        let mut times = vec![0.0];

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


