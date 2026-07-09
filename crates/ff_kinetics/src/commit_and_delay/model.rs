use serde::{Serialize, Deserialize};
use std::sync::Arc;
use std::path::Path;
use std::convert::TryFrom;
use rand::rng;
use ndarray::Array2;
use ndarray::s;
use ff_structure::PairTable;
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;


use crate::RateModel;
use crate::LoopStructure;
use crate::LoopStructureSSA;
use crate::commit_and_delay::ExitMacrostateRegistry;

type MacrostateID = usize;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactiveMicroTrajectory {
    i: DotBracketVec,
    j: DotBracketVec,
    simu_time: f64,
    mean_time: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactiveTrajectoryEnsemble {
    start: MacrostateID,
    stop: MacrostateID,
    t_min: Option<f64>,
    t_max: Option<f64>,
    successes: Vec<ReactiveMicroTrajectory>,
}

impl ReactiveTrajectoryEnsemble {
    pub fn sort_trajectories(&mut self) {
        self.successes.sort_by(|a, b| a.simu_time.partial_cmp(&b.simu_time).unwrap());
    }

    pub fn split_ensemble(&self, num_splits: usize) -> Vec<Self> {
        if self.successes.is_empty() || num_splits == 0 {
            return vec![];
        }

        let eps = 1.000001;

        // compute min/max simu_time
        let mint = self
            .successes
            .iter()
            .map(|t| t.mean_time)
            .fold(f64::INFINITY, f64::min);
        let maxt = self
            .successes
            .iter()
            .map(|t| t.mean_time)
            .fold(f64::NEG_INFINITY, f64::max);

        if mint <= 0.0 || maxt <= 0.0 {
            eprintln!("Warning: non-positive times found, cannot use log spacing.");
            return vec![self.clone()];
        }

        // logarithmic bin edges
        let log_min = (mint / eps).ln();
        let log_max = (maxt * eps).ln();
        let step = (log_max - log_min) / (num_splits + 1) as f64;

        let mut ensembles = Vec::with_capacity(num_splits + 1);
        for k in 0..=num_splits {
            let t_low = (log_min + k as f64 * step).exp();
            let t_high = (log_min + (k + 1) as f64 * step).exp();

            let chunk: Vec<_> = self
                .successes
                .iter()
                .cloned()
                .filter(|traj| traj.simu_time >= t_low && traj.simu_time < t_high)
                .collect();

            if !chunk.is_empty() {
                let t_min = chunk.iter().map(|t| t.simu_time).fold(f64::INFINITY, f64::min);
                let t_max = chunk.iter().map(|t| t.simu_time).fold(f64::NEG_INFINITY, f64::max);

                ensembles.push(Self {
                    start: self.start,
                    stop: self.stop,
                    t_min: Some(t_min),
                    t_max: Some(t_max),
                    successes: chunk,
                });
            }
        }

        ensembles
    }

    //pub fn split_ensemble(&self, num_splits: usize) -> Vec<Self> {
    //    if self.successes.is_empty() {
    //        return vec![];
    //    }

    //    //TODO: test with num_splits = 0
    //    let chunk_size = self.successes.len().div_ceil(num_splits);
    //    self.successes
    //        .chunks(chunk_size)
    //        .map(|chunk| Self {
    //            start: self.start,
    //            stop: self.stop,
    //            t_min: chunk.first().map(|t| t.simu_time),
    //            t_max: chunk.last().map(|t| t.simu_time),
    //            successes: chunk.to_vec(),
    //        })
    //        .collect()
    //}

    pub fn len(&self) -> usize {
        self.successes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.successes.is_empty()
    }

    pub fn mean_time_stats(&self) -> Option<(f64, f64, f64)> {
        if self.successes.is_empty() {
            return None;
        }

        let mut min_t = f64::INFINITY;
        let mut max_t = f64::NEG_INFINITY;
        let mut sum_t = 0.0;

        for traj in &self.successes {
            let t = traj.mean_time;
            if t < min_t {
                min_t = t;
            }
            if t > max_t {
                max_t = t;
            }
            sum_t += t;
        }

        let mean_t = sum_t / self.successes.len() as f64;
        Some((min_t, max_t, mean_t))
    }

}

impl<'a, E: EnergyModel, R: RateModel> CommitAndDelay<'a, E, R> {
    pub fn merge(&mut self, other: Self) {
        for ((i, j), val) in other.trajectories.indexed_iter() {
            if let Some(ens) = val {
                if let Some(ref mut existing) = self.trajectories[(i, j)] {
                    existing.successes.extend_from_slice(&ens.successes);
                } else {
                    self.trajectories[(i, j)] = Some(ens.clone());
                }
            }
        }
        

        // recompute marginals for self
        self.recompute_marginals();
    }
}

pub struct CommitAndDelay<'a, E: EnergyModel, R: RateModel> {
    exit_registry: Arc<ExitMacrostateRegistry<'a, E, R>>,
    trajectories: Array2<Option<ReactiveTrajectoryEnsemble>>,
}

impl<'a, E: EnergyModel, R: RateModel> From<Arc<ExitMacrostateRegistry<'a, E, R>>> 
for CommitAndDelay<'a, E, R> {
    fn from(exit_registry: Arc<ExitMacrostateRegistry<'a, E, R>>
    ) -> Self {
        let n = exit_registry.len();
        Self {
            exit_registry,
            trajectories: Array2::from_elem((n, n), None),
        }
    }
}

impl<'a, E: EnergyModel, R: RateModel> CommitAndDelay<'a, E, R> {
    /// recompute (0, j) and (i, 0) marginals
    pub fn recompute_marginals(&mut self) {
        let nrows = self.trajectories.nrows();
        let ncols = self.trajectories.ncols();

        // reset marginals
        for i in 0..nrows {
            self.trajectories[(i, 0)] = None;
        }
        for j in 0..ncols {
            self.trajectories[(0, j)] = None;
        }

        // sum outflows per row
        for i in 1..nrows {
            let mut merged = ReactiveTrajectoryEnsemble {
                start: i,
                stop: 0,
                t_min: None,
                t_max: None,
                successes: Vec::new(),
            };
            for j in 1..ncols {
                if let Some(ref ens) = self.trajectories[(i, j)] {
                    merged.successes.extend_from_slice(&ens.successes);
                }
            }
            if !merged.successes.is_empty() {
                self.trajectories[(i, 0)] = Some(merged);
            }
        }

        // sum inflows per column
        for j in 1..ncols {
            let mut merged = ReactiveTrajectoryEnsemble {
                start: 0,
                stop: j,
                t_min: None,
                t_max: None,
                successes: Vec::new(),
            };
            for i in 1..nrows {
                if let Some(ref ens) = self.trajectories[(i, j)] {
                    merged.successes.extend_from_slice(&ens.successes);
                }
            }
            if !merged.successes.is_empty() {
                self.trajectories[(0, j)] = Some(merged);
            }
        }
    }


    pub fn simulate_from(&mut self, start_id: MacrostateID) {
        let sequence = self.exit_registry.parent_registry().sequence();
        let energy_model = self.exit_registry.parent_registry().energy_model();
        let rate_model = self.exit_registry.rate_model();

        let start_ms = self.exit_registry.exit_macrostates()
            .get(start_id)
            .expect("invalid macrostate index");

        let start_db = start_ms.get_random_microstate().unwrap();
        let pairings = PairTable::try_from(&start_db).unwrap();
        let loops = LoopStructure::try_from((&sequence[..], &pairings, energy_model)).unwrap();
        let mut simulator = LoopStructureSSA::from((loops, rate_model));

        let mut mean_time = 0.0;
        simulator.simulate(
            &mut rng(), 
            f64::MAX,
            |t, _tinc, flux, ls| {
                let stop_db = DotBracketVec::from(ls);
                let stop_id = self.exit_registry.parent_registry().classify(&stop_db);
                //println!("current: {} {} {}", stop_db, stop_id, t);
                if stop_id != 0usize {
                    let traj = ReactiveMicroTrajectory {
                        i: start_db.clone(),
                        j: stop_db,
                        simu_time: t,
                        mean_time,
                    };

                    //Marginal
                    self.trajectories
                        .get_mut((start_id, 0))
                        .unwrap()
                        .get_or_insert_with(|| ReactiveTrajectoryEnsemble {
                            start: start_id,
                            stop: stop_id,
                            t_min: None,
                            t_max: None,
                            successes: Vec::new(),
                        })
                        .successes
                        .push(traj.clone());

                    //Marginal
                    self.trajectories
                        .get_mut((0, stop_id))
                        .unwrap()
                        .get_or_insert_with(|| ReactiveTrajectoryEnsemble {
                            start: start_id,
                            stop: stop_id,
                            t_min: None,
                            t_max: None,
                            successes: Vec::new(),
                        })
                        .successes
                        .push(traj.clone());

                    self.trajectories
                        .get_mut((start_id, stop_id))
                        .unwrap()
                        .get_or_insert_with(|| ReactiveTrajectoryEnsemble {
                            start: start_id,
                            stop: stop_id,
                            t_min: None,
                            t_max: None,
                            successes: Vec::new(),
                        })
                        .successes
                        .push(traj);
                    return false;
                }
                mean_time += 1.0/flux;
                true
            },
        );
    }

    pub fn simulate_between(&mut self, start_id: MacrostateID, stop_id: MacrostateID) {
        let sequence = self.exit_registry.parent_registry().sequence();
        let energy_model = self.exit_registry.parent_registry().energy_model();
        let rate_model = self.exit_registry.rate_model();

        let start_ms = self.exit_registry.exit_macrostates()
            .get(start_id)
            .expect("invalid macrostate index");

        let start_db = start_ms.get_random_microstate().unwrap();
        let pairings = PairTable::try_from(&start_db).unwrap();
        let loops = LoopStructure::try_from((&sequence[..], &pairings, energy_model)).unwrap();
        let mut simulator = LoopStructureSSA::from((loops, rate_model));

        let mut mean_time = 0.0;
        let mut curr_db = start_db;
        let mut curr_id = start_id;
        let mut toggle = 0;
        simulator.simulate(
            &mut rng(), 
            f64::MAX,
            |t, _tinc, flux, ls| {
                let next_db = DotBracketVec::from(ls);
                let next_id = self.exit_registry.parent_registry().classify(&next_db);
                println!("current: {} {} {}", next_db, next_id, t);
                //TODO: think about this more..
                if next_id != toggle {
                    if next_id != 0 {
                        let traj = ReactiveMicroTrajectory {
                            i: curr_db.clone(),
                            j: next_db.clone(),
                            simu_time: t,
                            mean_time,
                        };
                        self.trajectories
                            .get_mut((curr_id, next_id))
                            .unwrap()
                            .get_or_insert_with(|| ReactiveTrajectoryEnsemble {
                                start: curr_id,
                                stop: next_id,
                                t_min: None,
                                t_max: None,
                                successes: Vec::new(),
                            })
                        .successes
                            .push(traj);
                        if next_id == stop_id {
                            return false;
                        } 
                        curr_db = next_db;
                        curr_id = next_id;
                        mean_time = 0.0;
                    }
                    toggle = next_id;

                }
                mean_time += 1.0/flux;
                true
            },
        );
    }

    pub fn gather_data(&self) {
        let split = 0;
        let mut total = 0;
        for ((i, j), value) in self.trajectories.indexed_iter() {
            if i == 0 {
                continue;
            } 
            if j == 0 {
                total = value.as_ref().unwrap().len();
                continue;
            }
            if let Some(ens) = value {
                let k_alpha = self.exit_registry.exit_macrostates()[i].k_alpha();
                for (sid, sms) in ens.split_ensemble(split).iter().enumerate().map(|(i, x)| (i + 1, x)) {
                    println!("M_{i} -> S{sid}_[{i},{j}] @ {:14.8e} /s", k_alpha * sms.len() as f64 / total as f64);
                    println!("S{sid}_[{i},{j}] -> M_{j} @ {:14.8e} /s", 1f64/sms.mean_time_stats().unwrap().2);

                }
                println!("M_{i} -> T_[{i},{j}] @ {:14.8e} /s", k_alpha * ens.len() as f64 / total as f64);
                println!("T_[{i},{j}] -> M_{j} @ {:14.8e} /s", 1f64/ens.mean_time_stats().unwrap().2);
            } else {
                println!("M_{i} -> M_{j} @ 0");
            }
        }
    }

    pub fn to_rate_matrix(&self) -> Array2<f64> {
        todo!("let's do only k_commit for now")
    }

    pub fn trajectories(&self) -> &Array2<Option<ReactiveTrajectoryEnsemble>> {
        &self.trajectories
    }
}

#[derive(Serialize, Deserialize)]
pub struct SerializableCommitAndDelay {
    trajectories: Array2<Option<ReactiveTrajectoryEnsemble>>,
}

impl<'a, E: EnergyModel, R: RateModel> CommitAndDelay<'a, E, R> {
    pub fn to_serializable(&self) -> SerializableCommitAndDelay {
        let view = self.trajectories.slice(s![1.., 1..]);
        SerializableCommitAndDelay {
            trajectories: view.to_owned(),
        }
    }

    pub fn from_serializable(
        serial: SerializableCommitAndDelay,
        exit_registry: Arc<ExitMacrostateRegistry<'a, E, R>>,
    ) -> Self {
        let nrows = serial.trajectories.nrows() + 1;
        let ncols = serial.trajectories.ncols() + 1;

        // Allocate full array (including margins)
        let mut full = Array2::from_elem((nrows, ncols), None);

        // Copy serialized values into the sub-matrix (1.., 1..)
        full.slice_mut(s![1.., 1..]).assign(&serial.trajectories);
        

        // Leave row/col 0 empty or recompute marginals later
        Self {
            exit_registry,
            trajectories: full,
        }
    }
    pub fn save_json(&self, path: &str) -> anyhow::Result<()> {
        let serial = self.to_serializable();
        serde_json::to_writer_pretty(std::fs::File::create(path)?, &serial)?;
        Ok(())
    }

    pub fn load_json<P: AsRef<Path>>(
        path: P,
        exit_registry: Arc<ExitMacrostateRegistry<'a, E, R>>,
    ) -> anyhow::Result<Self> {
        let reader = std::io::BufReader::new(std::fs::File::open(path.as_ref())?);
        let serial: SerializableCommitAndDelay = serde_json::from_reader(reader)?;

        Ok(Self::from_serializable(serial, exit_registry))
    }

    
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;
    use crate::Metropolis;
    use crate::macrostates::MacrostateRegistry;

    fn test_ms1() -> std::io::Cursor<&'static [u8]> {
        Cursor::new(b">lmin=lm1_bh=4.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....((((.((((........))))..))))....)))).
        ..(((....((((.((((........))))..))))....)))..
        .((((....((((.((((.(....).))))..))))....)))).
        .(.((....((((.((((........))))..))))....)).).
        .((((....((((.((((.(.....)))))..))))....)))).
        .((.(....((((.((((........))))..))))....).)).
        .((((....(((..((((........))))...)))....)))).
        .((((....((((..(((........)))...))))....)))).
        .((((....((((.(((..........)))..))))....)))).
        .((((.....(((.((((........))))..))).....)))).
        .(((.....((((.((((........))))..)))).....))).
        .(((.....((((..(((........)))...)))).....))).
        ...((....((((.((((........))))..))))....))...
        .((((....((((.((((........))).).))))....)))).
        ..((.....((((.((((........))))..)))).....))..
        .((((....((((.((((........)))).).)))....)))).
        .(((.....((((.((((.(.....)))))..)))).....))).
        .((......((((.((((........))))..))))......)).
        .(((......(((.((((........))))..)))......))).
        .(((.....((((.((((........))).).)))).....))).")
    }

    fn test_ms2() -> std::io::Cursor<&'static [u8]> {
        Cursor::new(b">lmin=lm2_bh=4.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((.((((...))))..)))).........
        .((((....)))).((((.(((.....)))..)))).........
        ..(((....)))..((((.((((...))))..)))).........
        .((((....)))).(((..((((...))))...))).........
        .((((....))))..(((.((((...))))..)))..........
        .((((....)))).((((..(((...)))...)))).........
        .(((......))).((((.((((...))))..)))).........
        .((((....))))..(((.(((.....)))..)))..........
        ..(((....)))..((((.(((.....)))..)))).........
        .(((.(...)))).((((.((((...))))..)))).........
        .((((....))))..(((..(((...)))...)))..........
        .(((......))).((((.(((.....)))..)))).........
        ...((....))...((((.((((...))))..)))).........
        .(((......))).((((..(((...)))...)))).........
        .((((....)))).((((.((((...))).).)))).........
        .((((....)))).((((.((((...))))..))).)........
        ..(((....)))...(((.((((...))))..)))..........
        .(((......)))..(((.((((...))))..)))..........
        .((((....)))).((((...((...))....)))).........
        ..((......))..((((.((((...))))..)))).........
        .((((....)))).((((..((.....))...)))).........")
    }
 
    fn test_ms3() -> std::io::Cursor<&'static [u8]> {
        Cursor::new(b">lmin=lm3_bh=3.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((........))))...............
        .((((....)))).((((.(....).))))...............
        .((((....))))..(((........)))................
        .((((....)))).((((.(.....)))))...............
        .(((......))).((((........))))...............
        ..(((....)))..((((........))))...............
        .(((......)))..(((........)))................
        .(((.(...)))).((((........))))...............")
    }

    #[test]
    fn test_commit_and_delay_minimal() {
        let energy_model = ViennaRNA::default();
        let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();
        let mut registry = MacrostateRegistry::from((&seq, &energy_model));

        registry.insert_from_reader(test_ms1(), "manual").unwrap();
        assert_eq!(registry.len(), 2);

        let rate_model = Metropolis::new(energy_model.temperature(), 1.0);
        let exitreg = ExitMacrostateRegistry::from((&registry, &rate_model));

        let mut cad = CommitAndDelay::from(Arc::new(exitreg));
        cad.simulate_from(1);
        assert_eq!(cad.trajectories.get((1, 1)).and_then(|opt| opt.as_ref()).unwrap().len(), 1);
        cad.simulate_from(1);
        assert_eq!(cad.trajectories.get((1, 1)).and_then(|opt| opt.as_ref()).unwrap().len(), 2);
        cad.simulate_from(1);
        assert_eq!(cad.trajectories.get((1, 1)).and_then(|opt| opt.as_ref()).unwrap().len(), 3);
    }

    #[test]
    fn test_commit_and_delay() {
        let energy_model = ViennaRNA::default();
        let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();
        let mut registry = MacrostateRegistry::from((&seq, &energy_model));

        registry.insert_from_reader(test_ms1(), "manual").unwrap();
        registry.insert_from_reader(test_ms2(), "manual").unwrap();
        registry.insert_from_reader(test_ms3(), "manual").unwrap();
        assert_eq!(registry.len(), 4);

        let rate_model = Metropolis::new(energy_model.temperature(), 1.0);
        let _exitreg = ExitMacrostateRegistry::from((&registry, &rate_model));

        //NOTE: too slow for a unittest at the moment.
        // let mut cad = CommitAndDelay::from(Arc::new(exitreg));
        // cad.simulate_between(3,1);
        // cad.simulate_between(1,3);
        // cad.simulate_between(2,1);
        // cad.simulate_between(1,2);
        // cad.simulate_between(3,2);
        // cad.simulate_between(2,3);
        // for row in cad.trajectories.rows() {
        //     let line = row
        //         .iter()
        //         .map(|el| el.as_ref().map_or("0".into(), |ens| ens.len().to_string()))
        //         .collect::<Vec<_>>()
        //         .join(" ");
        //     println!("{line}");
        // }
    }
}

