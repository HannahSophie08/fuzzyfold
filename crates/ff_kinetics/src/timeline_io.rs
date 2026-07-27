use std::fs;
use std::result;
use std::sync::Arc;
use serde::{Serialize, Deserialize};

use ff_energy::EnergyModel;

use crate::timeline::Timeline;
use crate::timeline::TimelineError;
use crate::macrostates::MacrostateRegistry;

#[derive(Serialize, Deserialize)]
pub struct SerializableTimeline {
    points: Vec<SerializableTimePoint>,
}

#[derive(Serialize, Deserialize)]
pub struct SerializableTimePoint {
    time: f64,
    ensemble: Vec<(usize, String, usize)>, // (length, macrostate name, count)
    counter: usize,
}

fn same_time(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

impl<E: EnergyModel> Timeline<E> {
    pub fn to_serializable(&self) -> SerializableTimeline {
        SerializableTimeline {
            points: self.points.iter().map(|tp| {
                let ensemble = tp.ensemble.iter()
                    .map(|(id, count)| {
                        let len = self.registry.macrostates()[*id].0;
                        let name = self.registry.macrostates()[*id].1.name().to_string();
                        (len, name, *count)
                    })
                    .collect();
                SerializableTimePoint {
                    time: tp.time,
                    ensemble,
                    counter: tp.counter,
                }
            }).collect()
        }
    }

    /// Load a timeline from a JSON file, checking against the provided registry
    pub fn from_file<P: AsRef<std::path::Path>>(
        path: P,
        times: &[f64],
        registry: Arc<MacrostateRegistry<E>>,
    ) -> result::Result<Self, TimelineError> {
        let data = fs::read_to_string(path)?;
        let serial: SerializableTimeline = serde_json::from_str(&data)?;

        // Sanity check: number of timepoints must match
        if serial.points.len() != times.len() {
            return Err(TimelineError::TimepointCountMismatch {
                found: serial.points.len(),
                expected: times.len(),
            });
        }

        let mut timeline = Timeline::new(times, Arc::clone(&registry));

        for (tp, serial_tp) in timeline.points.iter_mut().zip(serial.points) {
            if !same_time(tp.time, serial_tp.time) {
                return Err(TimelineError::TimeMismatch {
                    file_time: serial_tp.time,
                    expected_time: tp.time,
                });
            }

            for (length, name, count) in serial_tp.ensemble {
                // Look up macrostate by name in registry
                if let Some((idx, _m)) = registry.iter().find(|(_, (l, m))| *l == length && m.name() == name) {
                    *tp.ensemble.entry(idx).or_insert(0) += count;
                    tp.counter += count;
                } else {
                    return Err(TimelineError::MacrostateNotFound(name));
                }
            }
        }
        Ok(timeline)
    }
}