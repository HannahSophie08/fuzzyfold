use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use std::sync::Arc;
use std::path::PathBuf;
use rand::SeedableRng;
use rand::rngs::SmallRng;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::ViennaRNA;
use ff_kinetics::SSA;
use ff_kinetics::shift_policy;
use ff_kinetics::Arrhenius;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::MacrostateRegistry;
use ff_kinetics::timeline::Timeline;
use ff_energy::parameters::RNA_EXTENDED;
use ff_energy::parameters::RNA_TURNER_2004;
use ff_energy::parameters::DNA_MATHEWS_2004;
use ff_shared::kinetics_parsers::TimelineParameters;

//TODO: support shifts, rename to arrhenius

#[pyclass]
pub struct Simulator {
    energy_model: Arc<ViennaRNA>,
    rate_model: Arrhenius,
    is_rna: bool,
}


#[pymethods]
impl Simulator {
    #[new]
    #[pyo3(signature = (
        params = "rna_default",
        celsius=37.0,
        k0=1e5,
        k3ws=0.0,
        k4ws=0.0,
    ))]
    fn new(
        params: &str,
        celsius: f64,
        k0: f64,
        k3ws: f64,
        k4ws: f64,
    ) -> PyResult<Self> {
        let mut is_rna = true;
        let thermo = match params {
            "rna_default" => &RNA_TURNER_2004,
            "rna_extended" => &RNA_EXTENDED,
            "dna" => {
                is_rna = false;
                &DNA_MATHEWS_2004
            },
            _ => {
                return Err(PyValueError::new_err(
                    format!(
                        "Unknown parameter set '{}'. \
                         Valid options are: 'rna_default', 'rna_extended', 'dna'.",
                        params
                    )
                ));
            }
        };

        if k0 < 0.0 || k3ws < 0.0 || k4ws < 0.0 {
            return Err(PyValueError::new_err(
                "Rate constants must be non-negative",
            ));
        }

        let energy_model = ViennaRNA::from_thermo_params(thermo, celsius);
        let rate_model = Arrhenius::new(
            celsius,
            k0,
            Some(k3ws),
            Some(k4ws),
        );

        Ok(Self {
            energy_model: Arc::new(energy_model),
            rate_model,
            is_rna,
        })
    }

    #[pyo3(signature = (
            sequence,
            start=None,
            t_ext=None,
            t_end=1.0,
    ))]
    fn simulate(
        &self,
        sequence: &str,
        start: Option<&str>,
        t_ext: Option<f64>,
        t_end: f64,
    ) -> PyResult<SimulationIterator> {

        let (sequence, start_pt, times) = parse_inputs(self, sequence, start, t_ext, t_end)?;
        
        match (self.rate_model.k3ws().is_some(), self.rate_model.k4ws().is_some()) {
            (false, false) => build_iterator(
                sequence,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::NoShift,
                SSAKind::NoShift,
            ),

            (true, false) => build_iterator(
                sequence,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::ThreeWayOnly,
                SSAKind::ThreeWayOnly,
            ),

            (false, true) => build_iterator(
                sequence,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::FourWayOnly,
                SSAKind::FourWayOnly,
            ),

            (true, true) => build_iterator(
                sequence,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::ThreeAndFour,
                SSAKind::ThreeAndFour,
            ),
        }
   }
   
   #[pyo3(signature = (
            sequence,
            start=None,
            t_ext=None,
            t_end=1.0,
            t_lin=None,
            t_log=50,
            t_sep=None,
            num_sims=100,
    ))]
    fn simulate_timecourse(
        &self,
        py: Python<'_>,
        sequence: &str,
        start: Option<&str>,
        t_ext: Option<f64>,
        t_end: f64,
        t_lin: Option<usize>,
        t_log: usize,
        t_sep: Option<f64>,
        num_sims: usize,
    ) -> PyResult<Vec<(f64, FxHashMap<String, usize>)>> {

       let (sequence, start_pt, times) = parse_inputs(self, sequence, start, t_ext, t_end)?;

       let k3ws = self.rate_model.k3ws().is_some();
       let k4ws = self.rate_model.k4ws().is_some();
       let energy_model = Arc::clone(&self.energy_model);
       let rate_model = self.rate_model;

       let mut tl_params = TimelineParameters {
        t_ext,
        t_end,
        t_sep,
        t_lin,
        t_log,
       };

       let num_ext = sequence.len() - start_pt.len();
       let k0 = self.rate_model.k0().ok_or_else(|| PyValueError::new_err("rate model has no k0 set"))?;

       tl_params.validate(k0, num_ext).map_err(|e| PyValueError::new_err(e.to_string()))?;

       let output_times = tl_params.get_output_times(num_ext).map_err(|e| PyValueError::new_err(e.to_string()))?;

       let results: Result<Vec<Vec<String>>, String> = py.allow_threads (|| {
            let run_one = |_: usize| -> Result<Vec<String>, String> {
                let mut structures: Vec<String> = Vec::new();

                macro_rules! run_with_policy {
                    ($policy:expr) => {{
                        let walker = LoopNeighbors::try_from((
                            sequence.clone(), &start_pt, Arc::clone(&energy_model), $policy,
                        )).map_err(|e| format!("{:?}", e))?;

                        let mut ssa = SSA::from((walker, rate_model));
                        let mut rng = SmallRng::from_os_rng();
                        let mut t_idx = 0; 

                        ssa.co_simulate(&mut rng, &times, |t, tinc, _flux, w| {
                            while t_idx < output_times.len() && t + tinc >= output_times[t_idx] {
                                structures.push(w.to_string());
                                t_idx += 1;
                            }
                            true
                        });
                    }};

                }

                match (k3ws, k4ws) {
                    (false, false) => run_with_policy!(shift_policy::NoShift),
                    (true,  false) => run_with_policy!(shift_policy::ThreeWayOnly),
                    (false, true)  => run_with_policy!(shift_policy::FourWayOnly),
                    (true,  true)  => run_with_policy!(shift_policy::ThreeAndFour),
                }
                Ok(structures)
            };

            (0..num_sims).into_par_iter().map(run_one).collect()
       });
       
       let results = results.map_err(PyValueError::new_err)?;

       let mut counts: Vec<FxHashMap<String, usize>> = (0..output_times.len()).map(|_| FxHashMap::default()).collect();
       
       for structures in &results {
            for (t_idx, structure) in structures.iter().enumerate() {
                *counts[t_idx].entry(structure.clone()).or_insert(0) += 1;
            }    
        }

        Ok(output_times.iter().copied().zip(counts).collect())

    }

    #[pyo3(signature = (
            sequence,
            start=None,
            t_ext=None,
            t_end=1.0,
            t_lin=None,
            t_log=50,
            t_sep=None,
            num_sims=100,
            macrostates=vec![],
    ))]

    fn simulate_macrostates(
        &self,
        py: Python<'_>,
        sequence: &str,
        start: Option<&str>,
        t_ext: Option<f64>,
        t_end: f64,
        t_lin: Option<usize>,
        t_log: usize,
        t_sep: Option<f64>,
        num_sims: usize,
        macrostates: Vec<PathBuf>,
    ) -> PyResult<Vec<(f64, FxHashMap<String, f64>)>> {

       let (sequence, start_pt, times) = parse_inputs(self, sequence, start, t_ext, t_end)?;

       let k3ws = self.rate_model.k3ws().is_some();
       let k4ws = self.rate_model.k4ws().is_some();
       let energy_model = Arc::clone(&self.energy_model);
       let rate_model = self.rate_model;

       let mut tl_params = TimelineParameters {
        t_ext,
        t_end,
        t_sep,
        t_lin,
        t_log,
       };

       let num_ext = sequence.len() - start_pt.len();
       let k0 = self.rate_model.k0().ok_or_else(|| PyValueError::new_err("rate model has no k0 set"))?;

       tl_params.validate(k0, num_ext).map_err(|e| PyValueError::new_err(e.to_string()))?;

       let output_times = tl_params.get_output_times(num_ext).map_err(|e| PyValueError::new_err(e.to_string()))?;

       let mut ms = MacrostateRegistry::from((sequence.clone(), energy_model.clone()));
       ms.insert_files(&macrostates, t_ext.is_some()).map_err(|e| PyValueError::new_err(e.to_string()))?;
       let registry = Arc::new(ms);

       let timelines: Result<Vec<Timeline<ViennaRNA>>, String> = py.allow_threads (|| {
            let run_one = |_: usize| -> Result<Timeline<ViennaRNA>, String> {
                let thread_registry = Arc::clone(&registry);
                let mut timeline = Timeline::new(&output_times, thread_registry);

                macro_rules! run_with_policy {
                    ($policy:expr) => {{
                        let walker = LoopNeighbors::try_from((
                            sequence.clone(), &start_pt, Arc::clone(&energy_model), $policy,
                        )).map_err(|e| format!("{:?}", e))?;

                        let mut ssa = SSA::from((walker, rate_model));
                        let mut rng = SmallRng::from_os_rng();
                        let mut t_idx = 0; 

                        ssa.co_simulate(&mut rng, &times, |t, tinc, _flux, w| {
                            while t_idx < output_times.len() && t + tinc >= output_times[t_idx] {
                                let structure = w.current_structure();
                                timeline.assign_structure(t_idx, &structure);
                                t_idx += 1;
                            }
                            true
                        });
                    }};

                }

                match (k3ws, k4ws) {
                    (false, false) => run_with_policy!(shift_policy::NoShift),
                    (true,  false) => run_with_policy!(shift_policy::ThreeWayOnly),
                    (false, true)  => run_with_policy!(shift_policy::FourWayOnly),
                    (true,  true)  => run_with_policy!(shift_policy::ThreeAndFour),
                }
                Ok(timeline)
            };

            (0..num_sims).into_par_iter().map(run_one).collect()
       });

       let timelines = timelines.map_err(PyValueError::new_err)?;

       let mut master = Timeline::new(&output_times, Arc::clone(&registry));

       for timeline in timelines {
            master.merge(timeline);
       }

       Ok(timeline_to_occupancy(&master))  
        
    }
}

fn parse_inputs(
        sim: &Simulator,
        sequence: &str,
        start: Option<&str>,
        t_ext: Option<f64>,
        t_end: f64,
    ) -> PyResult<(NucleotideVec, PairTable, Vec<f64>)> {
        let sequence = match sim.is_rna {
            true => NucleotideVec::try_from_rna(sequence)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
            false => NucleotideVec::try_from_dna(sequence)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        };

        let start_db = match start {
            Some(s) => DotBracketVec::try_from(s)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
            None => DotBracketVec::try_from(".")
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        };

        if start_db.len() < sequence.len() && t_ext.is_none() {
            return Err(PyValueError::new_err(
                    "t_ext must be provided when start is shorter than sequence",
            ));
        }

        let times = if let Some(dt) = t_ext {
            let mut v = vec![dt; sequence.len() - start_db.len()];
            v.push(t_end);
            v
        } else {
            vec![t_end]
        };

        let start_pt = PairTable::try_from(&start_db)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok((sequence, start_pt, times))
}

fn timeline_to_occupancy(timeline: &Timeline<ViennaRNA>) -> Vec<(f64, FxHashMap<String, f64>)> {
    
    let macrostates = timeline.registry.macrostates();
    let mut name_order: Vec<&str> = Vec::new();
    let mut name_to_indices: FxHashMap<&str, Vec<usize>> = FxHashMap::default();
    for (idx, (_len, ms)) in macrostates.iter().enumerate() {
        let name = ms.name();
        name_to_indices
            .entry(name)
            .or_insert_with(|| {
                name_order.push(name);
                Vec::new()
            })
        .push(idx);
    }
    
    timeline.points.iter().map(|tp| {
        let mut occupancy: FxHashMap<String, f64> = FxHashMap::default();
        for name in &name_order {
            let indices = &name_to_indices[name];
            let count: usize = indices.iter()
                .map(|&i| tp.ensemble.get(&i).copied().unwrap_or(0))
                .sum();
            let occu = if tp.counter > 0 {
                count as f64 / tp.counter as f64
            } else { 0.0 };
            occupancy.insert(name.to_string(), occu);
        }
        (tp.time, occupancy)
    }).collect()
}

fn build_iterator<P>(
    seq: NucleotideVec,
    start_pt: &PairTable,
    energy_model: Arc<ViennaRNA>,
    rate_model: Arrhenius,
    times: Vec<f64>,
    policy: P,
    wrap: fn(SSA<LoopNeighbors<ViennaRNA, P>, Arrhenius>) -> SSAKind,
) -> PyResult<SimulationIterator>
where
    P: shift_policy::ShiftPolicy,
{
    let walker = LoopNeighbors::try_from((
        seq,
        start_pt,
        energy_model,
        policy,
    ))
    .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let ssa = wrap(SSA::from((walker, rate_model)));

    Ok(SimulationIterator {
        ssa,
        rng: SmallRng::from_os_rng(),
        times,
        elapsed: 0.0,
        finished: false,
    })
}

enum SSAKind {
    NoShift(SSA<LoopNeighbors<ViennaRNA, shift_policy::NoShift>, Arrhenius>),
    ThreeWayOnly(SSA<LoopNeighbors<ViennaRNA, shift_policy::ThreeWayOnly>, Arrhenius>),
    FourWayOnly(SSA<LoopNeighbors<ViennaRNA, shift_policy::FourWayOnly>, Arrhenius>),
    ThreeAndFour(SSA<LoopNeighbors<ViennaRNA, shift_policy::ThreeAndFour>, Arrhenius>),
}

#[pyclass]
pub struct SimulationIterator {
    ssa: SSAKind,
    rng: SmallRng,
    times: Vec<f64>,
    elapsed: f64,
    finished: bool,
}

#[pymethods]
impl SimulationIterator {

    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    fn __next__(
        mut slf: PyRefMut<Self>
    ) -> Option<(String, i32, f64, f64, f64)> {

        let this: &mut Self = &mut slf;

        if this.finished {
            return None;
        }

        let mut produced: Option<(String, i32, f64, f64, f64)> = None;

        let rng = &mut this.rng;
        let mut mytinc = 0.0;
        let mut first_pass = true;

        macro_rules! dispatch_ssa {
            ($ssa:expr) => {{
                $ssa.co_simulate(
                    rng,
                    &this.times,
                    |t, tinc, flux, w| {
                        if first_pass {
                            mytinc = tinc.min(this.times[0]);

                            produced = Some((
                                    w.to_string(),
                                    w.current_energy(),
                                    this.elapsed + t,
                                    mytinc,
                                    flux,
                            ));

                            this.elapsed += mytinc;
                            first_pass = false;
                            // advance the simulator to update the structure.
                            true
                        } else {
                            false
                        }
                    },
                    );
            }};
        }

        match &mut this.ssa {
            SSAKind::NoShift(ssa) => dispatch_ssa!(ssa),
            SSAKind::ThreeWayOnly(ssa) => dispatch_ssa!(ssa),
            SSAKind::FourWayOnly(ssa) => dispatch_ssa!(ssa),
            SSAKind::ThreeAndFour(ssa) => dispatch_ssa!(ssa),
        }

        if (this.times[0] - mytinc).abs() < f64::EPSILON {
            this.times.remove(0); 
            if this.times.is_empty() {
                this.finished = true;
            }
        } else {
            assert!(this.times[0] > mytinc);
            this.times[0] -= mytinc;
        }
        produced
    }
}


