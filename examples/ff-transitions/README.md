Note: `ff-transitions` is an experimental, pre-publication software. It has to be
actively enabled in Cargo.toml

# Analyze macro-state transitions.
`ff-transitions` performs multiple stochastic folding simulations between all
macro-states, in order to derive the success-probability of a transition and the
distribution of arrival times.

## Input Macro-states

To partition the overall secondary-structure ensemble into smaller ensembles of
interest, we define **macro-states** using files such as `dld1_lm*.ms`.

Example (`dld1_lm3.ms`):

```fasta
>LM3 (delta = 3.00)
UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
.((((....)))).((((........))))...............
.((((....)))).((((.(....).))))...............
.((((....))))..(((........)))................
.((((....)))).((((.(.....)))))...............
.(((......))).((((........))))...............
..(((....)))..((((........))))...............
.(((......)))..(((........)))................
.(((.(...)))).((((........))))...............
```

Here:
- The first line defines the **macro-state name** (`LM3`) and, optionally, some more description after a white-space (`delta = 3.00`).
- The second line specifies the **sequence**.
- The remaining lines list all **secondary structures** that belong to this macro-state. Note: Information after the secondary
structures, like corresponding free energies, are discarded. All structures are re-evaluated upon construction of the macro-state.

---

## Simulation setup

Simulate 1000 trajectories starting from each macro-state:

```bash
ff-transitions --macrostates dld1*.ms --num-sims 1000 --output dld1
```

---

The program uses `STDOUT` to report simulation parameters, display a progress
bar, and summarize the simulation results in a table. The main output is a
chemical reaction network (`--output.crn`), where each transition is
represented by a coupled pair of reactions. If sufficient data from
trajectories has been collected, then this CRN exhibits the same dynamics as
the corresponding stochastic simulations.

---

## Aggregating data from multiple simulations

The file `--output.dat` contains the results of trajectories, and is
automatically loaded when ff-transitions is called with the same `--output`
parameter. Thus, it is possible to *accumulate results incrementally* by
reloading previous simulation results. 

For example:

```bash
ff-transitions --macrostates dld1*.ms --num-sims 9000 --output dld1
```

This command updates `dld1.dat`, to include the results from the additional
9000 simulations. Running the same command again will automatically reload the
file, add another 9000 simulations.

---

Currently, `fuzzyfold` does not include a CRN simulator; meanwhile, the output
format is compatible with the [crnsimulator](https://pypi.org/project/crnsimulator) python package, e.g., simulate
with:

```
cat dld1.crn | crnsimulator --t0 1e-5 --t8 100 --t-log 100 --pyplot dld1.pdf --force --labels LM1 LM2 LM3 LM4 LM5 --p0 LM3=1
```


## TODO:
 - Store trajectory data more efficiently.
 - Build CRNs from forward and reverse trajectories.
 - Report table of mean folding times.
 - Provide Rust-based CRN simulation.


