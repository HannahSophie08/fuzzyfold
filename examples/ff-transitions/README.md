# Analyze macrostate transitions.
`ff-transitions` performs multiple stochastic folding simulations between all
macrostates, in order to derive the success-probability of a transition and the
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
- The first line defines the **macro-state name** (`LM3`) and, optionally, some more description after a white-space (`lmin=lm3_bh=3.0`).
- The second line specifies the **sequence**.
- The remaining lines list all **secondary structures** that belong to this macro-state.

---

## Simulation setup

Simulate 1000 trajectories starting from each macrostate:

```bash
ff-transitions --macrostates dld1*.ms --num-sims 1000 --output dld1
```

---

The program uses `STDOUT` to report simulation parameters, display a **progress
bar**, and summarize the simulation results in tabular format. The main output
is a chemical reaction network (`dld1.crn`), where each transition is represented by
a coupled pair of reactions. If sufficient data from trajectories has been
collected, then this CRN is guaranteed to follow the same dynamics as
stochastic simulations.

---

## Aggregating data from multiple simulations

The file `dld1.dat` contains the results of trajectories, and is automatically
loaded when ff-transitions is called with the same `--output dld1` parameter.
Thus, it is possible to *accumulate results incrementally* by reloading
previous simulation results. 

For example:

```bash
ff-transitions --macrostates dld1*.ms --num-sims 9000 --output dld1
```

This command updates `dld1.dat`, to include the results from the additional
9000 simulations. Running the same command again will automatically reload the
file, add another 9000 simulations.

Currently, `fuzzyfold` does not include a CRN simulator. The output format is 
compatible with the `crnsimulator` python package, e.g. simulate with:

cat dld1.crn | crnsimulator --t0 1e-5 --t8 100 --t-log 100 --pyplot dld1.pdf --force --labels LM1 LM2 LM3 LM4 LM5 --p0 LM3=1


## TODO
 - Include data from reverse trajectories.
 - Make storage more efficient.
 - Report table of mean folding times.


