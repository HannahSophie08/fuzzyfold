# Single RNA cotranscriptional folding trajectories from stochastic simulations

## Overview

The program `ff-co_trajectory` simulates co-transcriptional RNA folding and produces trajectories that show how the structure changes over time. The program extends the sequence nucleotide by nucleotide and performs a stochastic kinetic folding simulation for each transcript length. Each trajectory is represented as a sequence of structures with associated energies and transition times.

---

## Using FASTA input

You can start a trajectory simulation from a predefined **FASTA file**, which
contains a sequence and optionally an initial structure. By giving an initial structure the start length for the cotranscriptional simulation can be specified. When no initial structure is given, the simulations start from transcript length 1. 

For example, the file `dld3.fa` contains a designed sequence, and as no structure is given the simulation would start at transcript length 1:

```fasta
>dld3
UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
```

Here is an example with a starting length of three from an unfolded structure: 

```fasta
>dld3
UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
...
```

## Input Parameters

You can specify the transcription rate via the extension time (`--t-ext`), which defines the simulation time per transcript length. The posttranscriptional folding time (`--t-end`) defines the simulation time at full sequence length. To express simulation time in seconds, set the rate constant `k0` to `1e6`. To further familiarize yourself with the parameters:

```bash
ff-co_trajectory --help
```

You can run the trajectory simulation as follows:

```bash
cat dld3.fa | ff-co_trajectory --k0 1e6 --t-ext 0.02 --t-end 1.0
```
All visited structures are printed with their energies and transition times.

---

## Pausing sites

When you do not want an uniform extension time, but define pausing sites, you can define the positions and the duration of the pauses in two separate input parameters: `p-pos` are the positions of the pausing sites and `t-pos` the respective duration. To have a pause at position 2 for 1s and one at position 10 for 0.5s:

```bash
cat dld3.fa | ff-co_trajectory --k0 1e6 --t-ext 0.02 --t-end 1.0 --p-pos 2,10 --t-pos 1,0.5
```

## Example output

An example simulation output is shown below for a full sequence length:

```
GCGUUUCCAGGGUUUAGACGGACGGGUGUGACUCGCCCAGCCCCGACCUC   energy   arrival-time   waiting-time    mean-waiting
..................................................     0.00   0.00000000e0  5.48838308e-1    1.15990269e0
.......................(.......)..................     2.10  5.48838308e-1  7.42186429e-1   6.00461698e-1
.......................(.......).....(......).....     4.50   1.29102474e0  6.62519714e-1   3.57034800e-1
.....................................(......).....     2.40   1.95354445e0  8.70308936e-3   5.15932279e-1
..................................................     0.00   1.96224754e0   2.24827455e0    1.15990269e0
..................(......)........................     1.90   4.21052209e0   1.08942365e0   4.14255646e-1
.................((......)).......................     0.10   5.29994574e0  3.93454939e-1   5.63963354e-1
...............(.((......)).).....................    -0.20   5.69340068e0  1.06769149e-1   4.52730587e-1
.............(.(.((......)).).)...................     1.80   5.80016983e0  4.33243095e-1   3.16296938e-1
.............(((.((......)).)))...................    -2.40   6.23341293e0   1.28813977e0   8.21431625e-1
...........(.(((.((......)).))))..................    -1.60   7.52155270e0  2.45415452e-1   2.40731917e-1
.........(.(.(((.((......)).)))).)................    -1.50   7.76696815e0  6.94202383e-2   4.10851300e-1
.........(((.(((.((......)).))))))................    -5.60   7.83638839e0   3.47571920e0    3.26611240e0
.........(((..((.((......)).)).)))................    -2.90   1.13121076e1   1.55874453e0   4.13337266e-1
.........(((.(((.((......)).))))))................    -5.60   1.28708521e1   2.36328447e0    3.26611240e0
.........(((.(((.((......)).))))))(.....).........    -2.10   1.52341366e1  7.41356896e-3   4.84982503e-1
.........(((.(((.((......)).))))))((...)).........    -4.70   1.52415502e1   2.00829892e1    1.23961105e1
..(......(((.(((.((......)).))))))((...))......)..    -0.50   3.53245394e1   1.26546682e0   7.99897165e-1
.........(((.(((.((......)).))))))((...)).........    -4.70   3.65900062e1   4.58067121e1    1.23961105e1
```

Each line represents one **structure** visited during the trajectory, with columns:

| Header | Description |
|---------|-------------|
| **sequence** | Corresponding structure in dot-bracket notation. |
| **energy** | Free energy evaluation (kcal/mol). |
| **arrival-time** | Simulation time at which the structure was observed. |
| **waiting-time** | Time spent in the structure until transition. |
| **mean-waiting** | Mean waiting time in the structure (1/flux) |

---

## See also

For cotranscriptional ensemble-level analysis across many trajectories, see the example in
[`ff-co_timecourse`](../ff-co_timecourse/README.md), which aggregates population data
over multiple stochastic runs.


