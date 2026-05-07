//! # fuzzyfold
//!
//! The main entry point for the fuzzyfold nucleic acid folding package, with
//! interfaces to multiple crates that are organized as a workspace. At this
//! level, you can also find argument parsers for the binaries shipped with the
//! fuzzyfold package.
//!
//! This crate re-exports the main functionality from its submodules.


/// Exposing ff_structure as fuzzyfold::structure. A collection of practical data structures for
/// RNA structure representations. 
pub mod structure {
    pub use ::ff_structure::*;
}

/// Exposing ff_energy as fuzzyfold::energy. Handling of nucleotide sequences, nearest neighbor
/// loop decompositions and free energy evaluation models.
pub mod energy {
    pub use ::ff_energy::*;
}

/// Exposing ff_kinetics as fuzzyfold::kinetics. The main stochastic simulation framework,
/// introducing rate models, loop_structures, macrostates and time courses.
pub mod kinetics {
    pub use ::ff_kinetics::*;
}

/// Various flavors of handling sequence/structure input.
pub mod input_parsers;

/// Exposing the currently supported parameters of fuzzyfold's energy models.
pub mod energy_parsers;

/// Exposing the currently supported parameters of fuzzyfold's rate models and simulation parameters.
pub mod kinetics_parsers;

pub mod category;
