//we need those
use ahash::AHashSet;
use ahash::AHashMap;
use ff_energy::NucleotideVec;
use rand::Rng;

use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_energy::EnergyModel;

use crate::LoopStructure;
use crate::Macrostate;
use crate::MacrostateRegistry;
use crate::RateModel;
use crate::reaction::ApplyMove;

fn find_neighbors<'a, E: EnergyModel, R: RateModel>(
    dbr: &DotBracketVec,
    lss_opt: Option<&LoopStructure<'a, E>>,
    sequence: &'a NucleotideVec,
    energy_model: &'a E,
    rate_model: &'a R,
    origin: &AHashMap<DotBracketVec, (i32, f64)>,
    visited: &mut AHashSet<DotBracketVec>,
    neighbors: &mut AHashMap<DotBracketVec, (i32, f64)>,
) {
    // Stop if already processed
    if !visited.insert(dbr.clone()) {
        return;
    }

    // So the outer loop does not produce 
    // all loop structures unnecessarily.
    let mut lss = match lss_opt {
        Some(existing) => existing.clone(),
        None => {
            let pairings = PairTable::try_from(dbr).unwrap();
            LoopStructure::try_from((&sequence[..], &pairings, energy_model)).unwrap()
        }
    };

    let mut mdbr = dbr.clone();
    for (bp_move, delta) in lss.all_moves() {
        // Let's look up whether the move is 
        // worth applying based on DotBracketVec.
        mdbr.apply_move(bp_move);

        if origin.contains_key(&mdbr) {
            lss.apply_move(bp_move);
            find_neighbors(&mdbr, Some(&lss), sequence, energy_model, rate_model, origin, visited, neighbors);
            lss.undo_move(bp_move);
        } else {
            // Rate to step out of the macrostate => P(i|alpha) * k_{i->j}
            let rate = origin.get(dbr).unwrap().1 * rate_model.rate(delta);
            neighbors
                .entry(mdbr.clone())
                .and_modify(|(e, k)| {
                    debug_assert!(*e == lss.energy() + delta);
                    *k += rate;
                })
            .or_insert((lss.energy() + delta, rate));
        }
        mdbr.undo_move(bp_move);
    }
}


#[derive(Debug)]
pub struct ExitMacrostate<'a> {
    parent_macrostate: &'a Macrostate,
    ensemble: AHashMap<DotBracketVec, (i32, f64)>,
    k_alpha: f64,
}

impl<'a, E: EnergyModel, R: RateModel> From<(&'a Macrostate, &'a NucleotideVec, &'a E, &'a R)> for ExitMacrostate<'a> {
    fn from((parent_macrostate, sequence, energy_model, rate_model): (&'a Macrostate, &'a NucleotideVec, &'a E, &'a R)) -> Self {
        let mut visited = AHashSet::default();
        let mut ensemble = AHashMap::default();
        for dbr in parent_macrostate.ensemble().keys() {
            find_neighbors::<E, R>(
                dbr,
                None,
                sequence,
                energy_model,
                rate_model,
                parent_macrostate.ensemble(),
                &mut visited,
                &mut ensemble,
            );
        }

        let mut k_alpha = 0.0;
        for (dbv, (en, k_ij)) in &ensemble {
            let pt = PairTable::try_from(dbv)
                .expect("Invalid dot-bracket for energy evaluation");
            debug_assert_eq!(*en, energy_model.energy_of_structure(sequence, &pt));
            k_alpha += k_ij;
        }

        ExitMacrostate {
            parent_macrostate,
            ensemble,
            k_alpha,
        }
    }
}

impl<'a> ExitMacrostate<'a> {
    pub fn parent_macrostate(&self) -> &Macrostate {
        self.parent_macrostate
    }

    pub fn ensemble(&self) -> &AHashMap<DotBracketVec, (i32, f64)> {
        &self.ensemble
    }

    pub fn k_alpha(&self) -> f64 {
        self.k_alpha
    }
 
    /// Number of secondary structures.
    pub fn len(&self) -> usize {
        self.ensemble.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ensemble.is_empty()
    }

    /// Check if a secondary structure is contained in this macrostate.
    pub fn contains(&self, structure: &DotBracketVec) -> bool {
        self.ensemble.contains_key(structure)
    }

    /// Randomly pick a structure according to the exit probability.
    pub fn get_random_microstate(&self) -> Option<DotBracketVec> {
        if self.ensemble.is_empty() {
            return None
        }
                // Draw a random number in [0, total)
        let mut rng = rand::rng();
        let mut t = rng.random_range(0.0..self.k_alpha);

        // Walk through ensemble and subtract until threshold crosses 0
        for (dbv, &(_, k_ij)) in &self.ensemble {
            t -= k_ij;
            if t <= 0.0 {
                return Some(dbv.clone());
            }
        }
        eprintln!("WARNING: rounding error observed. This should be rare!");
        self.ensemble.keys().next().cloned()
    }
}


pub struct ExitMacrostateRegistry<'a, E: EnergyModel, R: RateModel> {
    parent_registry: &'a MacrostateRegistry<'a, E>,
    rate_model: &'a R,
    /// By convention: macrostates[0] = unassigned.
    exit_macrostates: Vec<ExitMacrostate<'a>>,
}

impl<'a, E: EnergyModel, R: RateModel>
    From<(&'a MacrostateRegistry<'a, E>, &'a R)>
    for ExitMacrostateRegistry<'a, E, R>
{
    fn from((parent_registry, rate_model): (&'a MacrostateRegistry<'a, E>, &'a R)) -> Self {
        let mut exit_macrostates = Vec::with_capacity(parent_registry.len());

        // Index 0 is unassigned, so just an empty placeholder.
        exit_macrostates.push(ExitMacrostate {
            parent_macrostate: &parent_registry.macrostates()[0],
            ensemble: AHashMap::default(),
            k_alpha: 0.0,
        });

        // Compute neighbors for each real macrostate
        for (i, ms) in parent_registry.iter().skip(1) {
            eprintln!("Calculating neighbors for macrostate #{i}: {}", ms.name());
            let exit_ms = ExitMacrostate::from((
                ms, 
                parent_registry.sequence(),
                parent_registry.energy_model(),
                rate_model,
            ));
            exit_macrostates.push(exit_ms);
        }

        ExitMacrostateRegistry {
            parent_registry, 
            rate_model,
            exit_macrostates,
        }
    }
}

impl<'a, E: EnergyModel, R: RateModel> ExitMacrostateRegistry<'a, E, R> {
    pub fn parent_registry(&self) -> &MacrostateRegistry<'a, E> {
        self.parent_registry
    }

    pub fn rate_model(&self) -> &R {
        self.rate_model
    }

    pub fn exit_macrostates(&self) -> &Vec<ExitMacrostate<'a>> {
        &self.exit_macrostates
    }

    /// Number of exit_macrostates, including the catch-all unassigned macrostate.
    pub fn len(&self) -> usize {
        self.exit_macrostates.len()
    }

    //NOTE: Useless: there is always one.
    pub fn is_empty(&self) -> bool {
        self.exit_macrostates.is_empty()
    }

    /// Iterate over all macrostates
    pub fn iter(&self) -> impl Iterator<Item = (usize, &ExitMacrostate)> {
        self.exit_macrostates.iter().enumerate()
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;
    use crate::Metropolis;

    #[test]
    fn test_exit_macrostate_init() {
        /*        
        >lmin=lm3_bh=3.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((........))))...............
        .((((....)))).((((.(....).))))...............
        .((((....))))..(((........)))................
        .((((....)))).((((.(.....)))))...............
        .(((......))).((((........))))...............
        ..(((....)))..((((........))))...............
        .(((......)))..(((........)))................
        .(((.(...)))).((((........))))...............
        */        

        let energy_model = ViennaRNA::default();

        let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();
        let db1 = DotBracketVec::try_from(".((((....)))).((((........))))...............").unwrap();
        let db2 = DotBracketVec::try_from(".((((....)))).((((.(....).))))...............").unwrap();
        let db3 = DotBracketVec::try_from(".((((....))))..(((........)))................").unwrap();
        let db4 = DotBracketVec::try_from(".((((....)))).((((.(.....)))))...............").unwrap();
        let db5 = DotBracketVec::try_from(".(((......))).((((........))))...............").unwrap();
        let db6 = DotBracketVec::try_from("..(((....)))..((((........))))...............").unwrap();
        let db7 = DotBracketVec::try_from(".(((......)))..(((........)))................").unwrap();
        let db8 = DotBracketVec::try_from(".(((.(...)))).((((........))))...............").unwrap();

        let macrostate = Macrostate::from_list(
            "lmin=lm3_bh=3.0",
            &seq, 
            &[db1, db2, db3, db4, db5, db6, db7, db8], 
            &energy_model
        );

        let rate_model = Metropolis::new(energy_model.temperature(), 1.0);
        let neighbors = ExitMacrostate::from((
            &macrostate, 
            &seq, 
            &energy_model, 
            &rate_model
        ));
        println!("Neighbors '{}':", neighbors.parent_macrostate().name());
        println!("  Ensemble size: {}", neighbors.len());
        assert_eq!(neighbors.len(), 345);

        let ensemble = neighbors.ensemble().clone();
        let mut ensemble: Vec<_> = ensemble.iter().collect();
        ensemble.sort_by_key(|(_, (energy, _))| *energy);
        for (dbr, (energy, prob)) in ensemble.iter() {
            println!("  {} -> E(s) = {energy}, P(s|alpha) = {prob:.4}", dbr);
        }

    }

}


