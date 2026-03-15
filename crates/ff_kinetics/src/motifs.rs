use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::collections::HashSet;
use std::io;
use std::path::PathBuf;
use ahash::AHashMap;
use rand::Rng;
use nohash_hasher::IntSet;

use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_structure::PairList;
use ff_structure::PairSet;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;


use crate::{K0, KB};

type NonPairSet = IntSet<usize>;

#[derive(Debug)]
pub struct Motif {
    name: String,
    forced_paired: PairSet,
    forced_unpaired: NonPairSet,
    motif_energy: Option<f64>,
    motif_prob: Option<f64>,
    dist_percentage: f64
}

impl Motif {
    /// Initialize a catch-all motif. 
    /// (This is the default motif when initializing a MotifRegistry.)
    pub fn new_catch_all(name: &str) -> Self {
        Motif { 
            name: name.to_owned(),
            forced_paired: PairSet::new(0),
            forced_unpaired: NonPairSet::default(),
            motif_energy: None,
            motif_prob: None,
            dist_percentage: 1.0
        }
    }

    pub fn from_list<'a, E: EnergyModel>(
        name: &str, 
        sequence: &'a NucleotideVec,
        structure: &DotBracketVec, 
        energy_model: &'a E,
        dist_percentage: f64
    ) -> Self {
        let pt = PairTable::try_from(structure)
                .expect("Invalid dot-bracket for energy evaluation");
        let ps = PairSet::from(&pt);
        let nps = structure.get_non_pair_set();

        let rt = KB * (K0 + energy_model.temperature());

        //Would need something like: energy_model.energy_of_constrained_structure(sequence, &pt)
                
        Self {
        name: name.to_owned(),
        forced_paired: ps,
        forced_unpaired: nps,
        motif_energy: Some(0.0), //Change this
        motif_prob: Some(0.0), //Change this
        dist_percentage: dist_percentage
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn n_paired(&self) -> usize {
        self.forced_paired.len()
    }

    pub fn n_unpaired(&self) -> usize {
        self.forced_unpaired.len()
    }

    pub fn energy(&self) -> Option<f64> {
        self.motif_energy
    }

    pub fn prob(&self) -> Option<f64> {
        self.motif_prob
    }

    /// Check if a secondary structure is contained in this motif.
    pub fn contains(&self, structure_paired: &PairSet, structure_unpaired: &NonPairSet) -> bool {
        // 1) Quick Check by lengths of sets
        // 1.1) Paired
        let n_paired = self.n_paired();
        let paired_diff = n_paired.abs_diff(structure_paired.len());
        let paired_threshold = (n_paired as f64 * self.dist_percentage).round() as usize;
        if paired_diff > paired_threshold {
            return false
        }
        // 1.2) Unpaired
        let n_unpaired = self.n_unpaired();
        let unpaired_diff = n_unpaired.abs_diff(structure_unpaired.len());
        let unpaired_threshold = (n_unpaired as f64 * self.dist_percentage).round() as usize;
        if unpaired_diff > unpaired_threshold {
            return false
        }
        // 2) Thorough Check by comparing length of intersect
        // 2.1) Paired
        let paired_intersect = self.forced_paired.intersect(&structure_paired);
        if n_paired - paired_intersect.len() > paired_threshold {
            return false
        }
        // 2.1) Unpaired
        let unpaired_intersect = self.forced_unpaired.intersection(&structure_unpaired).count();
        if n_unpaired - unpaired_intersect > unpaired_threshold {
            return false
        }
        true
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use ff_energy::{parameters::RNA_TURNER_2004, ViennaRNA};

    #[test]
    fn test_macrostatepl_init() {
        /*        
        >lmin=lm3_bh=3.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((........))))...............
        */        

        let energy_model = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0);
        let seq = NucleotideVec::try_from("UCAGUCUUGCUGCGCUGUAUCGAUACGGUUUCAGUUUUUAUUGC").unwrap();
        let db1 = DotBracketVec::try_from(".((((...)))).((((((....))))))...............").unwrap();
        let db2 = DotBracketVec::try_from(".((((...))))..(((((....)))))....(.....).....").unwrap(); // One Pair missing and one extra pair, should be within 0.1 distance
        let db3 = DotBracketVec::try_from(".((((...))))..(((((....)))))....((...)).....").unwrap(); // One Pair missing and two extra pair, should NOT be within 0.1 distance


        let motif = Motif::from_list(
            "lmin=lm3_bh=3.0",
            &seq,
            &db1,
            &energy_model,
            0.1
        );

        println!("Motif '{}':", motif.name());
        println!("n_paired={}, n_unpaired={}, E(motif) = {:.4}, P(motif) = {:.4}", motif.n_paired(), motif.n_unpaired(), motif.energy().unwrap(), motif.prob().unwrap());

        assert_eq!(motif.n_paired(), 10);
        assert_eq!(motif.n_unpaired(), 24);
        assert_eq!(motif.energy(), Some(0.0));
        assert_eq!(motif.prob(), Some(0.0));

        //Correct Structure
        let pt = PairTable::try_from(&db2).unwrap();
        let ps = PairSet::from(&pt);
        let nps = db2.get_non_pair_set();

        assert!(motif.contains(&ps, &nps));

        //Incorrect Structure
        let pt = PairTable::try_from(&db3).unwrap();
        let ps = PairSet::from(&pt);
        let nps = db3.get_non_pair_set();

        assert!(motif.contains(&ps, &nps));

    }


}

