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
use ff_structure::ConstrPosMap;
use ff_structure::PairSet;
use ff_structure::ConstrPos;
use ff_structure::NAIDX;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;


use crate::{K0, KB};

type NonPairSet = IntSet<usize>;

#[derive(Debug)]
pub struct Motif {
    name: String,
    constr_pos_map: ConstrPosMap,
    motif_energy: Option<f64>,
    motif_prob: Option<f64>,
    allowed_distance: NAIDX
}

impl Motif {
    /// Initialize a catch-all motif. 
    /// (This is the default motif when initializing a MotifRegistry.)
    pub fn new_catch_all(name: &str) -> Self {
        Motif { 
            name: name.to_owned(),
            constr_pos_map: ConstrPosMap::new(),
            motif_energy: None,
            motif_prob: None,
            allowed_distance: 0
        }
    }

    // Version 1: For `&DotBracketVec`
    pub fn from_list_dotbracket<'a, E: EnergyModel>(
        name: &str,
        sequence: &'a NucleotideVec,
        structure: &DotBracketVec,
        energy_model: &'a E,
        allowed_distance: NAIDX,
    ) -> Self {
        let rt = KB * (K0 + energy_model.temperature());

        Self {
            name: name.to_owned(),
            constr_pos_map: ConstrPosMap::try_from(structure).unwrap(),
            motif_energy: Some(0.0), // Change this to actual energy calculation
            motif_prob: Some(0.0),   // Change this to actual probability calculation
            allowed_distance: allowed_distance,
        }
    }

    // Version 2: For `&str`
    pub fn from_list_str<'a, E: EnergyModel>(
        name: &str,
        sequence: &'a NucleotideVec,
        structure: &str,
        energy_model: &'a E,
        allowed_distance: NAIDX,
    ) -> Self {
        let rt = KB * (K0 + energy_model.temperature());

        Self {
            name: name.to_owned(),
            constr_pos_map: ConstrPosMap::try_from(structure).unwrap(),
            motif_energy: Some(0.0), // Change this to actual energy calculation
            motif_prob: Some(0.0),   // Change this to actual probability calculation
            allowed_distance: allowed_distance,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn energy(&self) -> Option<f64> {
        self.motif_energy
    }

    pub fn prob(&self) -> Option<f64> {
        self.motif_prob
    }
        /// Check if a secondary structure is contained in this motif.
    pub fn contains(&self, structure: &PairTable) -> bool {
        for (key, value) in &self.constr_pos_map.0 {
            let entry: Option<NAIDX> = structure.get(key);

            let mut dist_counter = 0;
            match value {
                // Pair is mismatched -> add distance of 2
                ConstrPos::Pair(expected) => {
                    if entry != Some(*expected) {
                        dist_counter += 2;
                    }
                }
                // Unpaired Position is paired -> add distance of 1
                ConstrPos::X => {
                    if entry.is_some() {
                        dist_counter += 1;
                    }
                }
            }

            // If distance exceeds the allowed threshold, return false
            if dist_counter > self.allowed_distance {
                return false;
            }
        }
        true
    }

}



/// A registy to collect macrostate definitions.
pub struct MotifRegistry<'a, E: EnergyModel> {
    sequence: &'a NucleotideVec,
    energy_model: &'a E,
    /// By convention: macrostates[0] = unassigned.
    motifs: Vec<Motif>,
}

impl<'a, E: EnergyModel> From<(&'a NucleotideVec, &'a E)> for MotifRegistry<'a, E> {
    fn from((sequence, energy_model): (&'a NucleotideVec, &'a E)) -> Self {
        let motifs = vec![Motif::new_catch_all("Unassigned")];

        MotifRegistry {
            sequence,
            energy_model,
            motifs,
        }
    }
}

fn io_err(msg: &str, src: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("{} in {}", msg, src))
}

impl<'a, E: EnergyModel> MotifRegistry<'a, E> {

    /// High-level entry: read one or more macrostate files from disk.
    pub fn insert_from_file(&mut self, file: &PathBuf) -> io::Result<()> {
        let fh = File::open(file)?;
        let reader = BufReader::new(fh);
        self.insert_from_reader(reader, &file.display().to_string())?;
        Ok(())
    }

    pub fn insert_from_reader<R: BufRead>(&mut self, reader: R, source: &str) -> io::Result<()> {
        let mut lines = reader.lines();

        // Step 1: Read the sequence line first (before motifs)
        let seq_line = lines
            .next()
            .ok_or_else(|| io_err("Missing sequence line", source))??
            .trim()
            .to_string();

        // Step 2: Parse the sequence
        let file_seq = NucleotideVec::try_from(seq_line.as_str())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if &file_seq != self.sequence {
            return Err(io_err("Sequence does not match input sequence", source));
        }

        // Step 3: Initialize the motifs collection
        let mut motifs = Vec::new();

        // Step 4: Read the motif lines
        let mut warned_trailing = false;
        for (lineno, line) in lines.enumerate() {
            let line = line?;
            let line = line.trim();

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            // Split the motif line into the header (name, distance) and the dot-bracket structure
            let (motif_header, structure_str) = match line.split_once(char::is_whitespace) {
                Some((header, structure)) => (header, structure),
                None => {
                    return Err(io_err(
                        &format!("Invalid motif format at line {}", lineno + 3),
                        source,
                    ));
                }
            };

            // Step 5: Parse the motif header (name and distance)
            let mut header_parts = motif_header.split_whitespace();
            let motif_name = header_parts
                .next()
                .ok_or_else(|| io_err("Missing motif name", source))?
                .to_string();
            
            let distance_str = header_parts
                .next()
                .ok_or_else(|| io_err("Missing motif distance", source))?;

            // Attempt to parse the distance string as `NAIDX` (u16)
            let allowed_distance = distance_str
                .parse::<NAIDX>()
                .map_err(|_| io_err(&format!("Invalid motif distance: '{}'. Must be a valid number.", distance_str), source))?;

            // Step 6: Create the motif from the dot-bracket structure
            if structure_str.is_empty() {
                return Err(io_err(
                    &format!("Empty structure for motif at line {}", lineno + 3),
                    source,
                ));
            }

            // Step 7: Create the Motif
            let motif = Motif::from_list_str(
                &motif_name,
                &self.sequence,
                structure_str,
                self.energy_model,
                allowed_distance,
            );

            motifs.push(motif);

            // Step 8: Warn about trailing data (if any)
            if !warned_trailing {
                eprintln!(
                    "Warning: trailing data after dot-bracket structures is ignored in {}.",
                    source
                );
                warned_trailing = true;
            }
        }

        // Step 9: If no motifs are found, return an error
        if motifs.is_empty() {
            return Err(io_err("No motifs found", source));
        }

        // Step 10: Add the motifs to the macrostates
        self.motifs.extend(motifs);
        Ok(())
    }


    /// Try to classify a structure:
    /// - Returns Some(index) if exactly one macrostate contains the structure
    /// - Returns None if no macrostate matches
    /// - Panics if more than one macrostate matches (unimplemented)
    pub fn classify(&self, structure: &DotBracketVec) -> Vec<usize> {
        let mut matches: Vec<usize> = Vec::new();
        let structure_pt = PairTable::try_from(structure).unwrap();


        for (i, ms) in self.motifs.iter().enumerate() {
            if ms.contains(&structure_pt) {
                matches.push(i);
            }
        }

        match matches.len() {
            0 => vec![0usize],
            _ => matches
        }
    }

    pub fn sequence(&self) -> &NucleotideVec {
        self.sequence
    }

    pub fn energy_model(&self) -> &E {
        self.energy_model
    }

    pub fn macrostates(&self) -> &Vec<Motif> {
        &self.motifs
    }

    /// Number of macrostates, including the catch-all unassigned macrostate.
    pub fn len(&self) -> usize {
        self.motifs.len()
    }

    //NOTE: Useless: there is always one.
    pub fn is_empty(&self) -> bool {
        self.motifs.is_empty()
    }

    /// Iterate over all macrostates
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Motif)> {
        self.motifs.iter().enumerate()
    }


}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use ff_energy::{parameters::RNA_TURNER_2004, ViennaRNA};

    #[test]
    fn test_motif_init() {
        /*        
        >lmin=lm3_bh=3.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((........))))...............
        */        

        let energy_model = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0);
        let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();
        let db1 = DotBracketVec::try_from(".((((....)))).((((........))))...............").unwrap(); // REFERENCE
        let db2 = DotBracketVec::try_from(".((((....)))).((((........))))...............").unwrap(); // MATCHING
        let db3 = DotBracketVec::try_from(".((((....)))).((((........))))...............").unwrap(); // NON-MATCHING

        let dbs0 = "UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC";
        let dbs1 = ".((((xxxx)))).((((xxxxxxxx))))..............."; // REFERENCE
        let dbs2 = ".((((....)))).((((........))))...(......)...."; // MATCHING
        let dbs3 = ".((((....)))).(((((......)))))..............."; // NON-MATCHING


        // // Dot-Bracket Vec Based

        // let motif = Motif::from_list_dotbracket(
        //     "lmin=lm3_bh=3.0",
        //     &seq,
        //     &db1,
        //     &energy_model,
        //     0
        // );

        // println!("USING VECTOR");
        // println!("Motif '{}':", motif.name());
        // println!("E(motif) = {:.4}, P(motif) = {:.4}", motif.energy().unwrap(), motif.prob().unwrap());

        // assert_eq!(motif.energy(), Some(0.0));
        // assert_eq!(motif.prob(), Some(0.0));

        // //Correct Structure
        // let pt = PairTable::try_from(&db2).unwrap();

        // assert!(motif.contains(&pt));

        // //Incorrect Structure db3
        // let pt = PairTable::try_from(&db3).unwrap();

        // assert!(!motif.contains(&pt));



        // Dot-Bracket String Based
        let motif = Motif::from_list_str(
            "lmin=lm3_bh=3.0",
            &seq,
            &dbs1,
            &energy_model,
            0
        );

        println!("USING STRING");
        println!("Motif '{}':", motif.name());
        println!("E(motif) = {:.4}, P(motif) = {:.4}", motif.energy().unwrap(), motif.prob().unwrap());
        println!();


        //Correct Structure db2
        println!("{}", dbs0);
        println!("{}", dbs1);
        println!("{}", dbs2);
        println!();

        let pt = PairTable::try_from(dbs2).unwrap();
        for (k, v) in &motif.constr_pos_map.to_sorted_list() {
            println!("Pos: {} -> Motif: {:?}    Observed: {:?} ", k, v, pt.get(k));
        }

        assert!(motif.contains(&pt));


        //Correct Structure db3
        println!("{}", dbs0);
        println!("{}", dbs1);
        println!("{}", dbs3);
        println!();

        let pt = PairTable::try_from(dbs3).unwrap();
        for (k, v) in &motif.constr_pos_map.to_sorted_list() {
            println!("Pos: {} -> Motif: {:?}    Observed: {:?} ", k, v, pt.get(k));
        }

        assert!(!motif.contains(&pt));

    }

    #[test]
    fn test_motifregistry_init_and_classify() {
        let energy_model = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0);
        let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();

        let mut registry = MotifRegistry::from((&seq, &energy_model));

        // By convention: one unassigned macrostate
        assert_eq!(registry.len(), 1);
        assert_eq!(registry.macrostates()[0].name(), "Unassigned");

        let input = b"UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC\n
        >motif_name1 0\n
        .((((xxxx)))).((((xxxxxxxx))))...............\n
        >motif_name2 0\n
        .((((....)))).((((........))))...............\n
        >motif_name3 0\n
        .((((xxxx)))).((((........))))xxxxxxxxxxxxxxx\n";

        // registry.insert_from_reader(Cursor::new(input), "manual").unwrap();
        let file_path = PathBuf::from("/home/mescalin/dguerguerian/example_data/test_motif_file.fasta");
        registry.insert_from_file(&file_path)?;
        assert_eq!(registry.len(), 4);

        // Build a test macrostate with a few structures
        let s1 = DotBracketVec::try_from(".((((....)))).((((........)))).....(.....)...").unwrap();
        assert_eq!(registry.classify(&s1), vec![1]);
        let s2 = DotBracketVec::try_from(".((((....)))).((((((....))))))...............").unwrap();
        assert_eq!(registry.classify(&s2), vec![2]);
        let s3 = DotBracketVec::try_from(".((((....)))).((((........))))...............").unwrap();
        assert_eq!(registry.classify(&s3), vec![1, 2, 3]);

        // Unknown structure: should return 0 ("Unassigned")
        let s4 = DotBracketVec::try_from("..............").unwrap();
        assert_eq!(registry.classify(&s4), vec![0]);

        // Iteration test
        let all_names: Vec<_> = registry.iter().map(|(_, ms)| ms.name().to_string()).collect();
        assert!(all_names.contains(&"Unassigned".to_string()));
        assert!(all_names.contains(&"motif_name1".to_string()));
        assert!(all_names.contains(&"motif_name2".to_string()));
        assert!(all_names.contains(&"motif_name3".to_string()));
    }





}

