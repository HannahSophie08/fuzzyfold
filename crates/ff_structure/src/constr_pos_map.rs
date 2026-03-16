//! PairTable construction and helper traits.

use std::ops::{Deref, DerefMut};
use std::convert::TryFrom;
use crate::NAIDX;
use crate::StructureError;
use std::collections::HashMap;
use crate::{DotBracket, DotBracketVec};

/// As of v0.1.3 the PairTable field is private. A pair-table should
/// be constructed by From or TryFrom traits, but then be save to use.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConstrPos {
    Pair(usize),
    X,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstrPosMap(pub HashMap<usize, ConstrPos>);

impl ConstrPosMap {
    pub fn new() -> Self {
        ConstrPosMap(HashMap::new())
    }
}

impl Deref for ConstrPosMap {
    type Target = HashMap<usize, ConstrPos>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for ConstrPosMap {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl TryFrom<&str> for ConstrPosMap {
    type Error = StructureError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut stack = Vec::new();
        let mut map = HashMap::new();

        for (i, c) in s.chars().enumerate() {
            match c {
                '(' => stack.push(i),

                ')' => {
                    let j = stack.pop().ok_or(StructureError::UnmatchedClose(i))?;
                    map.insert(j, ConstrPos::Pair(i));
                }

                'x' => {
                    map.insert(i, ConstrPos::X);
                }

                '.' => {}

                _ => {
                    return Err(StructureError::InvalidToken(
                        format!("character '{}'", c),
                        "structure".to_string(),
                        i,
                    ));
                }
            }
        }

        if let Some(i) = stack.pop() {
            return Err(StructureError::UnmatchedOpen(i));
        }

        Ok(ConstrPosMap(map))
    }
}

//NEEDS NEW CONSTRAINED DotBracketVec
impl TryFrom<&DotBracketVec> for ConstrPosMap {
    type Error = StructureError;

    fn try_from(db: &DotBracketVec) -> Result<Self, Self::Error> {
        let mut stack = Vec::new();
        let mut map = HashMap::new();

        for (i, dot) in db.iter().enumerate() {
            match dot {
                DotBracket::Open => stack.push(i),
                DotBracket::Close => {
                    let j = stack.pop().ok_or(StructureError::UnmatchedClose(i))?;
                    map.insert(i, ConstrPos::Pair(j));
                }
                // Make sure this enum is correct for forced unpaired positions
                DotBracket::Unpaired => {
                    map.insert(i, ConstrPos::X);
                }
                DotBracket::Break => unreachable!("unexpected Break in single-stranded case"),
            }
        }

        if let Some(i) = stack.pop() {
            return Err(StructureError::UnmatchedOpen(i));
        }

        Ok(ConstrPosMap(map))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_pair_table() {
        let cpm: ConstrPosMap = ConstrPosMap::try_from("((xxx)).(xxx)").unwrap();

        for (k, v) in &*cpm {
            println!("{} -> {:?}", k, v);
        }

        assert_eq!(cpm.len(), 9);

        assert_eq!(cpm.get(&0), Some(&ConstrPos::Pair(6)));
        assert_eq!(cpm.get(&1), Some(&ConstrPos::Pair(5)));

        assert_eq!(cpm.get(&2), Some(&ConstrPos::X));
        assert_eq!(cpm.get(&3), Some(&ConstrPos::X));
        assert_eq!(cpm.get(&4), Some(&ConstrPos::X));

        assert_eq!(cpm.get(&5), None); // ) doesn't have to be checked
        assert_eq!(cpm.get(&6), None); // ) doesn't have to be checked

        assert_eq!(cpm.get(&7), None); // .

        assert_eq!(cpm.get(&8), Some(&ConstrPos::Pair(12)));

        assert_eq!(cpm.get(&9), Some(&ConstrPos::X));
        assert_eq!(cpm.get(&10), Some(&ConstrPos::X));
        assert_eq!(cpm.get(&11), Some(&ConstrPos::X));

        assert_eq!(cpm.get(&12), None); 
    }
}