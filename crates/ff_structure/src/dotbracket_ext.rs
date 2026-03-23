use std::fmt;
use std::ops::Deref;
use std::ops::DerefMut;
use std::convert::TryFrom;

use crate::PairTable;
use crate::MultiPairTable;
use crate::MultiStruct;
use crate::StrandPairTable;
use crate::StructureError;

// DotbracketExt stores a vector representing structural constraints:
// '()': fixed paired
// 'x': fixed unpaired
// '.': unspecified => either paired or unpaired\
// '+': break
//
// Used for motifs, where contraints can be defined, but not the full sequence has to be specified
//


#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
/// Characters allowed for 'DotBracketExtVector'
pub enum DotBracketExt {
    Unpaired, // 'x' (forced unpaired)
    Open,     // '(' (forced paired)
    Close,    // ')' (forced paired)
    Unspecified, // '.' (either paired or unpaired)
    Break,    // '+' or '&'
}

impl TryFrom<char> for DotBracketExt {
    type Error = StructureError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'x' => Ok(DotBracketExt::Unpaired),
            '(' => Ok(DotBracketExt::Open),
            ')' => Ok(DotBracketExt::Close),
            '.' => Ok(DotBracketExt::Unspecified),
            '+' | '&' => Ok(DotBracketExt::Break),
            _ => Err(StructureError::InvalidToken(c.to_string(), "dot-bracket".into(), 0)),
        }
    }
}

impl From<DotBracketExt> for char {
    fn from(db: DotBracketExt) -> Self {
        match db {
            DotBracketExt::Open => '(',
            DotBracketExt::Close => ')',
            DotBracketExt::Unpaired => 'x',
            DotBracketExt::Unspecified => '.',
            DotBracketExt::Break => '+',
        }
    }
}

/// DotBracketExtVec is a compact representation of constraints for a secondary structure. 
/// Note that the field is public, to allow unsafe modifications. Thus, DotBracketExtVecs 
/// can be malformed and should be converted using the TryFrom trait.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct DotBracketExtVec(pub Vec<DotBracketExt>);

impl Deref for DotBracketExtVec {
    type Target = [DotBracketExt];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for DotBracketExtVec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl TryFrom<&str> for DotBracketExtVec {
    type Error = StructureError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut vec = Vec::with_capacity(s.len());
        for (i, c) in s.chars().enumerate() {
            match DotBracketExt::try_from(c) {
                Ok(db) => vec.push(db),
                Err(StructureError::InvalidToken(tok, src, _)) => {
                    return Err(StructureError::InvalidToken(tok, src, i));
                }
                Err(e) => return Err(e),
            }
        }
        Ok(DotBracketExtVec(vec))
    }
}


impl From<&PairTable> for DotBracketExtVec {
    fn from(pt: &PairTable) -> Self {
        let mut result: Vec<DotBracketExt> = Vec::with_capacity(pt.len());
        for (i, &j_opt) in pt.iter().enumerate() {
            match j_opt {
                None => result.push(DotBracketExt::Unpaired),
                Some(j) if (j as usize) > i => result.push(DotBracketExt::Open),
                Some(j) if (j as usize) < i => result.push(DotBracketExt::Close),
                Some(j) if (j as usize) == i => {
                    unreachable!("PairTable construction prevents self-pairing! ({})", i);
                }
                _ => unreachable!(),
            }
        }
        DotBracketExtVec(result)
    }
}

impl From<&MultiPairTable> for DotBracketExtVec {
    fn from(mpt: &MultiPairTable) -> Self {
        let mut db = Vec::with_capacity(mpt.len());

        for (i, entry) in mpt.iter().enumerate() {
            match entry {
                MultiStruct::Unpaired => {
                    db.push(DotBracketExt::Unpaired);
                }
                MultiStruct::StrandBreak => {
                    db.push(DotBracketExt::Break);
                }
                MultiStruct::Paired(j) => {
                    let j = *j as usize;
                    if i < j {
                        db.push(DotBracketExt::Open);
                    } else {
                        db.push(DotBracketExt::Close);
                    }
                }
            }
        }
        DotBracketExtVec(db)
    }
}

impl From<&StrandPairTable> for DotBracketExtVec {

    fn from(pt: &StrandPairTable) -> Self {
        let mut result: Vec<DotBracketExt> = Vec::with_capacity(pt.len() + pt.num_strands());

        for (si, strand) in pt.iter().enumerate() {
            for (di, &pair) in strand.iter().enumerate() {
                match pair {
                    None => result.push(DotBracketExt::Unpaired),
                    Some((sj, dj)) => {
                        let sj = sj as usize;
                        let dj = dj as usize;
                        if (sj, dj) > (si, di) {
                            result.push(DotBracketExt::Open);
                        } else if (sj, dj) < (si, di) {
                            result.push(DotBracketExt::Close);
                        } else {
                            panic!("Invalid self-pairing at strand {si}, domain {di}");
                        }
                    }
                }
            }
            // NOTE:: pushes strand break at the end of the DotBracketVec,
            // intentionally!
            result.push(DotBracketExt::Break);
        }

        DotBracketExtVec(result)
    }
}

impl fmt::Display for DotBracketExtVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for db in &self.0 {
            write!(f, "{}", char::from(*db))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dot_bracket_from_char() {
        assert_eq!(DotBracketExt::try_from('x').unwrap(), DotBracketExt::Unpaired);
        assert_eq!(DotBracketExt::try_from('(').unwrap(), DotBracketExt::Open);
        assert_eq!(DotBracketExt::try_from(')').unwrap(), DotBracketExt::Close);
        assert_eq!(DotBracketExt::try_from('.').unwrap(), DotBracketExt::Unspecified);
    }

    #[test]
    fn test_char_from_dot_bracket() {
        assert_eq!(char::from(DotBracketExt::Unpaired), 'x');
        assert_eq!(char::from(DotBracketExt::Open), '(');
        assert_eq!(char::from(DotBracketExt::Close), ')');
        assert_eq!(char::from(DotBracketExt::Unspecified), '.');
    }

    #[test]
    fn test_dot_bracket_from_invalid_char() {
        let res = DotBracketExt::try_from('y');
        assert!(matches!(res, Err(StructureError::InvalidToken(_, src, _)) if src == "dot-bracket"));
    }

    #[test]
    fn test_dot_bracket_vec_from_str() {
        let dbv = DotBracketExtVec::try_from("(x).").unwrap();
        assert_eq!(format!("{}", dbv), "(x).");
        assert_eq!(dbv.len(), 4);
        assert_eq!(dbv[0], DotBracketExt::Open);
        assert_eq!(dbv[1], DotBracketExt::Unpaired);
        assert_eq!(dbv[2], DotBracketExt::Close);
        assert_eq!(dbv[3], DotBracketExt::Unspecified);
    }

    #[test]
    fn test_dot_bracket_vec_from_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        let dbv = DotBracketExtVec::from(&pt);
        assert_eq!(format!("{}", dbv), "((xx))");
    }

    #[test]
    fn test_dot_bracket_vec_from_multi_pair_table_hack() {
        let pt = StrandPairTable::try_from("((..))+").unwrap();
        let dbv = DotBracketExtVec::from(&pt);
        assert_eq!(format!("{}", dbv), "((xx))+");
    }

    #[test]
    fn test_dot_bracket_vec_from_multi_pair_table() {
        let pt = StrandPairTable::try_from("((..)+)").unwrap();
        let dbv = DotBracketExtVec::from(&pt);
        assert_eq!(format!("{}", dbv), "((xx)+)+");
    }

}
