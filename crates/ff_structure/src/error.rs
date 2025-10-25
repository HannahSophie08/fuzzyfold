use std::fmt;

#[derive(Debug)]
pub enum StructureError {
    UnmatchedOpen(usize),          // '(' at this position was never closed
    UnmatchedMultiOpen((usize, usize)),          // '(' at this position was never closed
    UnmatchedClose(usize),         // ')' at this position has no matching '('
    UnmatchedMultiClose((usize, usize)),          // ')' at this position has no matching '('
    InvalidToken(String, String, usize),   // invalid char and position
    InvalidPairTable(usize),   // invalid char and position
}

impl fmt::Display for StructureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StructureError::UnmatchedOpen(i) => {
                write!(f, "Unmatched '(' at position {}", i)
            }
            StructureError::UnmatchedMultiOpen((si, di)) => {
                write!(f, "Unmatched '(' at strand {}, domain {}", si, di)
            }
            StructureError::UnmatchedClose(i) => {
                write!(f, "Unmatched ')' at position {}", i)
            }
            StructureError::UnmatchedMultiClose((si, di)) => {
                write!(f, "Unmatched ')' at strand {}, domain {}", si, di)
            }
            StructureError::InvalidToken(tok, src, i) => {
                write!(f, "Invalid {} in {} at position {}", tok, src, i)
            }
            StructureError::InvalidPairTable(i) => {
                write!(f, "Invalid entry at pair table position {}", i)
            }
        }
    }
}

impl std::error::Error for StructureError {}


