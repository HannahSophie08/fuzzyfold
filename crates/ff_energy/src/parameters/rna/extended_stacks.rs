//! # Pseudouridine Parameters
//!
//! **Authors:** Graham A. Hudson, Richard J. Bloomingdale, and Brent M. Znosko  
//! **Title:** *Thermodynamic contribution and nearest-neighbor parameters of
//! pseudouridine–adenosine base pairs in oligoribonucleotides*  
//! **Journal:** RNA 19:1474–1482  
//! **Year:** 2013  
//! **DOI:** 10.1261/rna.039610.113
//! 
//! # Inosine Parameters 
//!   
//! **Authors:** Daniel J. Wright, Jamie L. Rice, Dawn M. Yanker, and Brent M. Znosko
//! **Title:** Nearest Neighbor Parameters for Inosine-Uridine Pairs in RNA Duplexes
//! **Journal:** Biochemistry, 46, 4625-4634
//! **Year:** 2007 
//! **DOI:** 10.1021/bi0616910
//!  
//! **Authors:** Daniel J. Wright, Christopher R. Force and Brent M. Znosko
//! **Title:** Stability of RNA duplexes containing inosine-cytosine pairs
//! **Journal:** Nucleic Acids Research, Vol. 46, No. 22 12099-12108
//! **Year:** 2018
//! **DOI:** 10.1093/nar/gky907
//!
//! ## Description
//!
//! This module provides stacking parameters including **P** for
//! pseudouridine (Ψ) and **I** for insosine.
//!
//! ## Implementation Notes
//!
//! - All unspecified pseudouridine interactions are treated as **U**.
//! - [AP][AP] = [AP][AU]
//! - [AP][GU] = [AU][GU]
//! - [GP] = [GU] (not even shown in table.)
//! - All unspecified insosine interactions are treated as **G***.
//! - [IU][IU] = [IU][IG]
//! - [IU][GU] = [GU][GU]
//! - [AP][IU] = [AP][GU]
//! - [IU][AP] = [IU][AU]

      

use crate::parameters::parameterset::ExtendedStackParams;
use crate::parameters::parameterset::E;

pub static STACKPARAMS_EN37: [[i32; E]; E] = [
    /* [cl] [ri]:  AU     UA     CG     GC     GU     UG     AP     PA     IU     UI     CI     IC*/
    /* [AU] */ [ -110,   -90,  -210,  -220,  -140,   -60,  -280,  -274,   -50,   -41,  -157,  -102],                  
    /* [UA] */ [  -90,  -130,  -210,  -240,  -130,  -100,  -162,  -210,    43,    37,   -96,  -118],
    /* [CG] */ [ -210,  -210,  -240,  -330,  -210,  -140,  -277,  -220,  -122,   -77,  -186,  -223],
    /* [GC] */ [ -220,  -240,  -330,  -340,  -250,  -150,  -329,  -249,  -103,  -134,  -262,  -189],
    /* [GU] */ [ -140,  -130,  -210,  -250,   130,   -50,  -140,  -130,   130,   -50,  -210,  -250],
    /* [UG] */ [  -60,  -100,  -140,  -150,   -50,    30,   -60,  -100,   -50,    30,  -140,  -150],
    /* [AP] */ [ -280,  -162,  -277,  -329,  -140,   -60,  -280,  -162,  -140,   -60,  -277,  -329],
    /* [PA] */ [ -274,  -210,  -220,  -249,  -130,  -100,  -274,  -210,  -130,  -100,  -220,  -249],
    /* [IU] */ [  -50,    43,  -122,  -103,   130,   -50,   -50,    43,   358,   266,  -122,  -103],           
    /* [UI] */ [  -41,    37,   -77,  -134,   -50,    30,   -41,    37,   266,   223,   -77,  -134],
    /* [CI] */ [ -157,   -96,  -186,  -262,  -210,  -140,  -157,   -96,  -210,  -140,  -186,  -262],           
    /* [IC] */ [ -102,  -118,  -223,  -189,  -250,  -150,  -102,  -118,  -250,  -150,  -223,  -189],
];

pub static STACKPARAMS_ENTH: [[i32; E]; E] = [
    /* [cl] [ri]:  AU     UA     CG     GC     GU     UG     AP     PA     IU     UI     CI     IC */
    /* [AU] */ [ -940,  -680, -1050, -1140,  -880,  -320, -2208, -2694, -1583, -1168, -1420,  -770],
    /* [UA] */ [ -680,  -770, -1040, -1240, -1280,  -700, -2081, -1247,  -822, -1008, -1180, -1530],
    /* [CG] */ [-1050, -1040, -1060, -1340, -1210,  -560, -1623, -1119, -1338, -1199, -1270, -1450],
    /* [GC] */ [-1140, -1240, -1340, -1490, -1260,  -830, -2407, -1729, -1156,  -981, -1680, -1060],
    /* [GU] */ [ -880, -1280, -1210, -1260, -1460, -1350,  -880, -1280, -1460, -1350, -1210, -1260],
    /* [UG] */ [ -320,  -700,  -560,  -830, -1350,  -930,  -320,  -700, -1350,  -930,  -560,  -830],
    /* [AP] */ [-2208, -2081, -1623, -2407,  -880,  -320, -2208, -2081,  -880,  -320, -1623, -2407],
    /* [PA] */ [-2694, -1247, -1119, -1729, -1280,  -700, -2694, -1247, -1280,  -700, -1119, -1729],
    /* [IU] */ [-1583,  -822, -1338, -1156, -1460, -1350, -1583,  -822,  1700,   953, -1338, -1156],           
    /* [UI] */ [-1168, -1008, -1199,  -981, -1350,  -930, -1168, -1008,   953,   841, -1199,  -981],
    /* [CI] */ [-1420, -1180, -1270, -1680, -1210,  -560, -1420, -1180, -1210,  -560, -1270, -1680],           
    /* [IC] */ [ -770, -1530, -1450, -1060, -1260,  -830,  -770, -1530, -1260,  -830, -1450, -1060],
];

/// The parameters embedded into whatever dimension is currently 
/// used by the ViennaRNA-style stacking tables.
pub const STACK_EN37: ExtendedStackParams = {
    let mut full: ExtendedStackParams = [[None; E]; E];

    let mut i = 0;
    while i < E {
        let mut j = 0;
        while j < E {
            full[i][j] = Some(STACKPARAMS_EN37[i][j]);
            j += 1;
        }
        i += 1;
    }

    full
};

/// The parameters embedded into whatever dimension is currently 
/// used by the ViennaRNA-style stacking tables.
pub const STACK_ENTH: ExtendedStackParams = {
    let mut full: ExtendedStackParams = [[None; E]; E];

    let mut i = 0;
    while i < E {
        let mut j = 0;
        while j < E {
            full[i][j] = Some(STACKPARAMS_ENTH[i][j]);
            j += 1;
        }
        i += 1;
    }

    full
};
