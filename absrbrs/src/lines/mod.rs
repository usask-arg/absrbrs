//! Spectroscopic line lists and isotopologue metadata (spec §5.1, §5.6).

mod isotopologue_table;
mod isotopologues;
mod list;

pub use isotopologues::{Isotopologue, isotopologue, isotopologues, molecule_id, molecule_name};
pub use list::{LineColumns, LineField, LineList};
