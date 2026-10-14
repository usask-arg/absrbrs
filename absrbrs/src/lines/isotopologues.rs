//! Isotopologue metadata in HITRAN's numbering: masses and natural abundances (from HAPI's
//! table, MIT licence; see absrbrs/LICENSES/HAPI-MIT.txt).

use super::isotopologue_table::ISOTOPOLOGUES;

/// One isotopologue: HITRAN molecule and (local) isotopologue numbers, and its properties.
///
/// Only the crate's table constructs these; more fields may be added.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct Isotopologue {
    /// HITRAN molecule number (1 = H₂O, 2 = CO₂, …).
    pub molecule: u16,
    /// HITRAN local isotopologue number (1, 2, …; `.par` codes 0, A, B are 10, 11, 12).
    pub isotopologue: u16,
    /// HITRAN global isotopologue id.
    pub global_id: u16,
    /// Formula, e.g. `H2(16O)`.
    pub formula: &'static str,
    /// Molecule name, e.g. `H2O`.
    pub molecule_name: &'static str,
    /// Natural terrestrial abundance (included in HITRAN line intensities).
    pub abundance: f64,
    /// Molar mass, g/mol.
    pub mass: f64,
}

/// The isotopologue (molecule, isotopologue), if the table has it.
pub fn isotopologue(molecule: u16, isotopologue: u16) -> Option<&'static Isotopologue> {
    ISOTOPOLOGUES
        .binary_search_by_key(&(molecule, isotopologue), |i| (i.molecule, i.isotopologue))
        .ok()
        .map(|k| &ISOTOPOLOGUES[k])
}

/// The HITRAN molecule number of `name` (exact, case-sensitive: "CO2", "H2O", …).
pub fn molecule_id(name: &str) -> Option<u16> {
    ISOTOPOLOGUES
        .iter()
        .find(|i| i.molecule_name == name)
        .map(|i| i.molecule)
}

/// The name of HITRAN molecule `molecule`.
pub fn molecule_name(molecule: u16) -> Option<&'static str> {
    ISOTOPOLOGUES
        .iter()
        .find(|i| i.molecule == molecule)
        .map(|i| i.molecule_name)
}

/// Every isotopologue in the table, sorted by (molecule, isotopologue).
pub fn isotopologues() -> &'static [Isotopologue] {
    &ISOTOPOLOGUES
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_sorted_and_unique() {
        assert!(
            isotopologues()
                .windows(2)
                .all(|w| (w[0].molecule, w[0].isotopologue) < (w[1].molecule, w[1].isotopologue))
        );
        assert_eq!(isotopologues().len(), 156);
    }

    #[test]
    fn known_entries() {
        let h2o = isotopologue(1, 1).unwrap();
        assert_eq!(
            (h2o.formula, h2o.molecule_name, h2o.abundance, h2o.mass),
            ("H2(16O)", "H2O", 0.9973173, 18.010565)
        );
        let co2 = isotopologue(2, 12).unwrap();
        assert_eq!((co2.formula, co2.mass), ("(13C)(17O)2", 47.001618));
        assert!(isotopologue(2, 99).is_none());
    }

    #[test]
    fn molecule_names() {
        assert_eq!(molecule_id("CO2"), Some(2));
        assert_eq!(molecule_id("O2"), Some(7));
        assert_eq!(molecule_id("co2"), None);
        assert_eq!(molecule_name(1), Some("H2O"));
        assert_eq!(molecule_name(0), None);
    }
}
