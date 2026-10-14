//! Line lists as sorted structure-of-arrays columns (spec §5.1).

use crate::{Error, Result};

/// A per-line floating-point parameter of a spectral line.
///
/// The discriminants index the column storage of a [`LineList`], so the declaration order is
/// part of the contract: use `field as usize` to index per-field arrays.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum LineField {
    /// Line centre ν₀, cm⁻¹.
    Nu,
    /// Intensity at 296 K, cm⁻¹/(molecule cm⁻²).
    SRef,
    /// Lower-state energy, cm⁻¹.
    ELower,
    /// Air-broadened half-width at 296 K and 1 atm, cm⁻¹ atm⁻¹.
    GammaAir,
    /// Self-broadened half-width at 296 K and 1 atm, cm⁻¹ atm⁻¹.
    GammaSelf,
    /// Temperature exponent of the air-broadened half-width (dimensionless).
    NAir,
    /// Temperature exponent of the self-broadened half-width (optional; default n_air).
    NSelf,
    /// Air pressure shift at 296 K and 1 atm, cm⁻¹ atm⁻¹.
    DeltaAir,
    /// Temperature coefficient of the air shift, cm⁻¹ atm⁻¹ K⁻¹ (optional; default 0).
    DeltapAir,
    /// Self pressure shift at 296 K and 1 atm, cm⁻¹ atm⁻¹ (optional; default δ_air).
    DeltaSelf,
}

impl LineField {
    /// Every field, in declaration (storage) order.
    pub const ALL: [LineField; 10] = [
        LineField::Nu,
        LineField::SRef,
        LineField::ELower,
        LineField::GammaAir,
        LineField::GammaSelf,
        LineField::NAir,
        LineField::NSelf,
        LineField::DeltaAir,
        LineField::DeltapAir,
        LineField::DeltaSelf,
    ];

    /// The field's name, as used in the Python API and in error messages.
    pub fn name(self) -> &'static str {
        match self {
            LineField::Nu => "nu",
            LineField::SRef => "s_ref",
            LineField::ELower => "e_lower",
            LineField::GammaAir => "gamma_air",
            LineField::GammaSelf => "gamma_self",
            LineField::NAir => "n_air",
            LineField::NSelf => "n_self",
            LineField::DeltaAir => "delta_air",
            LineField::DeltapAir => "deltap_air",
            LineField::DeltaSelf => "delta_self",
        }
    }

    /// The field with the given [`name`](Self::name), if there is one.
    pub fn from_name(name: &str) -> Option<LineField> {
        LineField::ALL.into_iter().find(|f| f.name() == name)
    }

    /// Whether the field may be absent for a line (n_self, δ′_air, δ_self).
    pub fn is_optional(self) -> bool {
        matches!(
            self,
            LineField::NSelf | LineField::DeltapAir | LineField::DeltaSelf
        )
    }

    /// Checks one value of this field; the error says what is wrong with it.
    fn check_value(self, value: f64) -> std::result::Result<(), &'static str> {
        if self.is_optional() {
            // NaN means "absent for this line".
            return if value.is_infinite() {
                Err("is infinite")
            } else {
                Ok(())
            };
        }
        if !value.is_finite() {
            return Err("is not finite");
        }
        match self {
            LineField::Nu if value <= 0.0 => Err("is not positive"),
            LineField::SRef | LineField::GammaAir | LineField::GammaSelf if value < 0.0 => {
                Err("is negative")
            }
            _ => Ok(()),
        }
    }
}

/// Unsorted input columns for [`LineList::from_columns`].
///
/// All columns must have the same length. In the optional columns a NaN entry means "absent for
/// this line", and a `None` column means "absent for every line".
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LineColumns {
    /// HITRAN molecule number of each line.
    pub molecule: Vec<u16>,
    /// HITRAN local isotopologue number of each line (1 is the most abundant).
    pub isotopologue: Vec<u16>,
    /// Line centres ν₀, cm⁻¹.
    pub nu: Vec<f64>,
    /// Intensities at 296 K, cm⁻¹/(molecule cm⁻²).
    pub s_ref: Vec<f64>,
    /// Lower-state energies, cm⁻¹.
    pub e_lower: Vec<f64>,
    /// Air-broadened half-widths, cm⁻¹ atm⁻¹.
    pub gamma_air: Vec<f64>,
    /// Self-broadened half-widths, cm⁻¹ atm⁻¹.
    pub gamma_self: Vec<f64>,
    /// Temperature exponents of the air-broadened half-widths.
    pub n_air: Vec<f64>,
    /// Air pressure shifts, cm⁻¹ atm⁻¹.
    pub delta_air: Vec<f64>,
    /// Temperature exponents of the self-broadened half-widths (default n_air).
    pub n_self: Option<Vec<f64>>,
    /// Self pressure shifts, cm⁻¹ atm⁻¹ (default δ_air).
    pub delta_self: Option<Vec<f64>>,
    /// Temperature coefficients of the air shifts, cm⁻¹ atm⁻¹ K⁻¹ (default 0).
    pub deltap_air: Option<Vec<f64>>,
}

/// An immutable list of spectral lines, stored as columns and sorted by wavenumber.
///
/// Every line index (in [`value`](Self::value), [`has_value`](Self::has_value) and the
/// positions within [`column`](Self::column)) refers to this sorted order. The seven required
/// fields are always present; the optional ones resolve to their defaults, see
/// [`value`](Self::value).
#[derive(Clone, Debug, PartialEq)]
pub struct LineList {
    molecule: Vec<u16>,
    isotopologue: Vec<u16>,
    /// Indexed by `LineField as usize`; the required fields are always `Some`.
    columns: [Option<Vec<f64>>; 10],
}

/// Checks that a column has `expected` entries and that every entry is valid for `field`.
fn validate_column(field: LineField, values: &[f64], expected: usize) -> Result<()> {
    let name = field.name();
    if values.len() != expected {
        return Err(Error::InvalidInput(format!(
            "line column `{name}` has {} entries, expected {expected}",
            values.len()
        )));
    }
    for (i, &v) in values.iter().enumerate() {
        if let Err(why) = field.check_value(v) {
            return Err(Error::InvalidInput(format!(
                "line column `{name}`: value {v} at index {i} {why}"
            )));
        }
    }
    Ok(())
}

/// Gathers `values[i]` for each `i` in `indices`.
fn gather<T: Copy>(values: &[T], indices: &[usize]) -> Vec<T> {
    indices.iter().map(|&i| values[i]).collect()
}

impl LineList {
    /// Builds a list from unsorted columns, sorting the lines by wavenumber (stable).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] if the columns differ in length, if a required value is not
    /// finite, if ν₀ ≤ 0, if S_ref or a half-width is negative, or if an optional value is
    /// infinite. The message names the column and the first offending index.
    pub fn from_columns(columns: LineColumns) -> Result<LineList> {
        let n = columns.nu.len();
        for (name, len) in [
            ("molecule", columns.molecule.len()),
            ("isotopologue", columns.isotopologue.len()),
        ] {
            if len != n {
                return Err(Error::InvalidInput(format!(
                    "line column `{name}` has {len} entries, expected {n}"
                )));
            }
        }
        let LineColumns {
            molecule,
            isotopologue,
            nu,
            s_ref,
            e_lower,
            gamma_air,
            gamma_self,
            n_air,
            delta_air,
            n_self,
            delta_self,
            deltap_air,
        } = columns;
        let mut storage: [Option<Vec<f64>>; 10] = Default::default();
        let fields = [
            (LineField::Nu, Some(nu)),
            (LineField::SRef, Some(s_ref)),
            (LineField::ELower, Some(e_lower)),
            (LineField::GammaAir, Some(gamma_air)),
            (LineField::GammaSelf, Some(gamma_self)),
            (LineField::NAir, Some(n_air)),
            (LineField::NSelf, n_self),
            (LineField::DeltaAir, Some(delta_air)),
            (LineField::DeltapAir, deltap_air),
            (LineField::DeltaSelf, delta_self),
        ];
        for (field, values) in fields {
            if let Some(values) = &values {
                validate_column(field, values, n)?;
            }
            storage[field as usize] = values;
        }
        let mut list = LineList {
            molecule,
            isotopologue,
            columns: storage,
        };
        list.sort_by_nu();
        Ok(list)
    }

    /// Same as [`from_columns`](Self::from_columns); the entry point for the file readers.
    // The readers arrive in later tasks.
    #[allow(dead_code)]
    pub(crate) fn from_unsorted(columns: LineColumns) -> Result<LineList> {
        LineList::from_columns(columns)
    }

    /// Stably sorts every column by wavenumber.
    fn sort_by_nu(&mut self) {
        let nu = self.nu();
        if nu.is_sorted_by(|a, b| a.total_cmp(b).is_le()) {
            return;
        }
        let mut order: Vec<usize> = (0..nu.len()).collect();
        order.sort_by(|&a, &b| nu[a].total_cmp(&nu[b]));
        self.reorder(&order);
    }

    /// Replaces every column by the rows listed in `indices`, in that order.
    fn reorder(&mut self, indices: &[usize]) {
        self.molecule = gather(&self.molecule, indices);
        self.isotopologue = gather(&self.isotopologue, indices);
        for column in self.columns.iter_mut().flatten() {
            *column = gather(column, indices);
        }
    }

    /// The wavenumbers, which are always present.
    fn nu(&self) -> &[f64] {
        self.columns[LineField::Nu as usize]
            .as_deref()
            .unwrap_or_default()
    }

    /// The number of lines.
    pub fn len(&self) -> usize {
        self.molecule.len()
    }

    /// Whether the list has no lines.
    pub fn is_empty(&self) -> bool {
        self.molecule.is_empty()
    }

    /// The HITRAN molecule number of each line.
    pub fn molecule(&self) -> &[u16] {
        &self.molecule
    }

    /// The HITRAN local isotopologue number of each line.
    pub fn isotopologue(&self) -> &[u16] {
        &self.isotopologue
    }

    /// The stored column for a field, or `None` for an optional field that is absent for every
    /// line. Individual NaN entries of an optional column mean "absent for this line".
    pub fn column(&self, field: LineField) -> Option<&[f64]> {
        self.columns[field as usize].as_deref()
    }

    /// The value of a field for line `i`, with the defaults of absent optional fields applied:
    /// n_self is n_air, δ_self is δ_air and δ′_air is 0.
    ///
    /// # Panics
    ///
    /// If `i >= self.len()`.
    pub fn value(&self, field: LineField, i: usize) -> f64 {
        if let Some(v) = self.stored(field, i) {
            return v;
        }
        match field {
            LineField::NSelf => self.value(LineField::NAir, i),
            LineField::DeltaSelf => self.value(LineField::DeltaAir, i),
            // δ′_air; the required fields are always stored.
            _ => 0.0,
        }
    }

    /// Whether line `i` has its own value for a field, rather than the default of
    /// [`value`](Self::value). Always true for the required fields.
    ///
    /// # Panics
    ///
    /// If `i >= self.len()`.
    pub fn has_value(&self, field: LineField, i: usize) -> bool {
        self.stored(field, i).is_some()
    }

    /// The stored, non-NaN value of a field for line `i`.
    fn stored(&self, field: LineField, i: usize) -> Option<f64> {
        assert!(
            i < self.len(),
            "line index {i} out of range for a list of {} lines",
            self.len()
        );
        self.column(field).map(|c| c[i]).filter(|v| !v.is_nan())
    }

    /// The lines that match every given filter: the molecule number, one of the isotopologue
    /// numbers, and a wavenumber in the inclusive range. `None` places no constraint. The
    /// result is sorted, like the list itself.
    pub fn select(
        &self,
        molecule: Option<u16>,
        isotopologues: Option<&[u16]>,
        wavenumber: Option<(f64, f64)>,
    ) -> LineList {
        let nu = self.nu();
        let keep: Vec<usize> = (0..self.len())
            .filter(|&i| {
                molecule.is_none_or(|m| self.molecule[i] == m)
                    && isotopologues.is_none_or(|iso| iso.contains(&self.isotopologue[i]))
                    && wavenumber.is_none_or(|(lo, hi)| (lo..=hi).contains(&nu[i]))
            })
            .collect();
        let mut selected = self.clone();
        selected.reorder(&keep);
        selected
    }

    /// A copy of the list with one field replaced. Replacing an optional field makes its
    /// column present (NaN entries still mean "absent for this line"); replacing the
    /// wavenumbers re-sorts the lines.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] if `values` has the wrong length or an invalid entry, by the same
    /// rules as [`from_columns`](Self::from_columns).
    pub fn with_column(&self, field: LineField, values: Vec<f64>) -> Result<LineList> {
        validate_column(field, &values, self.len())?;
        let mut list = self.clone();
        list.columns[field as usize] = Some(values);
        if field == LineField::Nu {
            list.sort_by_nu();
        }
        Ok(list)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn columns(nu: Vec<f64>) -> LineColumns {
        let n = nu.len();
        LineColumns {
            molecule: vec![2; n],
            isotopologue: (1..=n as u16).collect(),
            s_ref: (0..n).map(|i| 1e-21 * (i + 1) as f64).collect(),
            e_lower: vec![100.0; n],
            gamma_air: vec![0.07; n],
            gamma_self: vec![0.09; n],
            n_air: vec![0.75; n],
            delta_air: vec![-0.004; n],
            nu,
            ..Default::default()
        }
    }

    #[test]
    fn sorts_by_wavenumber_keeping_rows_together() {
        let lines = LineList::from_columns(columns(vec![3.0, 1.0, 2.0])).unwrap();
        assert_eq!(lines.column(LineField::Nu).unwrap(), &[1.0, 2.0, 3.0]);
        assert_eq!(lines.isotopologue(), &[2, 3, 1]);
        assert_eq!(
            lines.column(LineField::SRef).unwrap(),
            // As `columns()` computes them: 3.0 * 1e-21 is not the literal 3e-21.
            &[2.0 * 1e-21, 3.0 * 1e-21, 1e-21]
        );
        assert_eq!(lines.len(), 3);
    }

    #[test]
    fn optional_fields_default_per_line() {
        let mut c = columns(vec![1.0, 2.0]);
        c.n_self = Some(vec![0.6, f64::NAN]);
        let lines = LineList::from_columns(c).unwrap();
        assert_eq!(lines.value(LineField::NSelf, 0), 0.6);
        assert_eq!(lines.value(LineField::NSelf, 1), 0.75); // NaN: absent, so n_air
        assert!(lines.has_value(LineField::NSelf, 0) && !lines.has_value(LineField::NSelf, 1));
        assert_eq!(lines.value(LineField::DeltaSelf, 0), -0.004); // column absent: delta_air
        assert_eq!(lines.value(LineField::DeltapAir, 1), 0.0);
        assert!(lines.column(LineField::DeltaSelf).is_none());
    }

    #[test]
    fn rejects_bad_columns() {
        let mut c = columns(vec![1.0, 2.0]);
        c.gamma_air.pop();
        assert!(matches!(
            LineList::from_columns(c),
            Err(Error::InvalidInput(_))
        ));
        let mut c = columns(vec![1.0, 2.0]);
        c.s_ref[1] = f64::NAN;
        assert!(matches!(
            LineList::from_columns(c),
            Err(Error::InvalidInput(_))
        ));
        let mut c = columns(vec![1.0, 0.0]);
        c.nu[1] = 0.0;
        assert!(matches!(
            LineList::from_columns(c),
            Err(Error::InvalidInput(_))
        ));
        let mut c = columns(vec![1.0, 2.0]);
        c.n_self = Some(vec![0.5, f64::INFINITY]);
        assert!(matches!(
            LineList::from_columns(c),
            Err(Error::InvalidInput(_))
        ));
    }

    #[test]
    fn select_filters() {
        let mut c = columns(vec![1.0, 2.0, 3.0, 4.0]);
        c.molecule = vec![1, 2, 2, 2];
        let lines = LineList::from_columns(c).unwrap();
        let s = lines.select(Some(2), Some(&[3, 4]), Some((2.5, 4.0)));
        assert_eq!(s.column(LineField::Nu).unwrap(), &[3.0, 4.0]);
        assert_eq!(lines.select(None, None, None), lines);
        assert!(lines.select(Some(9), None, None).is_empty());
    }

    #[test]
    fn with_column_replaces_and_resorts() {
        let lines = LineList::from_columns(columns(vec![1.0, 2.0])).unwrap();
        let g = lines
            .with_column(LineField::GammaAir, vec![0.01, 0.02])
            .unwrap();
        assert_eq!(g.column(LineField::GammaAir).unwrap(), &[0.01, 0.02]);
        let r = lines.with_column(LineField::Nu, vec![5.0, 4.0]).unwrap();
        assert_eq!(r.column(LineField::Nu).unwrap(), &[4.0, 5.0]);
        assert_eq!(r.isotopologue(), &[2, 1]);
        let s = lines
            .with_column(LineField::NSelf, vec![0.5, f64::NAN])
            .unwrap();
        assert_eq!(s.value(LineField::NSelf, 1), 0.75);
        assert!(lines.with_column(LineField::NAir, vec![1.0]).is_err());
    }

    #[test]
    fn field_names_round_trip() {
        for f in LineField::ALL {
            assert_eq!(LineField::from_name(f.name()), Some(f));
        }
        assert_eq!(LineField::from_name("sw"), None);
        assert_eq!(LineField::ALL.iter().filter(|f| f.is_optional()).count(), 3);
    }
}
