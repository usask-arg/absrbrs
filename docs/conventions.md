# Conventions

These hold for the Rust core and the Python package alike.

## Units

| Quantity | Unit |
|---|---|
| Wavenumber | cm⁻¹ |
| Temperature | K |
| Pressure | Pa |
| Volume mixing ratio | dimensionless |
| Cross section | m² per molecule |

The Python package will offer cm² per molecule on request.

## Line parameters

Line parameters follow HITRAN units and conventions:

- the reference temperature is T_ref = 296 K;
- line intensities S_ref are in cm⁻¹/(molecule cm⁻²) and include the natural isotopic
  abundance, so cross sections are per molecule of the species, with every isotopologue at
  natural abundance;
- widths and shifts are in cm⁻¹ atm⁻¹, and widths are half widths at half maximum.

## Partial pressures

For an absorber with volume mixing ratio VMR at total pressure p, the self-broadening partial
pressure is p_self = VMR · p, and the air-broadening partial pressure is p_air = (1 − VMR) · p.

## Wavenumber grids

Grids passed to the Rust core must be strictly increasing; they need not be uniform. The Python
package accepts grids in any order and returns results in that order.
