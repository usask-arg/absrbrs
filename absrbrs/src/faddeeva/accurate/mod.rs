//! The Accurate tier: a Rust port of Steven G. Johnson's Faddeeva Package (`Faddeeva.cc`,
//! revision of 12 May 2015), MIT licence; see `absrbrs/LICENSES/Faddeeva-MIT.txt`.

mod erfcx;
mod w_im;

pub(crate) use erfcx::erfcx;
pub(crate) use w_im::w_im;
