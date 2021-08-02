//! `mzpeaks` implements the building blocks and machinery for representing peaks
//! in a mass spectrum.
//!
//! It's meant to be used as a building block for other tools and does not provide
//! any I/O machinery for peak lists
//!
//! ```rust
//! use mzpeaks::{CentroidPeak, PeakSet, PeakCollection, MassErrorType};
//!
//! let peaks = PeakSet::new(vec![
//!     CentroidPeak::new(186.04, 522.0, 0),
//!     CentroidPeak::new(204.07, 9800.0, 1),
//!     CentroidPeak::new(205.07, 150.0, 2)
//! ]);
//!
//! assert_eq!(peaks.search(204.05, 0.02, MassErrorType::Absolute).unwrap(), 1);
//!
//! let peak = match peaks.has_peak(204.05, 0.02, MassErrorType::Absolute) {
//!     Some(p) => p,
//!     None => panic!("Failed to retrieve peak!")
//! };
//!
//! assert!((peak.mz - 204.07).abs() < 1e-6);
//!```

pub mod coordinate;
pub mod macros;
pub mod mass_error;
pub mod peak;
pub mod peak_set;
pub mod prelude;
#[cfg(test)]
mod test_data;

pub use crate::coordinate::{
    CoordinateDimension, CoordinateLike, IndexType, IndexedCoordinate, MZLocated, Mass,
    MassLocated, MZ,
};
pub use crate::mass_error::MassErrorType;
pub use crate::peak::{
    CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak, IntensityMeasurement,
};
pub use crate::peak_set::{
    DeconvolutedPeakSet, MZPeakSetType, MassPeakSetType, PeakCollection, PeakSet,
};
