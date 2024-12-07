//! `mzpeaks` implements the building blocks and machinery for representing peaks
//! in a mass spectrum.
//!
//! It's meant to be used as a building block for other tools and does not provide
//! any I/O machinery for peak lists. For that, consider [`mzdata`](https://crates.io/crates/mzdata)
//!
//! ```rust
//! use mzpeaks::{CentroidPeak, PeakSet, PeakCollection, Tolerance};
//!
//! let peaks = PeakSet::new(vec![
//!     CentroidPeak::new(186.04, 522.0, 0),
//!     CentroidPeak::new(204.07, 9800.0, 1),
//!     CentroidPeak::new(205.07, 150.0, 2)
//! ]);
//!
//! assert_eq!(peaks.search(204.05, Tolerance::Da(0.02)).unwrap(), 1);
//!
//! let peak = match peaks.has_peak(204.05, Tolerance::Da(0.02)) {
//!     Some(p) => p,
//!     None => panic!("Failed to retrieve peak!")
//! };
//!
//! assert!((peak.mz - 204.07).abs() < 1e-6);
//!```

pub mod coordinate;
#[macro_use]
pub mod macros;
pub mod feature;
pub mod feature_map;
pub mod mass_error;
pub mod peak;
pub mod peak_set;
pub mod prelude;
#[cfg(test)]
mod test_data;

pub use crate::coordinate::{
    CoordinateLike, CoordinateLikeMut, CoordinateRange, CoordinateRangeParseError, IndexType,
    IndexedCoordinate, IonMobility, MZLocated, Mass, MassLocated, Time, MZ,
};
pub use crate::mass_error::{Tolerance, ToleranceParsingError};
pub use crate::peak::{
    CentroidLike, CentroidPeak, DeconvolutedCentroidLike, DeconvolutedPeak, IntensityMeasurement,
    IntensityMeasurementMut, KnownCharge, KnownChargeMut,
};
pub use crate::peak_set::{
    DeconvolutedPeakSet, MZPeakSetType, MassPeakSetType, PeakCollection, PeakSet,
};
