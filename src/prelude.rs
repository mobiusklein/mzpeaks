//! A prelude to bring into scope all the traits of this library.

pub use crate::coordinate::{CoordinateLike, IndexedCoordinate, MZLocated, MassLocated};
pub use crate::peak::{CentroidLike, DeconvolutedCentroidLike, IntensityMeasurement, KnownCharge, KnownChargeMut, IntensityMeasurementMut};
pub use crate::feature::{TimeInterval, FeatureLike, FeatureLikeMut};
pub use crate::peak_set::PeakCollection;
pub use crate::mass_error::Tolerance;