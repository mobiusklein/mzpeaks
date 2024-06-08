//! A prelude to bring into scope all the traits of this library.

pub use crate::coordinate::{CoordinateLike, IndexedCoordinate, MZLocated, MassLocated};
pub use crate::feature::{FeatureLike, FeatureLikeMut, TimeInterval, SplittableFeatureLike, TimeArray};
pub use crate::feature_map::{FeatureMapLike, FeatureMapLikeMut};
pub use crate::mass_error::Tolerance;
pub use crate::peak::{
    CentroidLike, DeconvolutedCentroidLike, IntensityMeasurement, IntensityMeasurementMut,
    KnownCharge, KnownChargeMut,
};
pub use crate::peak_set::{PeakCollection, PeakCollectionMut};
