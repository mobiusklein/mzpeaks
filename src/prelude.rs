//! A prelude to bring into scope all the traits of this library.

pub use crate::coordinate::{
    CoordinateLike, CoordinateLikeMut, CoordinateSystem as _, HasProximity, IndexedCoordinate,
    IonMobilityLocated, MZLocated, MassLocated, Span1D, Span2D, TimeLocated,
};
pub use crate::feature::{
    BuildFromPeak, FeatureLike, FeatureLikeMut, NDFeatureLike, NDFeatureLikeMut,
    AsPeakIter, SplittableFeatureLike, TimeArray, TimeInterval, PeakSeries
};
pub use crate::feature_map::{FeatureMapLike, FeatureMapLikeMut};
pub use crate::mass_error::Tolerance;
pub use crate::peak::{
    CentroidLike, DeconvolutedCentroidLike, IntensityMeasurement, IntensityMeasurementMut,
    KnownCharge, KnownChargeMut,
};
pub use crate::peak_set::{PeakCollection, PeakCollectionMut};
