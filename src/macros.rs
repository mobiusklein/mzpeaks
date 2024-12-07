/// A set of code generation macros to make a type behave as [`CentroidLike`](crate::CentroidLike)
/// or [`DeconvolutedCentroidLike`](crate::DeconvolutedCentroidLike).

#[macro_export]
macro_rules! implement_mz_coord {
    ($t:ty) => {
        impl<T: $crate::CentroidLike> PartialEq<T> for $t {
            #[inline]
            fn eq(&self, other: &T) -> bool {
                if (self.mz - other.coordinate()).abs() > 1e-3
                    || (self.intensity - other.intensity()).abs() > 1e-3
                {
                    return false;
                }
                true
            }
        }

        impl<T: $crate::CentroidLike> std::cmp::PartialOrd<T> for $t {
            #[inline]
            fn partial_cmp(&self, other: &T) -> Option<std::cmp::Ordering> {
                Some(self.mz.total_cmp(&other.coordinate()).then_with(|| self.intensity.total_cmp(&other.intensity())))
            }
        }

        impl $crate::CoordinateLike<$crate::MZ> for $t {
            #[inline]
            fn coordinate(&self) -> f64 {
                self.mz
            }
        }

        impl $crate::CoordinateLikeMut<$crate::MZ> for $t {
            #[inline]
            fn coordinate_mut(&mut self) -> &mut f64 {
                &mut self.mz
            }
        }

        impl $crate::IntensityMeasurement for $t {
            #[inline]
            fn intensity(&self) -> f32 {
                self.intensity
            }
        }

        impl $crate::IntensityMeasurementMut for $t {
            #[inline]
            fn intensity_mut(&mut self) -> &mut f32 {
                &mut self.intensity
            }
        }
    };
}

#[macro_export]
macro_rules! implement_mass_coord {
    ($t:ty) => {
        impl $crate::CoordinateLike<$crate::Mass> for $t {
            #[inline]
            fn coordinate(&self) -> f64 {
                self.neutral_mass
            }
        }

        impl<T: $crate::DeconvolutedCentroidLike> PartialEq<T> for $t {
            #[inline]
            fn eq(&self, other: &T) -> bool {
                if self.charge != other.charge()
                    || (self.intensity - other.intensity()).abs() > 1e-3
                    || (self.neutral_mass - other.coordinate()).abs() > 1e-3
                {
                    return false;
                }
                true
            }
        }

        impl<T: $crate::DeconvolutedCentroidLike> std::cmp::PartialOrd<T> for $t {
            #[inline]
            fn partial_cmp(&self, other: &T) -> Option<std::cmp::Ordering> {
                match self.neutral_mass.total_cmp(&other.coordinate()) {
                    std::cmp::Ordering::Equal => Some(self.charge.cmp(&other.charge()).then_with(|| self.intensity().total_cmp(&other.intensity()))),
                    x => Some(x)
                }
            }
        }

        impl $crate::IntensityMeasurement for $t {
            #[inline]
            fn intensity(&self) -> f32 {
                self.intensity
            }
        }

        impl $crate::peak::KnownCharge for $t {
            #[inline]
            fn charge(&self) -> i32 {
                self.charge
            }
        }

        impl $crate::CoordinateLikeMut<$crate::Mass> for $t {
            #[inline]
            fn coordinate_mut(&mut self) -> &mut f64 {
                &mut self.neutral_mass
            }
        }

        impl $crate::KnownChargeMut for $t {
            #[inline]
            fn charge_mut(&mut self) -> &mut i32 {
                &mut self.charge
            }
        }

        impl $crate::IntensityMeasurementMut for $t {
            #[inline]
            fn intensity_mut(&mut self) -> &mut f32 {
                &mut self.intensity
            }
        }

    };
}

#[macro_export]
macro_rules! implement_centroid_conversion {
    ($t:ty) => {
        impl From<$t> for $crate::CentroidPeak {
            fn from(peak: $t) -> Self {
                peak.as_centroid()
            }
        }

        impl From<$t> for $crate::peak::MZPoint {
            fn from(peak: $t) -> Self {
                Self {
                    mz: peak.coordinate(),
                    intensity: peak.intensity(),
                    ..Self::default()
                }
            }
        }

        impl From<$crate::CentroidPeak> for $t {
            fn from(peak: $crate::CentroidPeak) -> Self {
                let mut inst = Self {
                    mz: peak.coordinate(),
                    intensity: peak.intensity(),
                    ..Self::default()
                };
                inst.set_index(peak.index);
                inst
            }
        }
    };
}

#[macro_export]
macro_rules! implement_deconvoluted_centroid_conversion {
    ($t:ty) => {
        impl From<$t> for $crate::DeconvolutedPeak {
            fn from(peak: $t) -> Self {
                peak.as_centroid()
            }
        }

        impl From<$crate::DeconvolutedPeak> for $t {
            fn from(peak: $crate::DeconvolutedPeak) -> Self {
                let neutral_mass: f64 = $crate::CoordinateLike::<$crate::Mass>::coordinate(&peak);
                let mut inst = Self {
                    neutral_mass,
                    intensity: peak.intensity(),
                    charge: peak.charge(),
                    ..Self::default()
                };
                inst.set_index(peak.index);
                inst
            }
        }
    };
}

#[macro_export]
/// The first argument is the type, the second is whether
/// to strictly adhere to the indexing requirement, and the
/// third is whether or not to generate `From` specializations.
macro_rules! implement_deconvoluted_centroidlike_inner {
    ($t:ty, true, true) => {
        $crate::implement_mass_coord!($t);
        impl $crate::IndexedCoordinate<$crate::Mass> for $t {
            #[inline]
            fn get_index(&self) -> $crate::IndexType {
                self.index
            }

            #[inline]
            fn set_index(&mut self, index: $crate::IndexType) {
                self.index = index
            }
        }
        $crate::implement_deconvoluted_centroid_conversion!($t);
    };

    ($t:ty, true, false) => {
        $crate::implement_mass_coord!($t);
        impl $crate::IndexedCoordinate<$crate::Mass> for $t {
            #[inline]
            fn get_index(&self) -> $crate::IndexType {
                self.index
            }

            #[inline]
            fn set_index(&mut self, index: $crate::IndexType) {
                self.index = index
            }
        }
    };

    ($t:ty, false, true) => {
        $crate::implement_mass_coord!($t);
        impl $crate::coordinate::IndexedCoordinate<$crate::Mass> for $t {
            fn get_index(&self) -> $crate::IndexType {
                0
            }
            fn set_index(&mut self, _index: $crate::IndexType) {}
        }

        $crate::implement_deconvoluted_centroid_conversion!($t);
    };
}

#[macro_export]
/// The first argument is the type, the second is whether
/// to strictly adhere to the indexing requirement.
macro_rules! implement_deconvoluted_centroidlike {
    ($t:ty, true) => {
        $crate::implement_deconvoluted_centroidlike_inner!($t, true, true);
    };
    ($t:ty, false) => {
        $crate::implement_deconvoluted_centroidlike_inner!($t, false, true);
    };
}

#[macro_export]
/// The first argument is the type, the second is whether
/// to strictly adhere to the indexing requirement, and the
/// third is whether or not to generate `From` specializations.
macro_rules! implement_centroidlike_inner {
    ($t:ty, true, true) => {
        $crate::implement_mz_coord!($t);
        impl $crate::IndexedCoordinate<$crate::MZ> for $t {
            #[inline]
            fn get_index(&self) -> $crate::IndexType {
                self.index
            }

            #[inline]
            fn set_index(&mut self, index: $crate::IndexType) {
                self.index = index
            }
        }

        $crate::implement_centroid_conversion!($t);
    };
    ($t:ty, true, false) => {
        $crate::implement_mz_coord!($t);
        impl $crate::IndexedCoordinate<$crate::MZ> for $t {
            #[inline]
            fn get_index(&self) -> $crate::IndexType {
                self.index
            }

            fn set_index(&mut self, index: $crate::IndexType) {
                self.index = index
            }
        }
    };

    ($t:ty, false, true) => {
        $crate::implement_mz_coord!($t);
        impl $crate::coordinate::IndexedCoordinate<MZ> for $t {
            #[inline]
            fn get_index(&self) -> $crate::IndexType {
                0
            }

            fn set_index(&mut self, _index: $crate::IndexType) {}
        }

        $crate::implement_centroid_conversion!($t);
    };
}

#[macro_export]
/// The first argument is the type, the second is whether
/// to strictly adhere to the indexing requirement.
macro_rules! implement_centroidlike {
    ($t:ty, true) => {
        $crate::implement_centroidlike_inner!($t, true, true);
    };
    ($t:ty, false) => {
        $crate::implement_centroidlike_inner!($t, false, true);
    };
}
