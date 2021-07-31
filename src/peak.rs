//! A peak is the most atomic unit of a (processed) mass spectrum. It represents
//! a location in one (or more) coordinate spaces with an measured intensity value
//!
//! All peak-like types implement [`PartialEq`], [`PartialOrd`], and [`CoordinateLike`]
//!

use std::cmp;
use std::fmt;
use std::hash;

use crate::coordinate::{CoordinateLike, IndexType, IndexedCoordinate, Mass, MZ};

/// An intensity measurement is an entity that has a measured intensity
/// of whateger it is.
pub trait IntensityMeasurement {
    fn intensity(&self) -> f32;
}

/// A [`CentroidLike`] entity is indexed in m/z coordinate space and
/// is an [`IntensityMeasurement`]
pub trait CentroidLike: IndexedCoordinate<MZ> + IntensityMeasurement {
    fn as_centroid(&self) -> CentroidPeak {
        CentroidPeak {
            mz: self.coordinate(),
            intensity: self.intensity(),
            index: self.get_index(),
        }
    }
}

/// A known charge has a determined charge state value
pub trait KnownCharge {
    fn charge(&self) -> i32;
}

/// A [`DeconvolutedCentroid`] entity is indexed in the neutral mass
/// coordiante space, has known charge state and an aggregated intensity
/// measurement. Any [`DeconvolutedCentroid`] can be converted into
/// a [`DeconvolutedPeak`]
pub trait DeconvolutedCentroid:
    IndexedCoordinate<Mass> + IntensityMeasurement + KnownCharge
{
    fn as_centroid(&self) -> DeconvolutedPeak {
        DeconvolutedPeak {
            neutral_mass: self.coordinate(),
            intensity: self.intensity(),
            charge: self.charge(),
            index: self.get_index(),
        }
    }
}

/// Represent a single m/z coordinate with an
/// intensity and an index. Nearly the most basic
/// peak representation for peak-picked data.
#[derive(Default, Clone, Debug)]
pub struct CentroidPeak {
    pub mz: f64,
    pub intensity: f32,
    pub index: IndexType,
}

impl CentroidPeak {
    pub fn new(mz: f64, intensity: f32, index: IndexType) -> CentroidPeak {
        CentroidPeak {
            mz,
            intensity,
            index,
        }
    }
}

impl fmt::Display for CentroidPeak {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "CentroidPeak({}, {}, {})",
            self.mz, self.intensity, self.index
        )
    }
}

impl hash::Hash for CentroidPeak {
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        let mz_val: i64 = self.coordinate().round() as i64;
        mz_val.hash(state);
    }
}

impl<T: CentroidLike> cmp::PartialOrd<T> for CentroidPeak {
    #[inline]
    fn partial_cmp(&self, other: &T) -> Option<cmp::Ordering> {
        self.mz.partial_cmp(&other.coordinate())
    }
}

impl<T: CentroidLike> cmp::PartialEq<T> for CentroidPeak {
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

impl CoordinateLike<MZ> for CentroidPeak {
    #[inline]
    fn coordinate(&self) -> f64 {
        self.mz
    }
}

impl IndexedCoordinate<MZ> for CentroidPeak {
    #[inline]
    fn get_index(&self) -> IndexType {
        self.index
    }

    #[inline]
    fn set_index(&mut self, index: IndexType) {
        self.index = index;
    }
}

impl IntensityMeasurement for CentroidPeak {
    #[inline]
    fn intensity(&self) -> f32 {
        self.intensity
    }
}

impl<T: IndexedCoordinate<MZ> + IntensityMeasurement> CentroidLike for T {}

#[derive(Default, Clone, Debug)]
/// Represent a single neutral mass coordinate with an
/// intensity, a known charge and an index.
pub struct DeconvolutedPeak {
    pub neutral_mass: f64,
    pub intensity: f32,
    pub charge: i32,
    pub index: IndexType,
}

impl fmt::Display for DeconvolutedPeak {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "DeconvolutedPeak({}, {}, {}, {})",
            self.neutral_mass, self.intensity, self.charge, self.index
        )
    }
}

impl hash::Hash for DeconvolutedPeak {
    fn hash<H: hash::Hasher>(&self, state: &mut H) {
        let neutral_mass: i64 = self.neutral_mass.round() as i64;
        neutral_mass.hash(state);
    }
}

impl<T: DeconvolutedCentroid> cmp::PartialOrd<T> for DeconvolutedPeak {
    #[inline]
    fn partial_cmp(&self, other: &T) -> Option<cmp::Ordering> {
        self.neutral_mass.partial_cmp(&other.coordinate())
    }
}

impl<T: DeconvolutedCentroid> cmp::PartialEq<T> for DeconvolutedPeak {
    #[inline]
    fn eq(&self, other: &T) -> bool {
        if (self.neutral_mass - other.coordinate()).abs() > 1e-3
            || self.charge != other.charge()
            || (self.intensity - other.intensity()).abs() > 1e-3
        {
            return false;
        }
        true
    }
}

impl CoordinateLike<Mass> for DeconvolutedPeak {
    #[inline]
    fn coordinate(&self) -> f64 {
        self.neutral_mass
    }
}

impl IndexedCoordinate<Mass> for DeconvolutedPeak {
    #[inline]
    fn get_index(&self) -> IndexType {
        self.index
    }

    #[inline]
    fn set_index(&mut self, index: IndexType) {
        self.index = index;
    }
}

impl CoordinateLike<MZ> for DeconvolutedPeak {
    #[inline]
    fn coordinate(&self) -> f64 {
        let charge_carrier: f64 = 1.007276;
        let charge = self.charge as f64;
        (self.neutral_mass - charge_carrier * charge) / charge
    }
}

impl IntensityMeasurement for DeconvolutedPeak {
    #[inline]
    fn intensity(&self) -> f32 {
        self.intensity
    }
}

impl KnownCharge for DeconvolutedPeak {
    #[inline]
    fn charge(&self) -> i32 {
        self.charge
    }
}

impl<T: IndexedCoordinate<Mass> + IntensityMeasurement + KnownCharge> DeconvolutedCentroid for T {}

#[derive(Debug, Clone, Default)]
pub struct MZPoint {
    pub mz: f64,
    pub intensity: f32,
}

impl MZPoint {
    pub fn new(mz: f64, intensity: f32) -> MZPoint {
        MZPoint { mz, intensity }
    }
}

impl<T: CentroidLike> cmp::PartialOrd<T> for MZPoint {
    #[inline]
    fn partial_cmp(&self, other: &T) -> Option<cmp::Ordering> {
        self.mz.partial_cmp(&other.coordinate())
    }
}

impl<T: CentroidLike> cmp::PartialEq<T> for MZPoint {
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

impl CoordinateLike<MZ> for MZPoint {
    fn coordinate(&self) -> f64 {
        self.mz
    }
}

impl IndexedCoordinate<MZ> for MZPoint {
    fn get_index(&self) -> IndexType {
        0
    }

    fn set_index(&mut self, _index: IndexType) {}
}

impl IntensityMeasurement for MZPoint {
    fn intensity(&self) -> f32 {
        self.intensity
    }
}

impl From<MZPoint> for CentroidPeak {
    fn from(peak: MZPoint) -> Self {
        peak.as_centroid()
    }
}

impl From<CentroidPeak> for MZPoint {
    fn from(peak: CentroidPeak) -> Self {
        Self {
            mz: peak.coordinate(),
            intensity: peak.intensity(),
        }
    }
}
