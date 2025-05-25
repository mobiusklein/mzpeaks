use std::fmt::Display;

use num_traits::{Float, FromPrimitive};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};


/// An enum over the different coordinate planes
#[derive(Debug, PartialEq, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Dimension {
    MZ(MZ),
    Mass(Mass),
    Time(Time),
    IonMobility(IonMobility),
    Dimensionless(Dimensionless),
}

macro_rules! dim_dispatch {
    ($d:ident, $f:tt) => {
        match $d {
            Dimension::MZ(_) => <MZ as CoordinateSystem>::$f(),
            Dimension::Mass(_) => <Mass as CoordinateSystem>::$f(),
            Dimension::Time(_) => <Time as CoordinateSystem>::$f(),
            Dimension::IonMobility(_) => <IonMobility as CoordinateSystem>::$f(),
            Dimension::Dimensionless(_) => <Dimensionless as CoordinateSystem>::$f(),
        }
    };
}

impl Dimension {
    pub const fn name(&self) -> &'static str {
        match self {
            Dimension::MZ(_) => "m/z",
            Dimension::Mass(_) => "neutral mass",
            Dimension::Time(_) => "time",
            Dimension::IonMobility(_) => "ion mobility",
            Dimension::Dimensionless(_) => "",
        }
    }

    pub fn minimum_value(&self) -> f64 {
        dim_dispatch!(self, minimum_value)
    }

    pub fn maximum_value(&self) -> f64 {
        dim_dispatch!(self, maximum_value)
    }
}

impl Display for Dimension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
/// The Mass To Charge Ratio (m/z) coordinate system
pub struct MZ();

impl MZ {
    /// Access the m/z of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<MZ>>(inst: &T) -> f64 {
        CoordinateLike::<MZ>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
/// The Mass coordinate system
pub struct Mass();

impl Mass {
    /// Access the neutral mass of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<Mass>>(inst: &T) -> f64 {
        CoordinateLike::<Mass>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
/// The Event Time coordinate system
pub struct Time();
impl Time {
    /// Access the elapsed time of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<Time>>(inst: &T) -> f64 {
        CoordinateLike::<Time>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
/// The Ion Mobility Time coordinate system
pub struct IonMobility();
impl IonMobility {
    /// Access the ion mobility time unit of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<IonMobility>>(inst: &T) -> f64 {
        CoordinateLike::<IonMobility>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Dimensionless();

#[allow(unused)]
impl Dimensionless {
    /// Access some unitless position coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<Dimensionless>>(inst: &T) -> f64 {
        CoordinateLike::<Dimensionless>::coordinate(inst)
    }
}

/// Describe a coordinate system as an object itself rather than as a type parameter
pub trait CoordinateSystem: Sized {
    #[inline]
    fn coordinate<T: CoordinateLike<Self>>(inst: &T) -> f64 {
        CoordinateLike::<Self>::coordinate(inst)
    }

    fn coordinate_mut<T: CoordinateLikeMut<Self>>(inst: &mut T) -> &mut f64 {
        CoordinateLikeMut::<Self>::coordinate_mut(inst)
    }

    fn name_of(&self) -> &'static str {
        Self::name()
    }

    fn dimension() -> Dimension;

    fn name() -> &'static str {
        Self::dimension().name()
    }

    fn minimum_value() -> f64 {
        0.0
    }

    fn maximum_value() -> f64 {
        f64::INFINITY
    }
}

impl CoordinateSystem for MZ {
    fn dimension() -> Dimension {
        Dimension::MZ(Self())
    }
}
impl CoordinateSystem for Mass {
    fn dimension() -> Dimension {
        Dimension::Mass(Self())
    }
}
impl CoordinateSystem for Time {
    fn dimension() -> Dimension {
        Dimension::Time(Self())
    }
}
impl CoordinateSystem for IonMobility {
    fn dimension() -> Dimension {
        Dimension::IonMobility(Self())
    }
}
impl CoordinateSystem for Dimensionless {
    fn dimension() -> Dimension {
        Dimension::Dimensionless(Self())
    }
}

/// Denote a type has a coordinate value on coordinate system `T`
pub trait CoordinateLike<T>: PartialOrd {
    /// The trait method for accessing the coordinate of the object on coordinate
    /// system `T`
    fn coordinate(&self) -> f64;
}

/// A [`CoordinateLike`] structure whose coordinate is mutable
pub trait CoordinateLikeMut<T>: CoordinateLike<T> {
    fn coordinate_mut(&mut self) -> &mut f64;
}

/// A named coordinate system membership for neutral mass
pub trait MassLocated: CoordinateLike<Mass> {
    #[inline]
    fn neutral_mass(&self) -> f64 {
        CoordinateLike::<Mass>::coordinate(self)
    }
}

/// A named coordinate system membership for m/z
pub trait MZLocated: CoordinateLike<MZ> {
    #[inline]
    fn mz(&self) -> f64 {
        CoordinateLike::<MZ>::coordinate(self)
    }
}

pub trait TimeLocated: CoordinateLike<Time> {
    #[inline]
    fn time(&self) -> f64 {
        CoordinateLike::<Time>::coordinate(self)
    }
}

pub trait IonMobilityLocated: CoordinateLike<IonMobility> {
    #[inline]
    fn ion_mobility(&self) -> f64 {
        CoordinateLike::<IonMobility>::coordinate(self)
    }
}

impl<T: CoordinateLike<C>, C> CoordinateLike<C> for &T {
    fn coordinate(&self) -> f64 {
        (*self).coordinate()
    }
}

impl<T: CoordinateLike<C>, C> CoordinateLike<C> for &mut T {
    fn coordinate(&self) -> f64 {
        CoordinateLike::<C>::coordinate(*self)
    }
}

impl<T: CoordinateLikeMut<C>, C> CoordinateLikeMut<C> for &mut T {
    fn coordinate_mut(&mut self) -> &mut f64 {
        CoordinateLikeMut::<C>::coordinate_mut(*self)
    }
}

impl<T: CoordinateLike<Mass>> MassLocated for T {}
impl<T: CoordinateLike<MZ>> MZLocated for T {}

impl<T: CoordinateLike<Time>> TimeLocated for T {}
impl<T: CoordinateLike<IonMobility>> IonMobilityLocated for T {}

/// A type alias for the index in an [`IndexedCoordinate`] structure
pub type IndexType = u32;

/// Indicate that an object may be indexed by coordinate system `T`
pub trait IndexedCoordinate<T>: CoordinateLike<T> {
    fn get_index(&self) -> IndexType;
    fn set_index(&mut self, index: IndexType);
}

impl<T: IndexedCoordinate<C>, C> IndexedCoordinate<C> for &T {
    fn get_index(&self) -> IndexType {
        (*self).get_index()
    }

    fn set_index(&mut self, _index: IndexType) {}
}

impl<T: IndexedCoordinate<C>, C> IndexedCoordinate<C> for &mut T {
    fn get_index(&self) -> IndexType {
        (**self).get_index()
    }

    fn set_index(&mut self, index: IndexType) {
        (**self).set_index(index)
    }
}

pub(crate) fn _isclose<T>(x: T, y: T, rtol: T, atol: T) -> bool
where
    T: Float,
{
    (x - y).abs() <= (atol + rtol * y.abs())
}

pub(crate) fn isclose<T>(x: T, y: T) -> bool
where
    T: Float + FromPrimitive,
{
    _isclose(x, y, T::from_f64(1e-5).unwrap(), T::from_f64(1e-8).unwrap())
}

pub trait HasProximity : PartialEq + PartialOrd + Copy {
    fn is_close(&self, other: &Self) -> bool {
        self == other
    }
}

macro_rules! impl_has_proximity {
    ($t:ty) => {
        impl $crate::coordinate::HasProximity for $t {
            fn is_close(&self, other: &Self) -> bool {
                isclose(*self, *other)
            }
        }
    };
}

impl_has_proximity!(f32);
impl_has_proximity!(f64);

macro_rules! impl_has_proximity_exact {
    ($t:ty) => {
        impl $crate::coordinate::HasProximity for $t {
            fn is_close(&self, other: &Self) -> bool {
                self == other
            }
        }
    };
}

impl_has_proximity_exact!(i8);
impl_has_proximity_exact!(i16);
impl_has_proximity_exact!(i32);
impl_has_proximity_exact!(i64);

impl<T: HasProximity> HasProximity for Option<T> {
    fn is_close(&self, other: &Self) -> bool {
        match (self, other) {
            (Some(x), Some(y)) => x.is_close(y),
            _ => false,
        }
    }
}

impl HasProximity for u8 {}
impl HasProximity for u16 {}
impl HasProximity for u32 {}
impl HasProximity for u64 {}

impl HasProximity for usize {}
impl HasProximity for isize {}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_is_close() {
        assert!(0.0.is_close(&0.0));
        assert!(Some(0.0).is_close(&Some(0.0)));
        assert!(!Some(0.0).is_close(&None));
        assert!(5.is_close(&5));
    }

    #[test]
    fn test_axes() {
        let dims = [Dimension::MZ(MZ()), Dimension::Mass(Mass()), Dimension::Time(Time()), Dimension::IonMobility(IonMobility()), Dimension::Dimensionless(Dimensionless())];
        for dim in dims {
            match dim {
                Dimension::MZ(x) => {
                    assert_eq!(x.name_of(), "m/z");
                    assert_eq!(dim.to_string(), "m/z");
                    assert_eq!(dim.minimum_value(), 0.0);
                    assert_eq!(dim.maximum_value(), f64::INFINITY);
                },
                Dimension::Mass(x) => {
                    assert_eq!(x.name_of(), "neutral mass");
                    assert_eq!(dim.to_string(), "neutral mass");
                    assert_eq!(dim.minimum_value(), 0.0);
                    assert_eq!(dim.maximum_value(), f64::INFINITY);
                },
                Dimension::IonMobility(x) => {
                    assert_eq!(x.name_of(), "ion mobility");
                    assert_eq!(dim.to_string(), "ion mobility");
                    assert_eq!(dim.minimum_value(), 0.0);
                    assert_eq!(dim.maximum_value(), f64::INFINITY);
                }
                Dimension::Time(x) => {
                    assert_eq!(x.name_of(), "time");
                    assert_eq!(dim.to_string(), "time");
                    assert_eq!(dim.minimum_value(), 0.0);
                    assert_eq!(dim.maximum_value(), f64::INFINITY);
                }
                Dimension::Dimensionless(x) => {
                    assert_eq!(x.name_of(), "");
                    assert_eq!(dim.to_string(), "");
                    assert_eq!(dim.minimum_value(), 0.0);
                    assert_eq!(dim.maximum_value(), f64::INFINITY);
                }
            }
        }
    }
}