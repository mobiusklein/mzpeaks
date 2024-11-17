use std::fmt::Display;

use num_traits::{Float, FromPrimitive};

/// An enum over the different coordinate planes
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Dimension {
    MZ(MZ),
    Mass(Mass),
    Time(Time),
    IonMobility(IonMobility),
    Dimensionless(Dimensionless),
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
}

impl Display for Dimension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq, PartialOrd)]
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

    fn dimension() -> Dimension;

    fn name() -> &'static str {
        Self::dimension().name()
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
