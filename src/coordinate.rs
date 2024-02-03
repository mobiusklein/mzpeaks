use std::{
    error::Error,
    fmt::Display,
    marker::PhantomData,
    num::ParseFloatError,
    ops::{Bound, Range, RangeBounds, RangeTo},
    str::FromStr,
};

#[derive(Default, Debug, Clone, Copy)]
/// The Mass To Charge Ratio (m/z) coordinate system
pub struct MZ {}

impl MZ {
    /// Access the m/z of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<MZ>>(inst: &T) -> f64 {
        CoordinateLike::<MZ>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy)]
/// The Mass coordinate system
pub struct Mass {}

impl Mass {
    /// Access the neutral mass of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<Mass>>(inst: &T) -> f64 {
        CoordinateLike::<Mass>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy)]
/// The Event Time coordinate system
pub struct Time {}
impl Time {
    /// Access the elapsed time of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<Time>>(inst: &T) -> f64 {
        CoordinateLike::<Time>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy)]
/// The Ion Mobility Time coordinate system
pub struct IonMobility {}
impl IonMobility {
    /// Access the ion mobility time unit of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<IonMobility>>(inst: &T) -> f64 {
        CoordinateLike::<IonMobility>::coordinate(inst)
    }
}

#[derive(Clone, Copy, Debug)]
pub enum CoordinateDimension {
    MZ(MZ),
    Mass(Mass),
    Time(Time),
    IonMobility(IonMobility),
}

/// Denote a type has a coordinate value on coordinate system `T`
pub trait CoordinateLike<T>: PartialOrd {
    /// The trait method for accessing the coordinate of the object on coordinate
    /// system `T`
    fn coordinate(&self) -> f64;
}

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

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct CoordinateRange<C> {
    pub start: Option<f64>,
    pub end: Option<f64>,
    coord: PhantomData<C>,
}

impl<C> CoordinateRange<C> {
    pub fn new(start: Option<f64>, end: Option<f64>) -> Self {
        Self {
            start,
            end,
            coord: PhantomData,
        }
    }

    pub fn contains<T: CoordinateLike<C>>(&self, point: &T) -> bool {
        let x = CoordinateLike::<C>::coordinate(point);
        RangeBounds::<f64>::contains(&self, &x)
    }
}

impl<C> Default for CoordinateRange<C> {
    fn default() -> Self {
        Self {
            start: None,
            end: None,
            coord: PhantomData,
        }
    }
}

#[derive(Debug)]
pub enum CoordinateRangeParseError {
    MalformedStart(ParseFloatError),
    MalformedEnd(ParseFloatError),
}

impl Display for CoordinateRangeParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CoordinateRangeParseError::MalformedStart(e) => {
                write!(f, "Failed to parse range start {e}")
            }
            CoordinateRangeParseError::MalformedEnd(e) => {
                write!(f, "Failed to parse range end {e}")
            }
        }
    }
}

impl Error for CoordinateRangeParseError {}

impl<C> FromStr for CoordinateRange<C> {
    type Err = CoordinateRangeParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut tokens = if s.contains(' ') {
            s.split(' ')
        } else if s.contains(':') {
            s.split(':')
        } else if s.contains('-') {
            s.split('-')
        } else {
            s.split(' ')
        };
        let start_s = tokens.next().unwrap();
        let start_t = if start_s == "" {
            None
        } else {
            match start_s.parse() {
                Ok(val) => Some(val),
                Err(e) => return Err(CoordinateRangeParseError::MalformedStart(e)),
            }
        };
        let end_s = tokens.next().unwrap();
        let end_t = if end_s == "" {
            None
        } else {
            match end_s.parse() {
                Ok(val) => Some(val),
                Err(e) => return Err(CoordinateRangeParseError::MalformedEnd(e)),
            }
        };
        Ok(CoordinateRange {
            start: start_t,
            end: end_t,
            coord: PhantomData,
        })
    }
}

impl<C> From<RangeTo<f64>> for CoordinateRange<C> {
    fn from(value: RangeTo<f64>) -> Self {
        Self::new(None, Some(value.end))
    }
}

impl<C> From<Range<f64>> for CoordinateRange<C> {
    fn from(value: Range<f64>) -> Self {
        Self::new(Some(value.start), Some(value.end))
    }
}

impl<C> RangeBounds<f64> for CoordinateRange<C> {
    fn start_bound(&self) -> Bound<&f64> {
        if let Some(start) = self.start.as_ref() {
            Bound::Included(start)
        } else {
            Bound::Unbounded
        }
    }

    fn end_bound(&self) -> Bound<&f64> {
        if let Some(end) = self.end.as_ref() {
            Bound::Included(end)
        } else {
            Bound::Unbounded
        }
    }
}

impl<C> RangeBounds<f64> for &CoordinateRange<C> {
    fn start_bound(&self) -> Bound<&f64> {
        (*self).start_bound()
    }

    fn end_bound(&self) -> Bound<&f64> {
        (*self).end_bound()
    }
}

impl<C> From<(f64, f64)> for CoordinateRange<C> {
    fn from(value: (f64, f64)) -> Self {
        Self::new(Some(value.0), Some(value.1))
    }
}

impl<C> From<CoordinateRange<C>> for Range<f64> {
    fn from(value: CoordinateRange<C>) -> Self {
        let start = match value.start {
            Some(x) => x,
            None => 0.0,
        };

        let end = match value.end {
            Some(x) => x,
            None => f64::INFINITY,
        };

        let r = start..end;
        r
    }
}