use std::{
    error::Error,
    fmt::Display,
    marker::PhantomData,
    num::ParseFloatError,
    ops::{Bound, Range, RangeBounds, RangeTo},
    str::FromStr,
};

use super::{CoordinateLike, HasProximity};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// An interval within a single dimension
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CoordinateRange<C> {
    pub start: Option<f64>,
    pub end: Option<f64>,
    #[cfg_attr(feature = "serde", serde(skip))]
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

    pub fn contains_raw(&self, x: &f64) -> bool {
        RangeBounds::<f64>::contains(&self, x)
    }

    pub fn overlaps<T: RangeBounds<f64>>(&self, interval: &T) -> bool {
        let interval_start = match interval.start_bound() {
            Bound::Included(x) => *x,
            Bound::Excluded(x) => *x,
            Bound::Unbounded => 0.0,
        };

        let interval_end = match interval.end_bound() {
            Bound::Included(y) => *y,
            Bound::Excluded(y) => *y,
            Bound::Unbounded => f64::INFINITY,
        };
        (self.end.unwrap_or(f64::INFINITY) >= interval_start
            && interval_end >= self.start.unwrap_or(0.0))
            || (self.end.is_close(&Some(interval_end))
                && self.start.is_close(&Some(interval_start)))
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

/** An inclusive interval over a single dimension
*/
pub trait Span1D {
    type DimType: HasProximity;

    fn start(&self) -> Self::DimType;
    fn end(&self) -> Self::DimType;

    fn contains(&self, i: &Self::DimType) -> bool {
        (self.start() <= *i && *i <= self.end()) || (self.start().is_close(i) || self.end().is_close(i))
    }

    fn is_close<T: Span1D<DimType = Self::DimType>>(&self, interval: &T) -> bool {
        self.start().is_close(&interval.start()) && self.end().is_close(&interval.end())
    }

    fn overlaps<T: Span1D<DimType = Self::DimType>>(&self, interval: &T) -> bool {
        (self.end() >= interval.start() && interval.end() >= self.start())
            || self.is_close(&interval)
    }

    fn is_contained_in_interval<T: Span1D<DimType = Self::DimType>>(&self, interval: &T) -> bool {
        (self.start() >= interval.start() && self.end() <= interval.end())
            || self.is_close(&interval)
    }

    fn contains_interval<T: Span1D<DimType = Self::DimType>>(&self, interval: &T) -> bool {
        (self.start() <= interval.start() && self.end() >= interval.end())
            || self.is_close(&interval)
    }
}

impl<T: Span1D> Span1D for &T {
    type DimType = T::DimType;

    fn start(&self) -> Self::DimType {
        (*self).start()
    }

    fn end(&self) -> Self::DimType {
        (*self).end()
    }
}

impl<C> Span1D for CoordinateRange<C> {
    type DimType = Option<f64>;

    fn start(&self) -> Self::DimType {
        self.start
    }

    fn end(&self) -> Self::DimType {
        self.end
    }
}

impl<T: HasProximity> Span1D for Range<T> {
    type DimType = T;

    fn start(&self) -> Self::DimType {
        self.start
    }

    fn end(&self) -> Self::DimType {
        self.end
    }
}

/// A basic [`Span1D`] implementation
#[derive(Debug, Default, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SimpleInterval<V: PartialOrd> {
    pub start: V,
    pub end: V,
}

impl<V: PartialOrd> SimpleInterval<V> {
    pub fn new(start: V, end: V) -> SimpleInterval<V> {
        SimpleInterval { start, end }
    }
}

impl<V: HasProximity> Span1D for SimpleInterval<V> {
    type DimType = V;

    fn start(&self) -> Self::DimType {
        self.start
    }

    fn end(&self) -> Self::DimType {
        self.end
    }
}

impl<V: PartialOrd> From<(V, V)> for SimpleInterval<V> {
    fn from(value: (V, V)) -> Self {
        Self::new(value.0, value.1)
    }
}

impl<V: PartialOrd> From<Range<V>> for SimpleInterval<V> {
    fn from(value: Range<V>) -> Self {
        Self::new(value.start, value.end)
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
        let start_t = if start_s.is_empty() {
            None
        } else {
            match start_s.parse() {
                Ok(val) => Some(val),
                Err(e) => return Err(CoordinateRangeParseError::MalformedStart(e)),
            }
        };
        let end_s = tokens.next().unwrap();
        let end_t = if end_s.is_empty() {
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
        let start = value.start.unwrap_or(0.0);
        let end = value.end.unwrap_or(f64::INFINITY);

        start..end
    }
}

#[cfg(test)]
mod test {
    use crate::Time;

    use super::*;

    #[test]
    fn test_conversion() {
        let time_range: CoordinateRange<Time> = (..5.0).into();
        assert_eq!(time_range.start(), None);
        assert_eq!(time_range.end(), Some(5.0));

        let time_range: CoordinateRange<Time> = (5.0..10.0).into();
        assert_eq!(time_range.start(), Some(5.0));
        assert_eq!(time_range.end(), Some(10.0));

        let time_range: CoordinateRange<Time> = ":5.0".parse().unwrap();
        assert_eq!(time_range.start(), None);
        assert_eq!(time_range.end(), Some(5.0));

        let time_range: CoordinateRange<Time> = "-5.0".parse().unwrap();
        assert_eq!(time_range.start(), None);
        assert_eq!(time_range.end(), Some(5.0));

        let time_range: CoordinateRange<Time> = " 5.0".parse().unwrap();
        assert_eq!(time_range.start(), None);
        assert_eq!(time_range.end(), Some(5.0));

        // let time_range: CoordinateRange<Time> = "5.0".parse().unwrap();
        // assert_eq!(time_range.start(), None);
        // assert_eq!(time_range.end(), Some(5.0));

        let time_range: CoordinateRange<Time> = (0.0, 5.0).into();
        assert_eq!(time_range.start(), Some(0.0));
        assert_eq!(time_range.end(), Some(5.0));

        let time_range: SimpleInterval<f64> = (0.0, 5.0).into();
        assert_eq!(time_range.start(), 0.0);
        assert_eq!(time_range.end(), 5.0);

        let time_range: SimpleInterval<f64> = (5.0..10.0).into();
        assert_eq!(time_range.start(), 5.0);
        assert_eq!(time_range.end(), 10.0);

        let empty = CoordinateRange::<Time>::default();
        assert_eq!(empty.start(), None);
        assert_eq!(empty.end(), None);

        assert_eq!(empty.start_bound(), Bound::Unbounded);
        assert_eq!(empty.end_bound(), Bound::Unbounded);
        let t: Range<f64> = empty.into();

        assert_eq!(t.start, 0.0);
        assert_eq!(t.end, f64::INFINITY);
    }

    #[test]
    fn test_interval() {
        let time_range: CoordinateRange<Time> = (..5.0).into();
        assert!(time_range.contains_raw(&2.5));
        assert!(!time_range.contains_raw(&7.5));

        let time_range2: CoordinateRange<Time> = (3.0..10.0).into();
        assert!(time_range.overlaps(&time_range2));

        assert!((5.0..7.0).is_contained_in_interval(&(3.0..10.0)));
        assert!((3.0..10.0).contains_interval(&(5.0..7.0)));
    }
}