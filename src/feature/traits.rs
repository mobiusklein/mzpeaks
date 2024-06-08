use std::{collections::VecDeque, ops::RangeBounds};

use crate::{coordinate::CoordinateLike, CoordinateRange, IntensityMeasurement};

/// Represent an interval of time
pub trait TimeInterval<T> {
    /// The earliest time point recorded
    fn start_time(&self) -> Option<f64>;

    /// The latest time point recorded
    fn end_time(&self) -> Option<f64>;

    /// The time point where the feature reaches its greatest intensity
    fn apex_time(&self) -> Option<f64>;

    /// Integrate the feature in the time dimension
    fn area(&self) -> f32;

    /// Represent the [`TimeInterval`] into a [`CoordinateRange`]
    fn as_range(&self) -> CoordinateRange<T> {
        CoordinateRange::new(self.start_time(), self.end_time())
    }

    /// Check if a time point is spanned by [`TimeInterval`]
    fn spans(&self, time: f64) -> bool {
        let range = self.as_range();
        range.contains_raw(&time)
    }

    /// Return an iterator over the time dimension
    fn iter_time(&self) -> impl Iterator<Item = f64>;

    /// Find the position in the interval closest to the requested time
    /// and the magnitude of the error
    fn find_time(&self, time: f64) -> (Option<usize>, f64) {
        let mut best_i = None;
        let mut best_err = f64::INFINITY;
        for (i, t) in self.iter_time().enumerate() {
            let err = (t - time).abs();
            let err_abs = err.abs();
            if err_abs < best_err {
                best_i = Some(i);
                best_err = err_abs;
            }
            // We passed the upper end of the error range since time is sorted
            // so stop searching
            if err > best_err {
                break;
            }
        }
        (best_i, best_err)
    }
}

impl<'a, T, U: TimeInterval<T>> TimeInterval<T> for &'a U {
    fn start_time(&self) -> Option<f64> {
        (*self).start_time()
    }

    fn end_time(&self) -> Option<f64> {
        (*self).end_time()
    }

    fn apex_time(&self) -> Option<f64> {
        (*self).apex_time()
    }

    fn area(&self) -> f32 {
        (*self).area()
    }

    fn iter_time(&self) -> impl Iterator<Item = f64> {
        (*self).iter_time()
    }
}

/// An expansion of [`TimeInterval`] which provides a contiguous slice over the time dimension
pub trait TimeArray<T> : TimeInterval<T> {
    /// A slice over the complete time dimension
    fn time_view(&self) -> &[f64];
}


impl<'a, T, U: TimeArray<T>> TimeArray<T> for &'a U {
    fn time_view(&self) -> &[f64] {
        (*self).time_view()
    }
}

/// Represents something that is located at a constrained but varying coordinate system `X` over a
/// sequentially ordered dimension `Y` with an abundance measure at each time point.
pub trait FeatureLike<X, Y>: IntensityMeasurement + TimeInterval<Y> + CoordinateLike<X> {
    /// The number of points in the feature
    fn len(&self) -> usize;
    /// Create an iterator that yields (x, y, intensity) references
    fn iter(&self) -> impl Iterator<Item = (&f64, &f64, &f32)>;
    /// Check if the feature has any points in it
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get an immutable reference to feature data at a specified index
    fn at(&self, index: usize) -> Option<(f64, f64, f32)> {
        self.iter().nth(index).map(|(x, y, z)| (*x, *y, *z))
    }

    /// Get an immutable reference to feature data at a specified time.Analogous
    /// to combining [`TimeInterval::find_time`] with [`FeatureLike::at`]
    fn at_time(&self, time: f64) -> Option<(f64, f64, f32)> {
        if let (Some(ix), _) = self.find_time(time) {
            self.at(ix)
        } else {
            None
        }
    }
}

impl<'a, X, Y, T: FeatureLike<X, Y> + TimeInterval<Y>> FeatureLike<X, Y> for &'a T {
    fn len(&self) -> usize {
        (*self).len()
    }

    fn iter(&self) -> impl Iterator<Item = (&f64, &f64, &f32)> {
        (*self).iter()
    }
}

/// A [`FeatureLike`] type that is also mutable
pub trait FeatureLikeMut<X, Y>: FeatureLike<X, Y> {
    /// Create an iterator that yields (x, y, intensity) mutable references
    fn iter_mut(&mut self) -> impl Iterator<Item = (&mut f64, &mut f64, &mut f32)>;
    /// Add a new peak-like reference to the feature at a given y "time" coordinate. If the "time"
    /// is not in sorted order, it should automatically re-sort.
    fn push<T: CoordinateLike<X> + IntensityMeasurement>(&mut self, pt: &T, time: f64);
    /// As [`FeatureLikeMut::push`], but instead add raw values instead of deriving them from
    /// a peak-like reference.
    fn push_raw(&mut self, x: f64, y: f64, z: f32);

    /// Get a mutable reference to feature data at a specified index
    fn at_mut(&mut self, index: usize) -> Option<(&mut f64, &mut f64, &mut f32)> {
        self.iter_mut().nth(index)
    }

    /// Get a mutable reference to feature data at a specified time. Analogous
    /// to combining [`TimeInterval::find_time`] with [`FeatureLikeMut::at_mut`]
    fn at_time_mut(&mut self, time: f64) -> Option<(&mut f64, &mut f64, &mut f32)> {
        if let (Some(ix), _) = self.find_time(time) {
            self.at_mut(ix)
        } else {
            None
        }
    }
}

pub(crate) trait CoArrayOps {
    fn weighted_average(&self, x: &[f64], w: &[f32]) -> f64 {
        let (acc, norm) = x
            .iter()
            .zip(w.iter())
            .fold((0.0, 0.0), |(acc, norm), (x, z)| {
                let norm = norm + (*z as f64);
                let acc = acc + *x * (*z as f64);
                (acc, norm)
            });
        if norm == 0.0 {
            return 0.0;
        }
        acc / norm
    }

    fn trapezoid_integrate(&self, y: &[f64], w: &[f32]) -> f32 {
        let mut it = y.iter().zip(w.iter());
        if let Some((first_y, first_z)) = it.next() {
            let (_y, _z, acc) =
                it.fold((first_y, first_z, 0.0), |(last_y, last_z, acc), (y, z)| {
                    let step = (last_z + z) / 2.0;
                    let dy = y - last_y;
                    (y, z, acc + (step as f64 * dy))
                });
            acc as f32
        } else {
            0.0
        }
    }

    fn idxmax(&self, w: &[f32]) -> Option<usize> {
        let pt = w
            .iter()
            .enumerate()
            .reduce(|(best_i, best), (current_i, current)| {
                if *current > *best {
                    (current_i, current)
                } else {
                    (best_i, best)
                }
            });
        pt.map(|(i, _)| i)
    }

    fn apex_of(&self, x: &[f64], w: &[f32]) -> Option<f64> {
        self.idxmax(w).and_then(|i| x.get(i).copied())
    }
}

/// A trait to split features at a given point or points, producing either copies
/// or borrows of the same feature data.
pub trait SplittableFeatureLike<'a, X, Y>: FeatureLike<X, Y> {
    type ViewType: FeatureLike<X, Y>;

    /// Split the feature at `index`, segmenting before and after it. The
    /// position at `index` should be retained in the second segment.
    fn split_at(&'a self, index: usize) -> (Self::ViewType, Self::ViewType);

    /// Split the feature at `time`.
    ///
    /// This may be implemented as simply `self.split_at(self.find_time(time).0.unwrap())`
    fn split_at_time(&'a self, point: f64) -> (Self::ViewType, Self::ViewType);

    /// Select the positions given by the range of positions in `bounds`
    fn slice<I: RangeBounds<usize> + Clone>(&'a self, bounds: I) -> Self::ViewType;

    /// Retain all positions where `mask` is `true`, keeping contiguous positions
    /// part of the same feature segments.
    ///
    /// This requires `mask.len() == self.len()`. It will panic otherwise.
    fn split_mask(&'a self, mask: &[bool]) -> Vec<Self::ViewType> {
        assert_eq!(
            self.len(),
            mask.len(),
            "Feature length ({}) != Mask length ({})",
            self.len(),
            mask.len()
        );
        let (mut spans, last_state, start) = mask.iter().copied().enumerate().fold(
            (VecDeque::<(usize, usize)>::new(), false, None),
            |(mut spans, last_state, start), (i, mark)| {
                if !last_state && mark {
                    (spans, mark, Some(i))
                } else if !mark && last_state {
                    spans.push_back((start.unwrap(), i));
                    (spans, mark, None)
                } else {
                    (spans, mark, start)
                }
            },
        );

        if last_state {
            spans.push_back((start.unwrap(), self.len()));
        }

        let mut segments = Vec::new();
        for (start, end) in spans {
            segments.push(self.slice(start..end));
        }
        segments
    }

    /// Given a function `f` that takes successive pairs of points of `(dim1, dim2, intensity)`
    /// and returns a `bool`, mask all out all positions where `f` returns `true`.
    ///
    /// See [`SplittableFeatureLike::split_mask`]
    fn split_when<F>(&'a self, mut f: F) -> Vec<Self::ViewType>
    where
        F: FnMut((f64, f64, f32), (f64, f64, f32)) -> bool,
    {
        let mut prev = self.at(0).unwrap_or_default();
        let mut mask = Vec::with_capacity(self.len());
        for cur in self.iter().map(|(x, y, z)| (*x, *y, *z)) {
            mask.push(!f(prev, cur));
            prev = cur;
        }
        self.split_mask(&mask)
    }

    /// Split the feature when there is a gap of size `max_gap_size` or more
    /// in the time dimension
    fn split_sparse(&'a self, max_gap_size: f64) -> Vec<Self::ViewType> {
        self.split_when(|(_, prev_time, _), (_, cur_time, _)| (cur_time - prev_time) > max_gap_size)
    }
}
