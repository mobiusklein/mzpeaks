//! Collections of peaks that are ordered and searchable by a coordinate,
//! support fast access and are growable. While the main behaviors are provided
//! through the [`PeakCollection`] generic trait, a (generic) full implementation
//! is given by [`PeakSetVec`].
//!
//! The two basic peak types, [`CentroidPeak`](crate::peak::CentroidPeak) and
//! [`DeconvolutedPeak`](crate::peak::DeconvolutedPeak) have fully specified
//! aliases of [`PeakSetVec`], [`PeakSet`] and [`DeconvolutedPeakSet`],
//! respectively.
//!
//! [`PeakCollection`] can be searched by its specified coordinate space, with
//! [`PeakCollection::search`], [`PeakCollection::has_peak`], [`PeakCollection::all_peaks_for`],
//! and [`PeakCollection::between`].
//!
use std::fmt;
use std::iter::{Extend, FromIterator};
use std::marker;
use std::ops;

#[cfg(feature="serde")]
use serde::{Serialize, Deserialize};

use crate::mass_error::MassErrorType;

use crate::coordinate::{CoordinateLike, IndexType, IndexedCoordinate, Mass, MZ};
use crate::peak::{CentroidPeak, DeconvolutedPeak};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// When adding a peak to a [`PeakCollection`], indicate
/// whether the addition required re-indexing the whole
/// collection.
pub enum OrderUpdateEvent {
    /// No change was required, the new peak was added to the end
    TailAppend,
    /// The collection was re-indexed, the addition occurred in the middle
    /// of the collection
    InsertResorted,
}

/// A trait for an ordered container of mass spectral peaks. The trait
/// interoperates with [`CoordinateLike`].
pub trait PeakCollection<T: CoordinateLike<C>, C>: ops::Index<usize>
where
    <Self as ops::Index<usize>>::Output: CoordinateLike<C>,
{
    /// Add `peak` to the collection, maintaining sort order and peak
    /// indexing.
    fn push(&mut self, peak: T) -> OrderUpdateEvent;

    /// Sort the collection, updating the peak indexing.
    fn sort(&mut self);

    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Implement index access
    fn get_item(&self, i: usize) -> &T;
    fn get_slice(&self, i: ops::Range<usize>) -> &[T];

    /// Most basic method for coordinate search, find the
    /// index in this collection whose coordinate value is nearest
    /// to `query`. The return signature is identical to
    /// [`slice::binary_search_by`][slice::binary_search_by]
    fn search_by(&self, query: f64) -> Result<usize, usize>;

    #[inline]
    fn _closest_peak(
        &self,
        query: f64,
        error_tolerance: f64,
        i: usize,
        error_type: MassErrorType,
    ) -> Option<usize> {
        let mut j = i;
        let mut best = j;
        let mut best_err = error_type.call(self.get_item(j).coordinate(), query).abs();
        let n = self.len();
        // search backwards
        while j > 0 && j < n {
            let err = error_type.call(self.get_item(j).coordinate(), query).abs();
            if err < best_err && err < error_tolerance {
                best_err = err;
                best = j;
            } else if err > best_err {
                break;
            }
            j -= 1;
        }
        j = i;
        // search forwards
        while j < n {
            let err = error_type.call(self.get_item(j).coordinate(), query).abs();
            if err < best_err && err < error_tolerance {
                best_err = err;
                best = j;
            } else if err > best_err {
                break;
            }
            j += 1;
        }
        if best_err > error_tolerance {
            return None;
        }
        Some(best)
    }

    #[inline]
    /// Find the nearest index for `query` within `error_tolerance` units of `error_type` in
    /// this peak collection, or `None`.
    fn search(&self, query: f64, error_tolerance: f64, error_type: MassErrorType) -> Option<usize> {
        let lower_bound = error_type.lower_bound(query, error_tolerance);
        match self.search_by(lower_bound) {
            Ok(j) => self._closest_peak(query, error_tolerance, j, error_type),
            Err(j) => self._closest_peak(query, error_tolerance, j, error_type),
        }
    }

    #[inline]
    /// Return the peak nearest to `query` within `error_tolerance` units of `error_type` in
    /// this peak collection, or `None`.
    fn has_peak(&self, query: f64, error_tolerance: f64, error_type: MassErrorType) -> Option<&T> {
        return match self.search(query, error_tolerance, error_type) {
            Some(j) => Some(self.get_item(j)),
            None => None,
        };
    }

    #[inline]
    /// Return a slice containing all peaks between `low` and `high` coordinates within
    /// `error_tolerance` units of `error_type`.
    fn between(
        &self,
        low: f64,
        high: f64,
        error_tolerance: f64,
        error_type: MassErrorType,
    ) -> &[T] {
        let lower_bound = error_type.lower_bound(low, error_tolerance);
        let upper_bound = error_type.upper_bound(high, error_tolerance);

        let n = self.len();
        if n == 0 {
            return &self.get_slice(0..0);
        }

        let mut lower_index = match self.search_by(lower_bound) {
            Ok(j) => j,
            Err(j) => j,
        };

        let mut upper_index = match self.search_by(upper_bound) {
            Ok(j) => j,
            Err(j) => j,
        };

        if lower_index < n {
            if self[lower_index].coordinate() < lower_bound {
                lower_index += 1;
            }
        }

        if upper_index < n && upper_index > 0 {
            if self[upper_index].coordinate() > upper_bound {
                upper_index -= 1;
            }
        }

        if upper_index < n {
            upper_index += 1;
        }

        let subset = &self.get_slice(lower_index..upper_index);
        subset
    }

    #[inline]
    /// Find all peaks which could match `query` within `error_tolerance` units of `error_type`
    fn all_peaks_for(&self, query: f64, error_tolerance: f64, error_type: MassErrorType) -> &[T] {
        let lower_bound = error_type.lower_bound(query, error_tolerance);
        let upper_bound = error_type.upper_bound(query, error_tolerance);

        let n = self.len();
        if n == 0 {
            return &self.get_slice(0..0);
        }

        let mut lower_index = match self.search_by(lower_bound) {
            Ok(j) => j,
            Err(j) => j,
        };

        let checkpoint = lower_index;

        while lower_index < n && lower_index != 0 {
            if self[lower_index - 1].coordinate() > lower_bound {
                lower_index -= 1;
            } else {
                break;
            }
        }

        let mut upper_index = checkpoint;

        while upper_index < n {
            if self[upper_index + 1].coordinate() < upper_bound {
                upper_index += 1;
            } else {
                break;
            }
        }

        let v = self.get_item(lower_index).coordinate();
        if v <= lower_bound || v >= upper_bound {
            lower_index += 1;
        }
        let c = lower_index..upper_index + 1;
        return &self.get_slice(c);
    }
}

// An experiment to implement traits with macros
macro_rules! impl_slicing {
    ($t:ty, $($args:tt)+) => {

        impl<$($args)+> std::ops::Index<std::ops::Range<usize>> for $t {
            type Output = [<Self as std::ops::Index<usize>>::Output];

            fn index(&self, index: std::ops::Range<usize>) -> &Self::Output {
                return self.get_slice(index)
            }
        }

        impl<$($args)+> std::ops::Index<std::ops::RangeFrom<usize>> for $t {
            type Output = [<Self as std::ops::Index<usize>>::Output];

            fn index(&self, index: std::ops::RangeFrom<usize>) -> &Self::Output {
                let idx = std::ops::Range { start: index.start, end: self.len() };
                <Self as std::ops::Index<std::ops::Range<usize>>>::index(self, idx)
            }
        }

        impl<$($args)+> std::ops::Index<std::ops::RangeTo<usize>> for $t {
            type Output = [<Self as std::ops::Index<usize>>::Output];

            fn index(&self, index: std::ops::RangeTo<usize>) -> &Self::Output {
                let idx = std::ops::Range { start: 0, end: index.end };
                <Self as std::ops::Index<std::ops::Range<usize>>>::index(self, idx)
            }
        }

    };
}

/// Represent a sorted list of processed mass spectral peaks. It is a
/// concrete implementation of [`PeakCollection`] based on a [`Vec`].
#[derive(Default, Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PeakSetVec<P: IndexedCoordinate<C>, C> {
    pub peaks: Vec<P>,
    phantom: marker::PhantomData<C>,
}

impl<P: IndexedCoordinate<C>, C> PeakSetVec<P, C> {
    pub fn new(mut peaks: Vec<P>) -> Self {
        Self::_sort(&mut peaks);
        Self {
            peaks,
            phantom: marker::PhantomData,
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            peaks: Vec::with_capacity(capacity),
            phantom: marker::PhantomData,
        }
    }

    pub fn empty() -> Self {
        Self::with_capacity(0)
    }

    pub fn from_iter<V: Iterator>(peaks: V, sort: bool) -> Self
    where
        V: Iterator<Item = P>,
    {
        let peaks: Vec<P> = peaks.collect();
        if sort {
            Self::new(peaks)
        } else {
            Self::wrap(peaks)
        }
    }

    pub fn wrap(peaks: Vec<P>) -> Self {
        Self {
            peaks,
            phantom: marker::PhantomData,
        }
    }

    fn _sort(peaks: &mut Vec<P>) {
        peaks.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (i, p) in peaks.iter_mut().enumerate() {
            p.set_index(i as IndexType);
        }
    }

    pub fn iter(&self) -> PeakSetIter<P, C> {
        PeakSetIter::new(self)
    }

    pub fn iter_mut(&mut self) -> PeakSetIterMut<P, C> {
        PeakSetIterMut::new(self)
    }

    fn _push(&mut self, peak: P) {
        self.peaks.push(peak);
    }
}

impl<P: IndexedCoordinate<C>, C> PeakCollection<P, C> for PeakSetVec<P, C> {
    fn sort(&mut self) {
        Self::_sort(&mut self.peaks);
    }

    fn push(&mut self, peak: P) -> OrderUpdateEvent {
        let n = self.len();
        match self.peaks.last() {
            Some(p) => {
                if p <= &peak {
                    self.peaks.push(peak);
                    let fin = &mut self.peaks[n];
                    fin.set_index(n as IndexType);
                    OrderUpdateEvent::TailAppend
                } else {
                    self.peaks.push(peak);
                    self.sort();
                    OrderUpdateEvent::InsertResorted
                }
            }
            None => {
                self.peaks.push(peak);
                let fin = &mut self.peaks[n];
                fin.set_index(n as IndexType);
                OrderUpdateEvent::TailAppend
            }
        }
    }

    #[inline]
    fn len(&self) -> usize {
        self.peaks.len()
    }

    #[inline]
    fn get_item(&self, i: usize) -> &P {
        &self[i]
    }

    #[inline]
    fn get_slice(&self, i: ops::Range<usize>) -> &[P] {
        &self.peaks[i]
    }

    #[inline]
    fn search_by(&self, query: f64) -> Result<usize, usize> {
        self.peaks
            .binary_search_by(|peak| peak.coordinate().partial_cmp(&query).unwrap())
    }
}

impl<P: IndexedCoordinate<C>, C> ops::Index<usize> for PeakSetVec<P, C> {
    type Output = P;

    fn index(&self, i: usize) -> &Self::Output {
        &(self.peaks[i])
    }
}

impl<P: IndexedCoordinate<C>, C> ops::IndexMut<usize> for PeakSetVec<P, C> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.peaks[index]
    }
}

impl_slicing!(PeakSetVec<P, C>, P: IndexedCoordinate<C>, C);

impl<P: IndexedCoordinate<C>, C> fmt::Display for PeakSetVec<P, C> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "PeakSetVec(<{} Peaks>)", self.len())?;
        Ok(())
    }
}

impl<P: IndexedCoordinate<C>, C> PartialEq for PeakSetVec<P, C> {
    fn eq(&self, other: &Self) -> bool {
        if self.len() != other.len() {
            false
        } else {
            for (a, b) in self.iter().zip(other.iter()) {
                if a != b {
                    return false
                }
            }
            true
        }
    }
}

impl<P: IndexedCoordinate<C>, C> From<Vec<P>> for PeakSetVec<P, C> {
    fn from(v: Vec<P>) -> PeakSetVec<P, C> {
        PeakSetVec::wrap(v)
    }
}

impl<P: IndexedCoordinate<C>, C> FromIterator<P> for PeakSetVec<P, C> {
    fn from_iter<T>(iter: T) -> Self
    where
        T: IntoIterator<Item = P>,
    {
        let mut result = Self::empty();
        result.extend(iter);
        result
    }
}

impl<P: IndexedCoordinate<C>, C> Extend<P> for PeakSetVec<P, C> {
    fn extend<T>(&mut self, iter: T)
    where
        T: IntoIterator<Item = P>,
    {
        let mut last_coord = 0.0;
        let last_index = self.len();
        if let Some(last_peak) = self.peaks.last() {
            last_coord = last_peak.coordinate()
        }

        let mut valid = true;
        for p in iter {
            let coord = p.coordinate();
            if coord < last_coord {
                valid = false
            } else {
                last_coord = coord;
            }
            self._push(p);
        }
        if valid {
            for i in last_index..self.len() {
                self[i].set_index(i as IndexType);
            }
        } else {
            self.sort()
        }
    }
}

impl<'a, P: IndexedCoordinate<C>, C> IntoIterator for &'a PeakSetVec<P, C> {
    type Item = &'a P;
    type IntoIter = PeakSetIter<'a, P, C>;

    fn into_iter(self) -> Self::IntoIter {
        PeakSetIter::new(self)
    }
}

impl<'a, P: IndexedCoordinate<C>, C> IntoIterator for &'a mut PeakSetVec<P, C> {
    type Item = &'a mut P;
    type IntoIter = PeakSetIterMut<'a, P, C>;

    fn into_iter(self) -> Self::IntoIter {
        PeakSetIterMut::new(self)
    }
}

// ---- Iterators -----

/// Reference Iterator over [`PeakSetVec`]
pub struct PeakSetIter<'a, P, C> {
    iter: std::slice::Iter<'a, P>,
    phantom: marker::PhantomData<C>,
}

impl<'a, P: IndexedCoordinate<C>, C> PeakSetIter<'a, P, C> {
    fn new(peaks: &'a PeakSetVec<P, C>) -> PeakSetIter<'a, P, C> {
        PeakSetIter {
            iter: peaks.peaks.iter(),
            phantom: marker::PhantomData,
        }
    }
}

impl<'a, P: IndexedCoordinate<C>, C> Iterator for PeakSetIter<'a, P, C> {
    type Item = &'a P;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

/// Mutable Reference Iterator over [`PeakSetVec`]
pub struct PeakSetIterMut<'a, P, C> {
    iter: std::slice::IterMut<'a, P>,
    phantom: marker::PhantomData<C>,
}

impl<'a, P: IndexedCoordinate<C>, C> PeakSetIterMut<'a, P, C> {
    fn new(peaks: &'a mut PeakSetVec<P, C>) -> PeakSetIterMut<'a, P, C> {
        PeakSetIterMut {
            iter: peaks.peaks.iter_mut(),
            phantom: marker::PhantomData,
        }
    }
}

impl<'a, P, C> Iterator for PeakSetIterMut<'a, P, C> {
    type Item = &'a mut P;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

// ----- Specializations -----

/// A [`PeakSetVec`] of [`CentroidPeak`](crate::peak::CentroidPeak) items
/// ordered by m/z
pub type PeakSet = PeakSetVec<CentroidPeak, MZ>;

/// A [`PeakSetVec`] of [`DeconvolutedPeak`](crate::peak::DeconvolutedPeak) items
/// ordered by neutral mass
pub type DeconvolutedPeakSet = PeakSetVec<DeconvolutedPeak, Mass>;

/// A partial specialization of [`PeakSetVec`] that requires that the ordering
/// coordinate is m/z
pub type MZPeakSetType<P> = PeakSetVec<P, MZ>;

/// A partial specialization of [`PeakSetVec`] that requires that the ordering
/// coordinate is neutral mass
pub type MassPeakSetType<D> = PeakSetVec<D, Mass>;

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_data;

    #[test]
    fn test_sequence_behavior() {
        let peaks = test_data::read_peaks_from_file("./test/data/test.txt").unwrap();

        assert_eq!(peaks.len(), 485);
        assert!((peaks[0].mz - 231.3888).abs() < 1e-3);

        for (i, peak) in peaks.iter().enumerate() {
            if i > 0 {
                assert!(peak.mz > peaks[i - 1].mz);
            }
        }

        let part = peaks.search(773.4414, 0.01, MassErrorType::Absolute);
        assert_eq!(part.expect("Match peak"), 300);
        let part = peaks.has_peak(773.4414, 10.0, MassErrorType::PPM);
        assert_eq!(part.expect("Match peak").index, 300);

        let part = peaks.all_peaks_for(773.4414, 10.0, MassErrorType::PPM);
        assert_eq!(part.len(), 1);
        assert_eq!(part[0].index, 300);

        let part = peaks.search(773.4414, 50.0, MassErrorType::Absolute);
        assert_eq!(part.expect("Match peak"), 300);

        let part = peaks.all_peaks_for(736.637, 10.0, MassErrorType::PPM);
        assert_eq!(part.len(), 1);
        let part = peaks.all_peaks_for(736.237, 10.0, MassErrorType::PPM);
        assert_eq!(part.len(), 0);

        let q = 1221.639893;
        let block = peaks.all_peaks_for(q, 0.5, MassErrorType::Absolute);
        assert_eq!(block.len(), 1);
    }
}
