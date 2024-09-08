//! Collections of peaks that are ordered and searchable by a coordinate,
//! support fast access and are growable. While the main behaviors are provided
//! through the [`PeakCollection`] generic trait, a (generic) full implementation
//! is given by [`PeakSetVec`].
//!
//! The two basic peak types, [`CentroidPeak`] and
//! [`DeconvolutedPeak`] have fully specified
//! aliases of [`PeakSetVec`], [`PeakSet`] and [`DeconvolutedPeakSet`],
//! respectively.
//!
//! [`PeakCollection`] can be searched by its specified coordinate space, with
//! [`PeakCollection::search`], [`PeakCollection::has_peak`], [`PeakCollection::all_peaks_for`],
//! and [`PeakCollection::between`].
//!
use std::fmt::{self, Display};
use std::iter::{Extend, FromIterator};
use std::marker::{self, PhantomData};
use std::ops::{self, Range};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::mass_error::Tolerance;

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
/// builds upon [`CoordinateLike`].
pub trait PeakCollection<T: CoordinateLike<C>, C>: ops::Index<usize, Output = T> {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Implement index access
    fn get_item(&self, i: usize) -> &T;
    fn get_slice(&self, i: ops::Range<usize>) -> &[T];

    fn iter(&self) -> impl Iterator<Item = &T>
    where
        T: 'static;

    /// Most basic method for coordinate search, find the
    /// index in this collection whose coordinate value is nearest
    /// to `query`. The return signature is identical to
    /// [`slice::binary_search_by`][slice::binary_search_by]
    fn search_by(&self, query: f64) -> Result<usize, usize>;

    #[inline]
    fn _closest_peak(&self, query: f64, error_tolerance: Tolerance, i: usize) -> Option<usize> {
        if i >= self.len() {
            return None;
        }
        let mut j = i;
        let mut best = j;
        let mut best_err = error_tolerance
            .call(self.get_item(j).coordinate(), query)
            .abs();
        let n = self.len();
        let tol = error_tolerance.tol();
        // search backwards
        while j > 0 && j < n {
            let err = error_tolerance
                .call(self.get_item(j).coordinate(), query)
                .abs();
            if err < best_err && err < tol {
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
            let err = error_tolerance
                .call(self.get_item(j).coordinate(), query)
                .abs();
            if err < best_err && err < tol {
                best_err = err;
                best = j;
            } else if err > best_err {
                break;
            }
            j += 1;
        }
        if best_err > tol {
            return None;
        }
        Some(best)
    }

    #[inline]
    /// Find the nearest index for `query` within `error_tolerance` in
    /// this peak collection, or `None`.
    fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        let lower_bound = error_tolerance.bounds(query).0;

        match self.search_by(lower_bound) {
            Ok(j) => self._closest_peak(query, error_tolerance, j),
            Err(j) => self._closest_peak(query, error_tolerance, j),
        }
    }

    #[inline]
    /// Return the peak nearest to `query` within `error_tolerance` in
    /// this peak collection, or `None`.
    fn has_peak(&self, query: f64, error_tolerance: Tolerance) -> Option<&T> {
        return match self.search(query, error_tolerance) {
            Some(j) => Some(self.get_item(j)),
            None => None,
        };
    }

    #[inline]
    /// Return a slice containing all peaks between `low` and `high` coordinates within
    /// `error_tolerance`.
    fn between(&self, low: f64, high: f64, error_tolerance: Tolerance) -> &[T] {
        let lower_bound = error_tolerance.bounds(low).0;
        let upper_bound = error_tolerance.bounds(high).1;

        let n = self.len();
        if n == 0 {
            return self.get_slice(0..0);
        }

        let mut lower_index = match self.search_by(lower_bound) {
            Ok(j) => j,
            Err(j) => j,
        };

        let mut upper_index = match self.search_by(upper_bound) {
            Ok(j) => j,
            Err(j) => j,
        };

        if lower_index < n && self[lower_index].coordinate() < lower_bound {
            lower_index += 1;
        }

        if upper_index < n && upper_index > 0 && self[upper_index].coordinate() > upper_bound {
            upper_index -= 1;
        }

        if upper_index < n {
            upper_index += 1;
        }

        if lower_index >= n {
            return self.get_slice(0..0);
        }

        let subset = self.get_slice(lower_index..upper_index);
        subset
    }

    #[inline]
    /// Find all peaks which could match `query` within `error_tolerance` units
    fn all_peaks_for(&self, query: f64, error_tolerance: Tolerance) -> &[T] {
        let (lower_bound, upper_bound) = error_tolerance.bounds(query);

        let n = self.len();
        if n == 0 {
            return self.get_slice(0..0);
        }

        let mut lower_index = match self.search_by(lower_bound) {
            Ok(j) => j,
            Err(j) => j.min(n - 1),
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

        while upper_index < n - 1 {
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
        return self.get_slice(c);
    }

    /// Given a **sorted** list of `query` values that are on the same coordinate system `C`, find
    /// the pairs of indices from `query` to `self` within `error_tolerance` units of error.
    ///
    /// This is more efficient than using [`PeakCollection::all_peaks_for`] on each query individually,
    /// as it consumes $`O(m\log_2{n})`$, where `self` contains `n` items and `queries` contains `m` iems.
    /// [`PeakCollection::search_sorted_all_indices`] merges two sorted lists, which has $`O(n + m)`$.
    fn search_sorted_all_indices<Q: CoordinateLike<C>>(
        &self,
        queries: &[Q],
        error_tolerance: Tolerance,
    ) -> Vec<(usize, usize)> {
        let mut checkpoint: usize = 0;
        let mut pairs: Vec<_> = Vec::new();
        let n = self.len();
        for (query_i, query) in queries.iter().enumerate() {
            let (lb, ub) = error_tolerance.bounds(query.coordinate());
            for (p, ref_i) in self.get_slice(checkpoint..n).iter().zip(checkpoint..n) {
                let coord = p.coordinate();
                if coord < lb {
                    checkpoint = ref_i;
                } else if coord > lb && coord < ub {
                    pairs.push((query_i, ref_i))
                } else if coord > ub {
                    break;
                }
            }
        }
        pairs
    }

    /// This alternative implementation of [`PeakCollection::search_sorted_all_indices`] that uses
    /// an iterator instead of creating intermediate storage. It should be equivalent in behavior
    /// but touch the heap less.
    fn search_sorted_all_indices_iter<'a, Q: CoordinateLike<C>>(
        &'a self,
        queries: &'a [Q],
        error_tolerance: Tolerance,
    ) -> SearchSortedIter<'a, C, T, Self, Q>
    where
        Self: Sized,
    {
        SearchSortedIter::new(self, queries, error_tolerance)
    }
}

#[derive(Debug)]
pub struct SearchSortedIter<
    'a,
    C,
    T: CoordinateLike<C>,
    P: PeakCollection<T, C>,
    Q: CoordinateLike<C>,
> {
    peaks: &'a P,
    queries: &'a [Q],
    error_tolerance: Tolerance,
    checkpoint: usize,
    query_i: usize,
    query_n: usize,
    query_lower_bound: f64,
    query_upper_bound: f64,
    self_n: usize,
    query_hit_range: <Range<usize> as IntoIterator>::IntoIter,
    _c: PhantomData<C>,
    _t: PhantomData<T>,
}

impl<'a, C, T: CoordinateLike<C>, P: PeakCollection<T, C>, Q: CoordinateLike<C>> Iterator
    for SearchSortedIter<'a, C, T, P, Q>
{
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        self.produce()
    }
}

impl<'a, C, T: CoordinateLike<C>, P: PeakCollection<T, C>, Q: CoordinateLike<C>>
    SearchSortedIter<'a, C, T, P, Q>
{
    pub fn new(peaks: &'a P, queries: &'a [Q], error_tolerance: Tolerance) -> Self {
        let self_n = peaks.len();
        let query_n = queries.len();

        let mut this = Self {
            peaks,
            queries,
            error_tolerance,
            checkpoint: 0,
            self_n,
            query_i: 0,
            query_n,
            query_lower_bound: 0.0,
            query_upper_bound: 0.0,
            query_hit_range: (0..0).into_iter(),
            _c: PhantomData,
            _t: PhantomData,
        };
        if query_n > 0 {
            this.update_bounds();
            this.query_hit_range = this.seek_next_matches().into_iter();
        }
        this
    }

    fn update_bounds(&mut self) {
        (self.query_lower_bound, self.query_upper_bound) = self
            .error_tolerance
            .bounds(self.queries[self.query_i].coordinate());
    }

    fn finalize_query(&mut self) -> bool {
        if self.query_i < self.query_n.saturating_sub(1) {
            self.query_i += 1;
            self.update_bounds();
            self.query_hit_range = self.seek_next_matches().into_iter();
            true
        } else {
            false
        }
    }

    fn seek_next_matches(&mut self) -> Range<usize> {
        let i = self.checkpoint;
        let n = self.self_n;
        let mut ref_start = None;
        let mut ref_end = 0;
        for (p, ref_i) in self.peaks.get_slice(i..n).iter().zip(i..n) {
            let coord = p.coordinate();
            if coord < self.query_lower_bound {
                self.checkpoint = ref_i;
            } else if coord >= self.query_lower_bound
                && coord < self.query_upper_bound
                && ref_start.is_none()
            {
                ref_start = Some(ref_i);
            } else if coord > self.query_upper_bound {
                ref_end = ref_i;
                break;
            }
        }
        ref_start.unwrap_or(ref_end)..ref_end
    }

    fn produce(&mut self) -> Option<(usize, usize)> {
        let pair = self.query_hit_range.next().map(|i| (self.query_i, i));
        if pair.is_none() {
            if !self.finalize_query() {
                return None;
            }
            self.query_hit_range.next().map(|i| (self.query_i, i))
        } else {
            pair
        }
    }
}

/// A [`PeakCollection`] that can have additional peaks added to it.
pub trait PeakCollectionMut<T: CoordinateLike<C>, C>: PeakCollection<T, C> {
    /// Add `peak` to the collection, maintaining sort order and peak
    /// indexing.
    fn push(&mut self, peak: T) -> OrderUpdateEvent;

    /// Sort the collection, updating the peak indexing.
    fn sort(&mut self);
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

        impl<$($args)+> std::ops::Index<std::ops::RangeFull> for $t {
            type Output = [<Self as std::ops::Index<usize>>::Output];

            fn index(&self, _: std::ops::RangeFull) -> &Self::Output {
                <Self as std::ops::Index<std::ops::Range<usize>>>::index(self, 0..self.len())
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
    /// Create a new [`PeakSetVec`] from an existing `Vec<P>` and sorts
    /// the newly created structure to ensure it is ordered by coordinate `C`
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

    pub fn from_iter<I: Iterator<Item = P>>(peaks: I, sort: bool) -> Self {
        let peaks: Vec<P> = peaks.collect();
        if sort {
            Self::new(peaks)
        } else {
            Self::wrap(peaks)
        }
    }

    /// Create a new [`PeakSetVec`] from an existing `Vec<P>`, but does not actively
    /// sort the collection. It is up to the caller to ensure that the provided `Vec`
    /// is sorted or that it will be sorted prior to any of its search functionality
    /// is used.
    pub fn wrap(peaks: Vec<P>) -> Self {
        Self {
            peaks,
            phantom: marker::PhantomData,
        }
    }

    fn _sort(peaks: &mut [P]) {
        peaks.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (i, p) in peaks.iter_mut().enumerate() {
            p.set_index(i as IndexType);
        }
    }

    /// Iterate over references to peaks
    pub fn iter(&self) -> PeakSetIter<P> {
        self.peaks.iter()
    }

    pub fn first(&self) -> Option<&P> {
        self.peaks.first()
    }

    pub fn last(&self) -> Option<&P> {
        self.peaks.last()
    }

    /// Iterate over mutable references to peaks
    pub fn iter_mut(&mut self) -> PeakSetIterMut<P> {
        self.peaks.iter_mut()
    }

    fn _push(&mut self, peak: P) {
        self.peaks.push(peak);
    }

    pub fn as_slice(&self) -> &[P] {
        self.peaks.as_slice()
    }

    pub fn as_mut_slice(&mut self) -> &mut [P] {
        self.peaks.as_mut_slice()
    }
}

impl<P: IndexedCoordinate<C>, C> PeakCollectionMut<P, C> for PeakSetVec<P, C> {
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
}

impl<P: IndexedCoordinate<C>, C> PeakCollection<P, C> for PeakSetVec<P, C> {
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

    fn iter(&self) -> impl Iterator<Item = &P> {
        self.iter()
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
                    return false;
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

impl<P: IndexedCoordinate<C>, C> IntoIterator for PeakSetVec<P, C> {
    type Item = P;
    type IntoIter = std::vec::IntoIter<P>;

    fn into_iter(self) -> Self::IntoIter {
        self.peaks.into_iter()
    }
}

impl<'a, P: IndexedCoordinate<C>, C> IntoIterator for &'a PeakSetVec<P, C> {
    type Item = &'a P;
    type IntoIter = PeakSetIter<'a, P>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a, P: IndexedCoordinate<C>, C> IntoIterator for &'a mut PeakSetVec<P, C> {
    type Item = &'a mut P;
    type IntoIter = PeakSetIterMut<'a, P>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

// ---- Iterators -----

pub type PeakSetIter<'a, P> = std::slice::Iter<'a, P>;
pub type PeakSetIterMut<'a, P> = std::slice::IterMut<'a, P>;

// ----- Specializations -----

/// A [`PeakSetVec`] of [`CentroidPeak`] items
/// ordered by m/z
pub type PeakSet = PeakSetVec<CentroidPeak, MZ>;

/// A [`PeakSetVec`] of [`DeconvolutedPeak`] items
/// ordered by neutral mass
pub type DeconvolutedPeakSet = PeakSetVec<DeconvolutedPeak, Mass>;

/// A partial specialization of [`PeakSetVec`] that requires that the ordering
/// coordinate is m/z
pub type MZPeakSetType<P> = PeakSetVec<P, MZ>;

/// A partial specialization of [`PeakSetVec`] that requires that the ordering
/// coordinate is neutral mass
pub type MassPeakSetType<D> = PeakSetVec<D, Mass>;

/// A borrowed view of a peak list that assumes that it is sorted by its coordinate
/// dimension ahead of time. Unlike [`PeakSetVec`], this collection does not attempt
/// to sort or re-index the peaks it contains.
#[derive(Debug, Clone, Copy)]
pub struct PeakSetView<'a, P: IndexedCoordinate<C>, C> {
    peaks: &'a [P],
    _c: PhantomData<C>,
}

impl<'a, P: IndexedCoordinate<C>, C> PeakSetView<'a, P, C> {
    /// # Safety
    /// This construction does not enforce a sorting on the elements of the peak
    /// list.
    pub unsafe fn wrap(peaks: &'a [P]) -> Self {
        Self {
            peaks,
            _c: PhantomData,
        }
    }

    pub fn first(&self) -> Option<&P> {
        self.peaks.first()
    }

    pub fn last(&self) -> Option<&P> {
        self.peaks.last()
    }

    fn is_sorted(peaks: &[P]) -> bool {
        let mut prev = 0.0;
        for p in peaks {
            let c = p.coordinate();
            if prev <= c {
                prev = c;
            } else {
                return false;
            }
        }
        true
    }

    pub fn iter(&self) -> PeakSetIter<'a, P> {
        self.peaks.iter()
    }
}

impl<'a, P: IndexedCoordinate<C>, C> ops::Index<usize> for PeakSetView<'a, P, C> {
    type Output = P;

    fn index(&self, i: usize) -> &Self::Output {
        &(self.peaks[i])
    }
}

impl_slicing!(PeakSetView<'a, P, C>, 'a, P: IndexedCoordinate<C>, C);

impl<'a, P: IndexedCoordinate<C>, C> PeakCollection<P, C> for PeakSetView<'a, P, C> {
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

    fn iter(&self) -> impl Iterator<Item = &P> {
        self.iter()
    }
}

impl<'a, P: IndexedCoordinate<C>, C> IntoIterator for PeakSetView<'a, P, C> {
    type Item = &'a P;

    type IntoIter = PeakSetIter<'a, P>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a, P: IndexedCoordinate<C>, C> IntoIterator for &PeakSetView<'a, P, C> {
    type Item = &'a P;

    type IntoIter = PeakSetIter<'a, P>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViewConversionError {
    Unsorted,
}

impl Display for ViewConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for ViewConversionError {}

impl<'a, P: IndexedCoordinate<C>, C> TryFrom<&'a [P]> for PeakSetView<'a, P, C> {
    type Error = ViewConversionError;

    fn try_from(value: &'a [P]) -> Result<Self, Self::Error> {
        if PeakSetView::is_sorted(value) {
            Ok(unsafe { PeakSetView::wrap(value) })
        } else {
            Err(ViewConversionError::Unsorted)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{peak::MZPoint, test_data};

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

        let part = peaks.search(773.4414, Tolerance::Da(0.01));
        assert_eq!(part.expect("Match peak"), 300);
        let part = peaks.has_peak(773.4414, Tolerance::PPM(10.0));
        assert_eq!(part.expect("Match peak").index, 300);

        let part = peaks.all_peaks_for(773.4414, Tolerance::PPM(10.0));
        assert_eq!(part.len(), 1);
        assert_eq!(part[0].index, 300);

        let part = peaks.search(773.4414, Tolerance::Da(50.0));
        assert_eq!(part.expect("Match peak"), 300);

        let part = peaks.all_peaks_for(736.637, Tolerance::PPM(10.0));
        assert_eq!(part.len(), 1);
        let part = peaks.all_peaks_for(736.237, Tolerance::PPM(10.0));
        assert_eq!(part.len(), 0);

        let q = 1221.639893;
        let block = peaks.all_peaks_for(q, Tolerance::Da(0.5));
        assert_eq!(block.len(), 1);

        let q = 2000.0;
        let block = peaks.all_peaks_for(q, Tolerance::PPM(10.));
        assert_eq!(block.len(), 0);

        let q = -2000.0;
        let block = peaks.all_peaks_for(q, Tolerance::PPM(10.));
        assert_eq!(block.len(), 0);

        let block = peaks.between(-2000f64, 2000f64, Tolerance::PPM(10.0));
        assert_eq!(block.len(), peaks.len());

        let block = peaks.between(0.0, 2000f64, Tolerance::PPM(10.0));
        assert_eq!(block.len(), peaks.len());

        let block = peaks.between(1313.0, 1316.0, "10.0ppm".parse().unwrap());
        assert_eq!(block.len(), 3);
    }

    #[test]
    fn test_edgecases() {
        let peaks = PeakSet::new(vec![CentroidPeak::new(500.0, 2., 0)]);

        let p = peaks.has_peak(500.0, Tolerance::Da(1.0));
        assert!(p.is_some());

        let p = peaks.all_peaks_for(500.0, Tolerance::Da(1.0));
        assert!(p.len() == 1);

        let peaks = PeakSet::new(vec![]);

        let p = peaks.has_peak(500.0, Tolerance::Da(1.0));
        assert!(p.is_none());

        let p = peaks.all_peaks_for(500.0, Tolerance::Da(1.0));
        assert!(p.len() == 0);
    }

    #[test]
    fn test_all_peaks_for() -> std::io::Result<()> {
        let peaks = test_data::read_peaks_from_file("./test/data/test.txt")?;

        let queries = vec![
            MZPoint::new(262.2675, 1000.0),
            MZPoint::new(566.3623, 1210.0),
            MZPoint::new(645.1916, 20.0),
        ];

        let indices = peaks.search_sorted_all_indices(&queries, Tolerance::PPM(20.0));
        let it = SearchSortedIter::new(&peaks, &queries, Tolerance::PPM(20.0));
        let iter_indices: Vec<_> = it.collect();

        assert_eq!(indices, iter_indices);
        Ok(())
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_serialize() -> std::io::Result<()> {
        use serde_json;
        use std::io;
        use std::io::prelude::*;

        let peaks = test_data::read_peaks_from_file("./test/data/test.txt")?;
        let mut buff = Vec::new();
        let buffer_writer = io::Cursor::new(&mut buff);
        let mut writer = io::BufWriter::new(buffer_writer);
        serde_json::to_writer_pretty(&mut writer, &peaks)?;
        writer.flush()?;
        let view = String::from_utf8_lossy(writer.get_ref().get_ref());
        let dup: PeakSet = serde_json::from_str(&view)?;
        assert_eq!(peaks.len(), dup.len());
        assert_eq!(peaks, dup);
        Ok(())
    }
}
