use core::slice;
use std::cmp::Ordering;
use std::marker::PhantomData;
use std::ops::{Bound, RangeBounds};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{coordinate::CoordinateLike, IntensityMeasurement};
use crate::{CentroidPeak, IonMobility, Time, MZ};

use super::traits::{CoArrayOps, FeatureLike, FeatureLikeMut, SplittableFeatureLike};
use super::util::{NonNan, EMPTY_X, EMPTY_Y, EMPTY_Z};
use super::{TimeArray, TimeInterval};

/// A basic implementation of [`FeatureLike`] and [`FeatureLikeMut`]
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Feature<X, Y> {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f32>,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<X, Y> CoArrayOps for Feature<X, Y> {}

impl<X, Y> Feature<X, Y> {
    pub fn new(x: Vec<f64>, y: Vec<f64>, z: Vec<f32>) -> Self {
        assert_eq!(x.len(), y.len());
        assert_eq!(x.len(), z.len());
        Self {
            x,
            y,
            z,
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    /// Create an empty [`Feature`], ready to be extended.
    pub fn empty() -> Self {
        Self {
            x: Vec::new(),
            y: Vec::new(),
            z: Vec::new(),
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    /// Create an empty [`Feature`] with pre-allocated capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            x: Vec::with_capacity(capacity),
            y: Vec::with_capacity(capacity),
            z: Vec::with_capacity(capacity),
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    /// Compute a weighted average over the X dimension
    pub(crate) fn coordinate_x(&self) -> f64 {
        self.weighted_average(&self.x, &self.z)
    }

    #[allow(unused)]
    /// Compute a weighted average over the Y dimension
    pub(crate) fn coordinate_y(&self) -> f64 {
        self.weighted_average(&self.y, &self.z)
    }

    /// The number of points in the feature
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Find the time where the feature achieves its maximum abundance
    pub(crate) fn apex_y(&self) -> Option<f64> {
        self.apex_of(&self.y, &self.z)
    }

    /// Sort the feature by the Y dimension
    pub(crate) fn sort_by_y(&mut self) {
        let mut indices: Vec<_> = (0..self.len()).collect();
        indices.sort_by_key(|i| NonNan::new(self.y[*i]));

        let mut xtmp: Vec<f64> = Vec::new();
        xtmp.resize(self.len(), 0.0);

        let mut ytmp: Vec<f64> = Vec::new();
        ytmp.resize(self.len(), 0.0);

        let mut ztmp: Vec<f32> = Vec::new();
        ztmp.resize(self.len(), 0.0);

        for (i, j) in indices.into_iter().enumerate() {
            xtmp[j] = self.x[i];
            ytmp[j] = self.y[i];
            ztmp[j] = self.z[i];
        }
        self.x = xtmp;
        self.y = ytmp;
        self.z = ztmp;
    }

    /// Add a new peak-like reference to the feature at a given y "time" coordinate. If the "time"
    /// is not in sorted order, it should automatically re-sort.
    pub fn push<T: CoordinateLike<X> + IntensityMeasurement>(&mut self, pt: &T, time: f64) {
        let x = pt.coordinate();
        let z = pt.intensity();
        self.push_raw(x, time, z);
    }

    /// As [`Feature::push`], but instead add raw values instead of deriving them from
    /// a peak-like reference.
    pub fn push_raw(&mut self, x: f64, y: f64, z: f32) {
        let needs_sort = !self.is_empty() && y < self.y.last().copied().unwrap();
        unsafe { self.push_raw_unchecked(x, y, z) };
        if needs_sort {
            self.sort_by_y();
        }
    }

    /// As [`Feature::push_raw`], but without the automatic sorting.
    ///
    /// # Safety
    /// This method does not enforce the sorting over Y dimension. Use it only if
    /// you do not need to maintain that invariant or intend to sort later.
    pub unsafe fn push_raw_unchecked(&mut self, x: f64, y: f64, z: f32) {
        if !self.is_empty() && y == *self.y.last().unwrap() {
            let last_x = self.x.last().unwrap();
            let last_z = self.z.last().unwrap();
            let new_x = (*last_x * (*last_z as f64) + x * z as f64) / (z + *last_z) as f64;
            *self.x.last_mut().unwrap() = new_x;
            *self.z.last_mut().unwrap() += z;
        } else {
            self.x.push(x);
            self.y.push(y);
            self.z.push(z);
        }
    }

    /// Check if the feature has any points in it
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    pub(crate) fn find_y(&self, y: f64) -> (Option<usize>, f64) {
        if self.is_empty() {
            return (None, y);
        }
        match self.y.binary_search_by(|yi| y.total_cmp(yi).reverse()) {
            Ok(i) => {
                let low = i.saturating_sub(3);
                (low..(low + 6).min(self.len()))
                    .map(|i| (Some(i), (self.y[i] - y).abs()))
                    .min_by(|(_, e), (_, d)| e.total_cmp(d))
                    .unwrap()
            }
            Err(i) => {
                let low = i.saturating_sub(3);
                (low..(low + 6).min(self.len()))
                    .map(|i| (Some(i), (self.y[i] - y).abs()))
                    .min_by(|(_, e), (_, d)| e.total_cmp(d))
                    .unwrap()
            }
        }
    }

    /// Create an iterator that yields (x, y, intensity) references
    pub fn iter(&self) -> Iter<'_, X, Y> {
        Iter::new(self)
    }

    /// Create an iterator that yields (x, y, intensity) mutable references
    pub fn iter_mut(&mut self) -> IterMut<'_, X, Y> {
        IterMut::new(self)
    }

    pub fn as_view(&self) -> FeatureView<'_, X, Y> {
        FeatureView::new(&self.x, &self.y, &self.z)
    }

    pub fn into_inner(self) -> (Vec<f64>, Vec<f64>, Vec<f32>) {
        (self.x, self.y, self.z)
    }

    fn integrate_y(&self) -> f32 {
        self.trapezoid_integrate(&self.y, &self.z)
    }

    /// Sum over the intensity dimension.
    ///
    /// This is not the **area under the curve**
    pub fn total_intensity(&self) -> f32 {
        self.z.iter().sum()
    }

    /// Integrate the feature in the Y dimension
    ///
    /// This uses trapezoid integration, and low quality features
    /// may be truncated.
    pub fn area(&self) -> f32 {
        self.trapezoid_integrate(&self.y, &self.z)
    }
}

impl<X, Y, P: CoordinateLike<X> + IntensityMeasurement> Extend<(P, f64)> for Feature<X, Y> {
    fn extend<T: IntoIterator<Item = (P, f64)>>(&mut self, iter: T) {
        for (x, t) in iter {
            self.push(&x, t)
        }
    }
}

impl<X, Y> Extend<(f64, f64, f32)> for Feature<X, Y> {
    fn extend<T: IntoIterator<Item = (f64, f64, f32)>>(&mut self, iter: T) {
        for (x, y, z) in iter {
            self.push_raw(x, y, z);
        }
    }
}

impl<X, Y, P: CoordinateLike<X> + IntensityMeasurement> FromIterator<(P, f64)> for Feature<X, Y> {
    fn from_iter<T: IntoIterator<Item = (P, f64)>>(iter: T) -> Self {
        let mut this = Self::empty();
        this.extend(iter);
        this
    }
}

impl<X, Y> FromIterator<(f64, f64, f32)> for Feature<X, Y> {
    fn from_iter<T: IntoIterator<Item = (f64, f64, f32)>>(iter: T) -> Self {
        let mut this = Self::empty();
        this.extend(iter);
        this
    }
}

impl<X, Y> PartialEq for Feature<X, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x
            && self.y == other.y
            && self.z == other.z
            && self._x == other._x
            && self._y == other._y
    }
}

impl<X, Y> PartialOrd for Feature<X, Y> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self == other {
            return Some(Ordering::Equal);
        }
        match self.coordinate_x().total_cmp(&other.coordinate_x()) {
            Ordering::Equal => {}
            x => return Some(x),
        };
        self.y.first().partial_cmp(&other.y.first())
    }
}

impl<X, Y> CoordinateLike<X> for Feature<X, Y> {
    fn coordinate(&self) -> f64 {
        self.coordinate_x()
    }
}

impl<X, Y> IntensityMeasurement for Feature<X, Y> {
    fn intensity(&self) -> f32 {
        self.total_intensity()
    }
}

impl<X, Y> TimeInterval<Y> for Feature<X, Y> {
    fn apex_time(&self) -> Option<f64> {
        self.apex_y()
    }

    fn area(&self) -> f32 {
        self.integrate_y()
    }

    fn end_time(&self) -> Option<f64> {
        self.y.last().copied()
    }

    fn start_time(&self) -> Option<f64> {
        self.y.first().copied()
    }

    fn iter_time(&self) -> impl Iterator<Item = f64> {
        self.y.iter().copied()
    }

    fn find_time(&self, time: f64) -> (Option<usize>, f64) {
        self.find_y(time)
    }
}

impl<Y> Feature<MZ, Y> {
    pub fn iter_peaks(&self) -> MZPeakIter<'_, Y> {
        MZPeakIter::new(self)
    }
}

impl<X, Y> TimeArray<Y> for Feature<X, Y> {
    fn time_view(&self) -> &[f64] {
        &self.y
    }
}

pub type LCMSFeature = Feature<MZ, Time>;
pub type IMSFeature = Feature<MZ, IonMobility>;

impl<X, Y> FeatureLike<X, Y> for Feature<X, Y>
where
    Feature<X, Y>: TimeInterval<Y>,
{
    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> impl Iterator<Item = (&f64, &f64, &f32)> {
        self.iter()
    }
}

impl<X, Y> FeatureLikeMut<X, Y> for Feature<X, Y>
where
    Feature<X, Y>: TimeInterval<Y>,
{
    fn iter_mut(&mut self) -> impl Iterator<Item = (&mut f64, &mut f64, &mut f32)> {
        self.iter_mut()
    }

    fn push<T: CoordinateLike<X> + IntensityMeasurement>(&mut self, pt: &T, time: f64) {
        self.push(pt, time)
    }

    fn push_raw(&mut self, x: f64, y: f64, z: f32) {
        self.push_raw(x, y, z)
    }
}

/// An iterator producing immutable references to feature data
/// as owned by [`Feature`].
pub struct Iter<'a, X, Y> {
    xiter: slice::Iter<'a, f64>,
    yiter: slice::Iter<'a, f64>,
    ziter: slice::Iter<'a, f32>,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<'a, X, Y> Iterator for Iter<'a, X, Y> {
    type Item = (&'a f64, &'a f64, &'a f32);

    fn next(&mut self) -> Option<Self::Item> {
        let x = self.xiter.next();
        let y = self.yiter.next();
        let z = self.ziter.next();
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((x, y, z)),
            _ => None,
        }
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        let (x, y, z) = (self.xiter.nth(n), self.yiter.nth(n), self.ziter.nth(n));
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((x, y, z)),
            _ => None,
        }
    }
}

impl<'a, X, Y> ExactSizeIterator for Iter<'a, X, Y> {
    fn len(&self) -> usize {
        self.xiter.len()
    }
}

impl<'a, X, Y> DoubleEndedIterator for Iter<'a, X, Y> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let x = self.xiter.next_back();
        let y = self.yiter.next_back();
        let z = self.ziter.next_back();
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((x, y, z)),
            _ => None,
        }
    }
}

impl<'a, X, Y> Iter<'a, X, Y> {
    pub fn new(source: &'a Feature<X, Y>) -> Self {
        Self {
            xiter: source.x.iter(),
            yiter: source.y.iter(),
            ziter: source.z.iter(),
            _x: PhantomData,
            _y: PhantomData,
        }
    }
}

/// An iterator over [`Feature`] that produces [`CentroidPeak`] instances and
/// an associated time point.
pub struct MZPeakIter<'a, Y> {
    source: Iter<'a, MZ, Y>,
}

impl<'a, Y> MZPeakIter<'a, Y> {
    pub fn new(source: &'a Feature<MZ, Y>) -> Self {
        let iter = source.iter();
        Self { source: iter }
    }
}

impl<'a, Y> Iterator for MZPeakIter<'a, Y> {
    type Item = (CentroidPeak, f64);

    fn next(&mut self) -> Option<Self::Item> {
        let xyz = self.source.next();
        match xyz {
            Some((x, y, z)) => Some((CentroidPeak::new(*x, *z, 0), *y)),
            _ => None,
        }
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        if let Some((x, y, z)) = self.source.nth(n) {
            Some((CentroidPeak::new(*x, *z, 0), *y))
        } else {
            None
        }
    }
}

impl<'a, Y> ExactSizeIterator for MZPeakIter<'a, Y> {
    fn len(&self) -> usize {
        self.source.len()
    }
}

impl<'a, Y> DoubleEndedIterator for MZPeakIter<'a, Y> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let xyz = self.source.next_back();
        match xyz {
            Some((x, y, z)) => Some((CentroidPeak::new(*x, *z, 0), *y)),
            _ => None,
        }
    }
}

/// An iterator producing mutable references to feature data
/// as owned by [`Feature`].
pub struct IterMut<'a, X, Y> {
    xiter: slice::IterMut<'a, f64>,
    yiter: slice::IterMut<'a, f64>,
    ziter: slice::IterMut<'a, f32>,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<'a, X, Y> Iterator for IterMut<'a, X, Y> {
    type Item = (&'a mut f64, &'a mut f64, &'a mut f32);

    fn next(&mut self) -> Option<Self::Item> {
        let x = self.xiter.next();
        let y = self.yiter.next();
        let z = self.ziter.next();
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((x, y, z)),
            _ => None,
        }
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        let x = self.xiter.nth(n);
        let y = self.yiter.nth(n);
        let z = self.ziter.nth(n);
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((x, y, z)),
            _ => None,
        }
    }
}

impl<'a, X, Y> ExactSizeIterator for IterMut<'a, X, Y> {
    fn len(&self) -> usize {
        self.xiter.len()
    }
}

impl<'a, X, Y> IterMut<'a, X, Y> {
    pub fn new(source: &'a mut Feature<X, Y>) -> Self {
        Self {
            xiter: source.x.iter_mut(),
            yiter: source.y.iter_mut(),
            ziter: source.z.iter_mut(),
            _x: PhantomData,
            _y: PhantomData,
        }
    }
}

impl<'a, X, Y> DoubleEndedIterator for IterMut<'a, X, Y> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let x = self.xiter.next_back();
        let y = self.yiter.next_back();
        let z = self.ziter.next_back();
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((x, y, z)),
            _ => None,
        }
    }
}

/// A consuming iterator for [`Feature`]
pub struct IntoIter<X, Y> {
    xiter: std::vec::IntoIter<f64>,
    yiter: std::vec::IntoIter<f64>,
    ziter: std::vec::IntoIter<f32>,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<X, Y> Iterator for IntoIter<X, Y> {
    type Item = (f64, f64, f32);

    fn next(&mut self) -> Option<Self::Item> {
        let x = self.xiter.next();
        let y = self.yiter.next();
        let z = self.ziter.next();
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((x, y, z)),
            _ => None,
        }
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        let x = self.xiter.nth(n);
        let y = self.yiter.nth(n);
        let z = self.ziter.nth(n);
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((x, y, z)),
            _ => None,
        }
    }
}

impl<X, Y> DoubleEndedIterator for IntoIter<X, Y> {
    fn next_back(&mut self) -> Option<Self::Item> {
        let x = self.xiter.next_back();
        let y = self.yiter.next_back();
        let z = self.ziter.next_back();
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((x, y, z)),
            _ => None,
        }
    }
}

impl<X, Y> ExactSizeIterator for IntoIter<X, Y> {
    fn len(&self) -> usize {
        self.xiter.len()
    }
}

impl<X, Y> IntoIter<X, Y> {
    pub fn new(source: Feature<X, Y>) -> Self {
        Self {
            xiter: source.x.into_iter(),
            yiter: source.y.into_iter(),
            ziter: source.z.into_iter(),
            _x: PhantomData,
            _y: PhantomData,
        }
    }
}

impl<X, Y> IntoIterator for Feature<X, Y> {
    type Item = <IntoIter<X, Y> as Iterator>::Item;

    type IntoIter = IntoIter<X, Y>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIter::new(self)
    }
}

impl<'a, X, Y> SplittableFeatureLike<'a, X, Y> for Feature<X, Y> {
    type ViewType = FeatureView<'a, X, Y>;

    fn split_at_time(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        if let Some(point) = self.find_time(point).0 {
            self.split_at(point)
        } else {
            let before = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            let after = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            (before, after)
        }
    }

    fn slice<I: RangeBounds<usize> + Clone>(&'a self, bounds: I) -> Self::ViewType {
        let start = bounds.start_bound();
        let end = bounds.end_bound();
        match (start, end) {
            (Bound::Included(i), Bound::Included(j)) => {
                Self::ViewType::new(&self.x[*i..=*j], &self.y[*i..=*j], &self.z[*i..=*j])
            }
            (Bound::Included(i), Bound::Excluded(j)) => {
                Self::ViewType::new(&self.x[*i..*j], &self.y[*i..*j], &self.z[*i..*j])
            }
            (Bound::Included(i), Bound::Unbounded) => {
                Self::ViewType::new(&self.x[*i..], &self.y[*i..], &self.z[*i..])
            }
            (Bound::Excluded(i), Bound::Included(j)) => {
                Self::ViewType::new(&self.x[*i..*j], &self.y[*i..*j], &self.z[*i..*j])
            }
            (Bound::Excluded(i), Bound::Excluded(j)) => {
                Self::ViewType::new(&self.x[*i..*j], &self.y[*i..*j], &self.z[*i..*j])
            }
            (Bound::Excluded(i), Bound::Unbounded) => {
                Self::ViewType::new(&self.x[*i..], &self.y[*i..], &self.z[*i..])
            }
            (Bound::Unbounded, Bound::Included(j)) => {
                Self::ViewType::new(&self.x[..=*j], &self.y[..=*j], &self.z[..=*j])
            }
            (Bound::Unbounded, Bound::Excluded(j)) => {
                Self::ViewType::new(&self.x[..*j], &self.y[..*j], &self.z[..*j])
            }
            (Bound::Unbounded, Bound::Unbounded) => {
                Self::ViewType::new(&self.x[..], &self.y[..], &self.z[..])
            }
        }
    }

    fn split_at(&'a self, index: usize) -> (Self::ViewType, Self::ViewType) {
        if !self.is_empty() {
            let before = Self::ViewType::new(&self.x[..index], &self.y[..index], &self.z[..index]);
            let after = Self::ViewType::new(&self.x[index..], &self.y[index..], &self.z[index..]);
            (before, after)
        } else {
            let before = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            let after = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            (before, after)
        }
    }
}

/// A non-owning version of [`Feature`]
#[derive(Debug, Clone, Copy)]
pub struct FeatureView<'a, X, Y> {
    x: &'a [f64],
    y: &'a [f64],
    z: &'a [f32],
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<'a, X, Y> CoArrayOps for FeatureView<'a, X, Y> {}

impl<'a, X, Y> FeatureView<'a, X, Y> {
    pub fn new(x: &'a [f64], y: &'a [f64], z: &'a [f32]) -> Self {
        Self {
            x,
            y,
            z,
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    pub fn empty() -> Self {
        Self::new(EMPTY_X, EMPTY_Y, EMPTY_Z)
    }

    fn coordinate_x(&self) -> f64 {
        self.weighted_average(self.x, self.z)
    }

    fn coordinate_y(&self) -> f64 {
        self.weighted_average(self.y, self.z)
    }

    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    fn apex_y(&self) -> Option<f64> {
        self.apex_of(self.y, self.z)
    }

    pub fn to_owned(&self) -> Feature<X, Y> {
        Feature::new(self.x.to_owned(), self.y.to_owned(), self.z.to_owned())
    }

    pub fn into_inner(self) -> (&'a [f64], &'a [f64], &'a [f32]) {
        (self.x, self.y, self.z)
    }

    fn find_y(&self, y: f64) -> (Option<usize>, f64) {
        if self.is_empty() {
            return (None, y);
        }
        match self.y.binary_search_by(|yi| y.total_cmp(yi).reverse()) {
            Ok(i) => {
                let low = i.saturating_sub(3);
                (low..(low + 6).min(self.len()))
                    .map(|i| (Some(i), (self.y[i] - y).abs()))
                    .min_by(|(_, e), (_, d)| e.total_cmp(d))
                    .unwrap()
            }
            Err(i) => {
                let low = i.saturating_sub(3);
                (low..(low + 6).min(self.len()))
                    .map(|i| (Some(i), (self.y[i] - y).abs()))
                    .min_by(|(_, e), (_, d)| e.total_cmp(d))
                    .unwrap()
            }
        }
    }

    pub fn iter(&self) -> Iter<'a, X, Y> {
        Iter {
            xiter: self.x.iter(),
            yiter: self.y.iter(),
            ziter: self.z.iter(),
            _x: PhantomData,
            _y: PhantomData,
        }
    }
}

impl<'a, X, Y> PartialEq for FeatureView<'a, X, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x
            && self.y == other.y
            && self.z == other.z
            && self._x == other._x
            && self._y == other._y
    }
}

impl<'a, X, Y> PartialOrd for FeatureView<'a, X, Y> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self == other {
            return Some(Ordering::Equal);
        }
        match self.coordinate_x().total_cmp(&other.coordinate_x()) {
            Ordering::Equal => {}
            x => return Some(x),
        };
        Some(self.coordinate_y().total_cmp(&other.coordinate_y()))
    }
}

impl<'a, X, Y> CoordinateLike<X> for FeatureView<'a, X, Y> {
    fn coordinate(&self) -> f64 {
        self.coordinate_x()
    }
}

impl<'a, X, Y> TimeInterval<Y> for FeatureView<'a, X, Y> {
    fn apex_time(&self) -> Option<f64> {
        self.apex_y()
    }

    fn area(&self) -> f32 {
        self.trapezoid_integrate(self.y, self.z)
    }

    fn end_time(&self) -> Option<f64> {
        self.y.last().copied()
    }

    fn start_time(&self) -> Option<f64> {
        self.y.first().copied()
    }

    fn iter_time(&self) -> impl Iterator<Item = f64> {
        self.y.iter().copied()
    }

    fn find_time(&self, time: f64) -> (Option<usize>, f64) {
        self.find_y(time)
    }
}

impl<'a, X, Y> IntensityMeasurement for FeatureView<'a, X, Y> {
    fn intensity(&self) -> f32 {
        self.z.iter().sum()
    }
}

impl<'a, X, Y> FeatureLike<X, Y> for FeatureView<'a, X, Y>
where
    FeatureView<'a, X, Y>: TimeInterval<Y>,
{
    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> impl Iterator<Item = (&f64, &f64, &f32)> {
        self.iter()
    }

    fn is_empty(&self) -> bool {
        self.is_empty()
    }
}

impl<'a, X, Y> SplittableFeatureLike<'a, X, Y> for FeatureView<'a, X, Y> {
    type ViewType = Self;

    fn split_at_time(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        if let Some(point) = self.find_time(point).0 {
            self.split_at(point)
        } else {
            let before = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            let after = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            (before, after)
        }
    }

    fn slice<I: RangeBounds<usize> + Clone>(&'a self, bounds: I) -> Self::ViewType {
        let start = bounds.start_bound();
        let end = bounds.end_bound();
        match (start, end) {
            (Bound::Included(i), Bound::Included(j)) => {
                Self::ViewType::new(&self.x[*i..=*j], &self.y[*i..=*j], &self.z[*i..=*j])
            }
            (Bound::Included(i), Bound::Excluded(j)) => {
                Self::ViewType::new(&self.x[*i..=*j], &self.y[*i..=*j], &self.z[*i..=*j])
            }
            (Bound::Included(i), Bound::Unbounded) => {
                Self::ViewType::new(&self.x[*i..], &self.y[*i..], &self.z[*i..])
            }
            (Bound::Excluded(i), Bound::Included(j)) => {
                Self::ViewType::new(&self.x[*i..*j], &self.y[*i..*j], &self.z[*i..*j])
            }
            (Bound::Excluded(i), Bound::Excluded(j)) => {
                Self::ViewType::new(&self.x[*i..=*j], &self.y[*i..=*j], &self.z[*i..=*j])
            }
            (Bound::Excluded(i), Bound::Unbounded) => {
                Self::ViewType::new(&self.x[*i..], &self.y[*i..], &self.z[*i..])
            }
            (Bound::Unbounded, Bound::Included(j)) => {
                Self::ViewType::new(&self.x[..=*j], &self.y[..=*j], &self.z[..=*j])
            }
            (Bound::Unbounded, Bound::Excluded(j)) => {
                Self::ViewType::new(&self.x[..*j], &self.y[..*j], &self.z[..*j])
            }
            (Bound::Unbounded, Bound::Unbounded) => {
                Self::ViewType::new(self.x, self.y, self.z)
            }
        }
    }

    fn split_at(&'a self, index: usize) -> (Self::ViewType, Self::ViewType) {
        if !self.is_empty() {
            let before = Self::ViewType::new(&self.x[..index], &self.y[..index], &self.z[..index]);
            let after = Self::ViewType::new(&self.x[index..], &self.y[index..], &self.z[index..]);
            (before, after)
        } else {
            let before = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            let after = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            (before, after)
        }
    }
}


impl<'a, X, Y> TimeArray<Y> for FeatureView<'a, X, Y> {
    fn time_view(&self) -> &[f64] {
        self.y
    }
}