//! A feature is a two dimensional mass spectrum concept for a measure over some time unit.
//! It represents something that is located at a constrained but varying coordinate system `X`
//! over a sequentially ordered dimension `Y` with an abundance measure at each time point.
//!

use core::slice;
use std::{cmp::Ordering, marker::PhantomData};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{
    coordinate::{CoordinateLike, IonMobility, Mass, Time, MZ},
    CentroidPeak, CoordinateRange, DeconvolutedPeak, IntensityMeasurement, KnownCharge,
    MassLocated,
};

#[derive(PartialEq)]
struct NonNan(f64);

impl NonNan {
    fn new(val: f64) -> Option<NonNan> {
        if val.is_nan() {
            None
        } else {
            Some(NonNan(val))
        }
    }
}

impl Eq for NonNan {}

impl PartialOrd for NonNan {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for NonNan {
    fn cmp(&self, other: &NonNan) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

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
}

trait CoArrayOps {
    fn weighted_average(&self, x: &[f64], w: &[f32]) -> f64 {
        let (acc, norm) = x
            .iter()
            .zip(w.iter())
            .fold((0.0, 0.0), |(acc, norm), (x, z)| {
                let norm = norm + *z as f64;
                let acc = acc + *x * *z as f64;
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

    /// Compute a weighted average over the X dimension
    fn coordinate_x(&self) -> f64 {
        self.weighted_average(&self.x, &self.z)
    }

    #[allow(unused)]
    /// Compute a weighted average over the Y dimension
    fn coordinate_y(&self) -> f64 {
        self.weighted_average(&self.y, &self.z)
    }

    /// The number of points in the feature
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Find the time where the feature achieves its maximum abundance
    fn apex_y(&self) -> Option<f64> {
        self.apex_of(&self.y, &self.z)
    }

    /// Sort the feature by the Y dimension
    fn sort_by_y(&mut self) {
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
            let new_x = *last_x * (*last_z as f64) + (x * z as f64) / (z + *last_z) as f64;
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

    fn find_y(&self, y: f64) -> (Option<usize>, f64) {
        if self.is_empty() {
            return (None, y);
        }
        match self.y.binary_search_by(|yi| y.total_cmp(yi)) {
            Ok(i) => {
                let low = i.saturating_sub(1);
                (low..(low + 3).min(self.len()))
                    .map(|i| (Some(i), (self.y[i] - y).abs()))
                    .min_by(|(_, e), (_, d)| e.total_cmp(d))
                    .unwrap()
            }
            Err(i) => {
                let low = i.saturating_sub(1);
                (low..(low + 3).min(self.len()))
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

pub struct MZPeakIter<'a, Y> {
    source: &'a Feature<MZ, Y>,
    xiter: slice::Iter<'a, f64>,
    yiter: slice::Iter<'a, f64>,
    ziter: slice::Iter<'a, f32>,
}

impl<'a, Y> MZPeakIter<'a, Y> {
    pub fn new(source: &'a Feature<MZ, Y>) -> Self {
        Self {
            source,
            xiter: source.x.iter(),
            yiter: source.y.iter(),
            ziter: source.z.iter(),
        }
    }
}

impl<'a, Y> Iterator for MZPeakIter<'a, Y> {
    type Item = (CentroidPeak, f64);

    fn next(&mut self) -> Option<Self::Item> {
        let x = self.xiter.next();
        let y = self.yiter.next();
        let z = self.ziter.next();
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((CentroidPeak::new(*x, *z, 0), *y)),
            _ => None,
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
        let x = self.xiter.next_back();
        let y = self.yiter.next_back();
        let z = self.ziter.next_back();
        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some((CentroidPeak::new(*x, *z, 0), *y)),
            _ => None,
        }
    }
}

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

#[derive(Debug, Default, Clone)]
pub struct ChargedFeature<X, Y> {
    pub feature: Feature<X, Y>,
    pub charge: i32,
}

impl<X, Y> CoordinateLike<X> for ChargedFeature<X, Y> {
    fn coordinate(&self) -> f64 {
        self.feature.coordinate_x()
    }
}

impl<X, Y> FeatureLikeMut<X, Y> for ChargedFeature<X, Y>
where
    Feature<X, Y>: FeatureLikeMut<X, Y>,
{
    fn iter_mut(&mut self) -> impl Iterator<Item = (&mut f64, &mut f64, &mut f32)> {
        <Feature<X, Y> as FeatureLikeMut<X, Y>>::iter_mut(&mut self.feature)
    }

    fn push<T: CoordinateLike<X> + IntensityMeasurement>(&mut self, pt: &T, time: f64) {
        <Feature<X, Y> as FeatureLikeMut<X, Y>>::push(&mut self.feature, pt, time)
    }

    fn push_raw(&mut self, x: f64, y: f64, z: f32) {
        <Feature<X, Y> as FeatureLikeMut<X, Y>>::push_raw(&mut self.feature, x, y, z)
    }
}

impl<X, Y> TimeInterval<Y> for ChargedFeature<X, Y>
where
    Feature<X, Y>: TimeInterval<Y>,
{
    fn apex_time(&self) -> Option<f64> {
        <Feature<X, Y> as TimeInterval<Y>>::apex_time(&self.feature)
    }

    fn area(&self) -> f32 {
        <Feature<X, Y> as TimeInterval<Y>>::area(&self.feature)
    }

    fn end_time(&self) -> Option<f64> {
        <Feature<X, Y> as TimeInterval<Y>>::end_time(&self.feature)
    }

    fn start_time(&self) -> Option<f64> {
        <Feature<X, Y> as TimeInterval<Y>>::start_time(&self.feature)
    }

    fn iter_time(&self) -> impl Iterator<Item = f64> {
        self.feature.iter_time()
    }

    fn find_time(&self, time: f64) -> (Option<usize>, f64) {
        self.feature.find_y(time)
    }
}

impl<X, Y> FeatureLike<X, Y> for ChargedFeature<X, Y>
where
    Feature<X, Y>: FeatureLike<X, Y>,
{
    fn len(&self) -> usize {
        <Feature<X, Y> as FeatureLike<X, Y>>::len(&self.feature)
    }

    fn iter(&self) -> impl Iterator<Item = (&f64, &f64, &f32)> {
        <Feature<X, Y> as FeatureLike<X, Y>>::iter(&self.feature)
    }
}

impl<X, Y> Extend<(f64, f64, f32)> for ChargedFeature<X, Y> {
    fn extend<T: IntoIterator<Item = (f64, f64, f32)>>(&mut self, iter: T) {
        <Feature<X, Y> as Extend<(f64, f64, f32)>>::extend(&mut self.feature, iter)
    }
}

impl<P: CoordinateLike<X> + IntensityMeasurement, X, Y> Extend<(P, f64)> for ChargedFeature<X, Y> {
    fn extend<T: IntoIterator<Item = (P, f64)>>(&mut self, iter: T) {
        <Feature<X, Y> as Extend<(P, f64)>>::extend(&mut self.feature, iter)
    }
}

impl<X, Y> IntensityMeasurement for ChargedFeature<X, Y> {
    fn intensity(&self) -> f32 {
        <Feature<X, Y> as IntensityMeasurement>::intensity(&self.feature)
    }
}

impl<X, Y> ChargedFeature<X, Y> {
    pub fn new(feature: Feature<X, Y>, charge: i32) -> Self {
        Self { feature, charge }
    }

    pub fn empty(charge: i32) -> Self {
        Self {
            feature: Feature::empty(),
            charge,
        }
    }

    pub fn iter(&self) -> Iter<'_, X, Y> {
        self.feature.iter()
    }

    pub fn iter_mut(&mut self) -> IterMut<'_, X, Y> {
        self.feature.iter_mut()
    }

    pub fn push<T: CoordinateLike<X> + IntensityMeasurement>(&mut self, pt: &T, time: f64) {
        self.feature.push(pt, time)
    }

    pub fn push_raw(&mut self, x: f64, y: f64, z: f32) {
        self.feature.push_raw(x, y, z)
    }

    /// # Safety
    /// This method does not enforce the sorting over Y dimension. Use it only if
    /// you do not need to maintain that invariant or intend to sort later.
    pub unsafe fn push_raw_unchecked(&mut self, x: f64, y: f64, z: f32) {
        self.feature.push_raw_unchecked(x, y, z)
    }

    pub fn len(&self) -> usize {
        self.feature.len()
    }

    pub fn is_empty(&self) -> bool {
        self.feature.is_empty()
    }
}

impl<Y> ChargedFeature<Mass, Y> {
    pub fn iter_peaks(&self) -> DeconvolutedPeakIter<'_, Y> {
        DeconvolutedPeakIter::new(self)
    }
}

impl<X, Y> PartialEq for ChargedFeature<X, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.feature == other.feature && self.charge == other.charge
    }
}

impl<X, Y> PartialOrd for ChargedFeature<X, Y> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match self.feature.partial_cmp(&other.feature) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        self.charge.partial_cmp(&other.charge)
    }
}

impl<Y> CoordinateLike<MZ> for ChargedFeature<Mass, Y> {
    fn coordinate(&self) -> f64 {
        let charge_carrier: f64 = 1.007276;
        let charge = self.charge as f64;
        (self.neutral_mass() + charge_carrier * charge) / charge
    }
}

impl<X, Y> KnownCharge for ChargedFeature<X, Y> {
    fn charge(&self) -> i32 {
        self.charge
    }
}

impl<X, Y> AsRef<Feature<X, Y>> for ChargedFeature<X, Y> {
    fn as_ref(&self) -> &Feature<X, Y> {
        &self.feature
    }
}

impl<X, Y> AsMut<Feature<X, Y>> for ChargedFeature<X, Y> {
    fn as_mut(&mut self) -> &mut Feature<X, Y> {
        &mut self.feature
    }
}

impl<X, Y, P: CoordinateLike<X> + IntensityMeasurement + KnownCharge> FromIterator<(P, f64)>
    for ChargedFeature<X, Y>
{
    fn from_iter<T: IntoIterator<Item = (P, f64)>>(iter: T) -> Self {
        let mut it = iter.into_iter();
        if let Some((peak, time)) = it.next() {
            let mut this = Self::empty(peak.charge());
            this.push(&peak, time);
            this.extend(it);
            this
        } else {
            Self::empty(0)
        }
    }
}

pub type DeconvolvedLCMSFeature = ChargedFeature<Mass, Time>;
pub type DeconvolvedIMSFeature = ChargedFeature<Mass, IonMobility>;

pub struct DeconvolutedPeakIter<'a, Y> {
    source: &'a ChargedFeature<Mass, Y>,
    point_iter: Iter<'a, Mass, Y>,
}

impl<'a, Y> DeconvolutedPeakIter<'a, Y> {
    pub fn new(source: &'a ChargedFeature<Mass, Y>) -> Self {
        Self {
            source,
            point_iter: source.iter(),
        }
    }
}

impl<'a, Y> Iterator for DeconvolutedPeakIter<'a, Y> {
    type Item = (DeconvolutedPeak, f64);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((mass, time, intensity)) = self.point_iter.next() {
            Some((
                DeconvolutedPeak::new(*mass, *intensity, self.source.charge, 0),
                *time,
            ))
        } else {
            None
        }
    }
}

impl<'a, Y> ExactSizeIterator for DeconvolutedPeakIter<'a, Y> {
    fn len(&self) -> usize {
        self.source.len()
    }
}

impl<'a, Y> DoubleEndedIterator for DeconvolutedPeakIter<'a, Y> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if let Some((mass, time, intensity)) = self.point_iter.next_back() {
            Some((
                DeconvolutedPeak::new(*mass, *intensity, self.source.charge, 0),
                *time,
            ))
        } else {
            None
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct SimpleFeature<X, Y> {
    pub label: f64,
    y: Vec<f64>,
    z: Vec<f32>,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<X, Y> CoArrayOps for SimpleFeature<X, Y> {}

impl<X, Y> SimpleFeature<X, Y> {
    pub fn empty(label: f64) -> Self {
        Self {
            label,
            y: Vec::new(),
            z: Vec::new(),
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    fn sort_by_y(&mut self) {
        let mut indices: Vec<_> = (0..self.len()).collect();
        indices.sort_by_key(|i| NonNan::new(self.y[*i]));

        let mut ytmp: Vec<f64> = Vec::new();
        ytmp.resize(self.len(), 0.0);

        let mut ztmp: Vec<f32> = Vec::new();
        ztmp.resize(self.len(), 0.0);

        for (i, j) in indices.into_iter().enumerate() {
            ytmp[j] = self.y[i];
            ztmp[j] = self.z[i];
        }
        self.y = ytmp;
        self.z = ztmp;
    }

    pub fn push<T: CoordinateLike<X> + IntensityMeasurement>(&mut self, pt: &T, time: f64) {
        self.push_raw(0.0, time, pt.intensity())
    }

    pub fn push_raw(&mut self, _x: f64, y: f64, z: f32) {
        let needs_sort = self.len() > 0 && self.y.last().copied().unwrap() > y;
        self.y.push(y);
        self.z.push(z);
        if needs_sort {
            self.sort_by_y();
        }
    }

    fn apex_y(&self) -> Option<f64> {
        self.apex_of(&self.y, &self.z)
    }

    fn find_y(&self, y: f64) -> (Option<usize>, f64) {
        if self.is_empty() {
            return (None, y);
        }
        match self.y.binary_search_by(|yi| y.total_cmp(yi)) {
            Ok(i) => {
                let low = i.saturating_sub(1);
                (low..(low + 3).min(self.len()))
                    .map(|i| (Some(i), (self.y[i] - y).abs()))
                    .min_by(|(_, e), (_, d)| e.total_cmp(d))
                    .unwrap()
            }
            Err(i) => {
                let low = i.saturating_sub(1);
                (low..(low + 3).min(self.len()))
                    .map(|i| (Some(i), (self.y[i] - y).abs()))
                    .min_by(|(_, e), (_, d)| e.total_cmp(d))
                    .unwrap()
            }
        }
    }
}

impl<X, Y> PartialEq for SimpleFeature<X, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.label == other.label
            && self.y == other.y
            && self.z == other.z
            && self._x == other._x
            && self._y == other._y
    }
}

impl<X, Y> PartialOrd for SimpleFeature<X, Y> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match self.label.partial_cmp(&other.label) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self.y.partial_cmp(&other.y) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        self.z.partial_cmp(&other.z)
    }
}

impl<X, Y> CoordinateLike<X> for SimpleFeature<X, Y> {
    fn coordinate(&self) -> f64 {
        self.label
    }
}

impl<X, Y> TimeInterval<Y> for SimpleFeature<X, Y> {
    fn start_time(&self) -> Option<f64> {
        self.y.first().copied()
    }

    fn end_time(&self) -> Option<f64> {
        self.y.last().copied()
    }

    fn apex_time(&self) -> Option<f64> {
        self.apex_y()
    }

    fn area(&self) -> f32 {
        let mut it = self.y.iter().zip(self.z.iter());
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

    fn iter_time(&self) -> impl Iterator<Item = f64> {
        self.y.iter().copied()
    }

    fn find_time(&self, time: f64) -> (Option<usize>, f64) {
        self.find_y(time)
    }
}

impl<X, Y> IntensityMeasurement for SimpleFeature<X, Y> {
    fn intensity(&self) -> f32 {
        self.z.iter().sum()
    }
}

impl<X, Y> FeatureLike<X, Y> for SimpleFeature<X, Y> {
    fn len(&self) -> usize {
        self.y.len()
    }

    fn iter(&self) -> impl Iterator<Item = (&f64, &f64, &f32)> {
        self.y
            .iter()
            .zip(self.z.iter())
            .map(|(y, z)| (&self.label, y, z))
    }
}

impl<X, Y> Extend<(f64, f64, f32)> for SimpleFeature<X, Y> {
    fn extend<T: IntoIterator<Item = (f64, f64, f32)>>(&mut self, iter: T) {
        for (_x, y, z) in iter {
            self.y.push(y);
            self.z.push(z);
        }
    }
}

impl<X, Y, P: CoordinateLike<X> + IntensityMeasurement> Extend<(P, f64)> for SimpleFeature<X, Y> {
    fn extend<T: IntoIterator<Item = (P, f64)>>(&mut self, iter: T) {
        for (x, t) in iter {
            self.push(&x, t)
        }
    }
}

impl<X, Y, P: CoordinateLike<X> + IntensityMeasurement + KnownCharge> FromIterator<(P, f64)>
    for SimpleFeature<X, Y>
{
    fn from_iter<T: IntoIterator<Item = (P, f64)>>(iter: T) -> Self {
        let mut it = iter.into_iter();
        if let Some((peak, time)) = it.next() {
            let mut this = Self::empty(peak.coordinate());
            this.push(&peak, time);
            this.extend(it);
            this
        } else {
            Self::empty(0.0)
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
        Feature::new(
            self.x.to_owned(),
            self.y.to_owned(),
            self.z.to_owned()
        )
    }

    fn find_y(&self, y: f64) -> (Option<usize>, f64) {
        if self.is_empty() {
            return (None, y);
        }
        match self.y.binary_search_by(|yi| y.total_cmp(yi)) {
            Ok(i) => {
                let low = i.saturating_sub(1);
                (low..(low + 3).min(self.len()))
                    .map(|i| (Some(i), (self.y[i] - y).abs()))
                    .min_by(|(_, e), (_, d)| e.total_cmp(d))
                    .unwrap()
            }
            Err(i) => {
                let low = i.saturating_sub(1);
                (low..(low + 3).min(self.len()))
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

/// A non-owning version of [`ChargedFeature`]
#[derive(Debug, Clone, Copy)]
pub struct ChargedFeatureView<'a, X, Y> {
    feature: FeatureView<'a, X, Y>,
    pub charge: i32,
}

impl<'a, X, Y> TimeInterval<Y> for ChargedFeatureView<'a, X, Y> {
    fn apex_time(&self) -> Option<f64> {
        <FeatureView<'a, X, Y> as TimeInterval<Y>>::apex_time(&self.feature)
    }

    fn area(&self) -> f32 {
        <FeatureView<'a, X, Y> as TimeInterval<Y>>::area(&self.feature)
    }

    fn end_time(&self) -> Option<f64> {
        <FeatureView<'a, X, Y> as TimeInterval<Y>>::end_time(&self.feature)
    }

    fn start_time(&self) -> Option<f64> {
        <FeatureView<'a, X, Y> as TimeInterval<Y>>::start_time(&self.feature)
    }

    fn iter_time(&self) -> impl Iterator<Item = f64> {
        <FeatureView<'a, X, Y> as TimeInterval<Y>>::iter_time(&self.feature)
    }

    fn find_time(&self, time: f64) -> (Option<usize>, f64) {
        <FeatureView<'a, X, Y> as TimeInterval<Y>>::find_time(&self.feature, time)
    }
}

impl<'a, X, Y> PartialEq for ChargedFeatureView<'a, X, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.feature == other.feature && self.charge == other.charge
    }
}

impl<'a, X, Y> PartialOrd for ChargedFeatureView<'a, X, Y> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match self.feature.partial_cmp(&other.feature) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        self.charge.partial_cmp(&other.charge)
    }
}

impl<'a, X, Y> CoordinateLike<X> for ChargedFeatureView<'a, X, Y> {
    fn coordinate(&self) -> f64 {
        self.feature.coordinate()
    }
}

impl<'a, X, Y> IntensityMeasurement for ChargedFeatureView<'a, X, Y> {
    fn intensity(&self) -> f32 {
        self.feature.intensity()
    }
}

impl<'a, X, Y> FeatureLike<X, Y> for ChargedFeatureView<'a, X, Y>
where
    ChargedFeatureView<'a, X, Y>: TimeInterval<Y>,
{
    fn len(&self) -> usize {
        self.feature.len()
    }

    fn iter(&self) -> impl Iterator<Item = (&f64, &f64, &f32)> {
        self.iter()
    }

    fn is_empty(&self) -> bool {
        self.is_empty()
    }
}

impl<'a, X, Y> ChargedFeatureView<'a, X, Y> {
    pub fn new(feature: FeatureView<'a, X, Y>, charge: i32) -> Self {
        Self { feature, charge }
    }

    pub fn len(&self) -> usize {
        self.feature.len()
    }

    pub fn is_empty(&self) -> bool {
        self.feature.is_empty()
    }

    pub fn iter(&self) -> Iter<'a, X, Y> {
        self.feature.iter()
    }

    pub fn to_owned(&self) -> ChargedFeature<X, Y> {
        ChargedFeature::new(
            self.feature.to_owned(), self.charge
        )
    }
}


/// A trait to split features at a given time point
pub trait SplittableFeatureLike<'a, X, Y>: FeatureLike<X, Y> {
    type ViewType: FeatureLike<X, Y>;

    fn split_at(&'a self, point: f64) -> (Self::ViewType, Self::ViewType);
}

const EMPTY_X: &[f64] = &[];
const EMPTY_Y: &[f64] = &[];
const EMPTY_Z: &[f32] = &[];

impl<'a, X, Y> SplittableFeatureLike<'a, X, Y> for Feature<X, Y> {
    type ViewType = FeatureView<'a, X, Y>;

    fn split_at(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        if let Some(point) = self.find_time(point).0 {
            let before = Self::ViewType::new(&self.x[..point], &self.y[..point], &self.z[..point]);
            let after = Self::ViewType::new(&self.x[point..], &self.y[point..], &self.z[point..]);
            (before, after)
        } else {
            let before = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            let after = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            (before, after)
        }
    }
}

impl<'a, X, Y> SplittableFeatureLike<'a, X, Y> for FeatureView<'a, X, Y> {
    type ViewType = Self;

    fn split_at(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        if let Some(point) = self.find_time(point).0 {
            let before = Self::ViewType::new(&self.x[..point], &self.y[..point], &self.z[..point]);
            let after = Self::ViewType::new(&self.x[point..], &self.y[point..], &self.z[point..]);
            (before, after)
        } else {
            let before = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            let after = Self::ViewType::new(EMPTY_X, EMPTY_Y, EMPTY_Z);
            (before, after)
        }
    }
}

impl<'a, X, Y> SplittableFeatureLike<'a, X, Y> for ChargedFeature<X, Y> {
    type ViewType = ChargedFeatureView<'a, X, Y>;

    fn split_at(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        let (before, after) = self.feature.split_at(point);
        (
            Self::ViewType::new(before, self.charge),
            Self::ViewType::new(after, self.charge),
        )
    }
}

impl<'a, X, Y> SplittableFeatureLike<'a, X, Y> for ChargedFeatureView<'a, X, Y> {
    type ViewType = Self;

    fn split_at(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        let (before, after) = self.feature.split_at(point);
        (
            Self::ViewType::new(before, self.charge),
            Self::ViewType::new(after, self.charge),
        )
    }
}



#[cfg(test)]
mod test {
    use super::*;
    use crate::{CentroidPeak, DeconvolutedPeak, MZLocated};

    #[test]
    fn test_build_raw() {
        let mut x = LCMSFeature::empty();

        let points = vec![
            (CentroidPeak::new(204.08, 3432.1, 0), 0.1),
            (CentroidPeak::new(204.07, 7251.9, 0), 0.2),
            (CentroidPeak::new(204.08, 5261.7, 0), 0.3),
        ];

        x.extend(points.iter().cloned());
        assert_eq!(x.len(), 3);

        let y: f32 = points.iter().map(|p| p.0.intensity()).sum();
        assert!((x.intensity() - y).abs() < 1e-6);

        let mz = 204.07545212;
        assert!((x.mz() - mz).abs() < 1e-6);

        let area = 1159.879980;
        assert!((x.area() - area).abs() < 1e-6);

        assert_eq!(x.iter().len(), 3);

        if let Some((pt, t)) = x.iter_peaks().next() {
            let (rpt, rt) = &points[0];
            assert_eq!(pt, rpt);
            assert_eq!(t, *rt);
        }

        assert_eq!(x, points.into_iter().collect());

        let (i, e) = x.find_time(0.3);
        assert_eq!(i, Some(2));
        assert_eq!(e, 0.0);

        let (i, e) = x.find_time(0.5);
        assert_eq!(i, Some(2));
        assert_eq!(e, 0.2);

        let (b, a) = x.split_at(0.2);
        assert_eq!(b.len(), 1);
        assert_eq!(a.len(), 2);

    }

    #[test]
    fn test_build_charged() {
        let mut x = DeconvolvedLCMSFeature::empty(1);

        let points = vec![
            (DeconvolutedPeak::new(203.08, 3432.1, 1, 0), 0.1),
            (DeconvolutedPeak::new(203.07, 7251.9, 1, 0), 0.2),
            (DeconvolutedPeak::new(203.08, 5261.7, 1, 0), 0.3),
        ];

        x.extend(points.iter().cloned());

        assert_eq!(x.len(), 3);

        let y: f32 = points.iter().map(|p| p.0.intensity()).sum();
        assert!((x.intensity() - y).abs() < 1e-6);

        let mass = 203.07545212;
        assert!((x.neutral_mass() - mass).abs() < 1e-6);

        let area = 1159.879980;
        assert!((x.area() - area).abs() < 1e-6);

        assert_eq!(x.iter().len(), 3);

        if let Some((pt, t)) = x.iter_peaks().next() {
            let (rpt, rt) = &points[0];
            assert_eq!(pt, rpt);
            assert_eq!(t, *rt);
        }

        assert_eq!(x, points.into_iter().collect());
    }
}
