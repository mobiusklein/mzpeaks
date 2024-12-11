use std::{
    cmp::Ordering,
    marker::PhantomData,
    ops::{Bound, RangeBounds},
};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{coordinate::CoordinateLike, IntensityMeasurement};

use super::util::{NonNan, EMPTY_Y, EMPTY_Z};
use super::{
    traits::{CoArrayOps, FeatureLike, SplittableFeatureLike, TimeInterval},
    FeatureLikeMut, TimeArray,
};

/// A feature-like type that doesn't have a variable first dimension, instead
/// using a constant value
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SimpleFeature<X, Y> {
    /// A constant value for the first dimension
    pub label: f64,
    y: Vec<f64>,
    z: Vec<f32>,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<X, Y> CoArrayOps for SimpleFeature<X, Y> {}

impl<X, Y> SimpleFeature<X, Y> {
    pub fn new(label: f64, y: Vec<f64>, z: Vec<f32>) -> Self {
        Self {
            label,
            y,
            z,
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    pub fn as_view(&self) -> SimpleFeatureView<'_, X, Y> {
        SimpleFeatureView::new(self.label, &self.y, &self.z)
    }

    pub fn into_inner(self) -> (f64, Vec<f64>, Vec<f32>) {
        (self.label, self.y, self.z)
    }

    /// Create an empty [`SimpleFeature`]
    pub fn empty(label: f64) -> Self {
        Self::with_capacity(0, label)
    }

    /// Createa a new empty [`SimpleFeature`] with pre-allocated capacity
    pub fn with_capacity(capacity: usize, label: f64) -> Self {
        Self {
            label,
            y: Vec::with_capacity(capacity),
            z: Vec::with_capacity(capacity),
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    fn sort_by_y(&mut self) {
        let n = self.len();
        let mut mask: Vec<_> = (0..n).collect();
        mask.sort_by_key(|i| NonNan::new(self.y[*i]));

        const TOMBSTONE: usize = usize::MAX;

        for idx in 0..n {
            if mask[idx] != TOMBSTONE {
                let mut current_idx = idx;
                loop {
                    let next_idx = mask[current_idx];
                    mask[current_idx] = TOMBSTONE;
                    if mask[next_idx] == TOMBSTONE {
                        break;
                    }
                    self.y.swap(current_idx, next_idx);
                    self.z.swap(current_idx, next_idx);
                    current_idx = next_idx;
                }
            }
        }
    }

    /// Push a new data point onto the feature and ensure the time ordering invariant is satisfied.
    pub fn push<T: CoordinateLike<X> + IntensityMeasurement>(&mut self, pt: &T, time: f64) {
        self.push_raw(pt.coordinate(), time, pt.intensity())
    }

    /// Push a new data point onto the feature and ensure the time ordering invariant is satisfied.
    pub fn push_raw(&mut self, _x: f64, y: f64, z: f32) {
        let needs_sort = self.len() > 0 && self.y.last().copied().unwrap() > y;
        if !self.is_empty() && y == *self.y.last().unwrap() {
            *self.z.last_mut().unwrap() += z;
        } else {
            self.y.push(y);
            self.z.push(z);
        }
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
}

impl<'a, X, Y> TimeArray<Y> for SimpleFeature<X, Y> {
    fn time_view(&self) -> &[f64] {
        &self.y
    }

    fn intensity_view(&self) -> &[f32] {
        &self.z
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

    fn iter(&self) -> impl Iterator<Item = (f64, f64, f32)> {
        self.y
            .iter()
            .zip(self.z.iter())
            .map(|(y, z)| (self.label, *y, *z))
    }
}

pub struct IntoIter<X, Y> {
    xval: f64,
    yiter: std::vec::IntoIter<f64>,
    ziter: std::vec::IntoIter<f32>,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<X, Y> Iterator for IntoIter<X, Y> {
    type Item = (f64, f64, f32);

    fn next(&mut self) -> Option<Self::Item> {
        let y = self.yiter.next();
        let z = self.ziter.next();
        match (y, z) {
            (Some(y), Some(z)) => Some((self.xval, y, z)),
            _ => None,
        }
    }
}

impl<X, Y> IntoIterator for SimpleFeature<X, Y> {
    type Item = (f64, f64, f32);

    type IntoIter = IntoIter<X, Y>;

    fn into_iter(self) -> Self::IntoIter {
        let xval = self.coordinate();
        IntoIter {
            xval,
            yiter: self.y.into_iter(),
            ziter: self.z.into_iter(),
            _x: PhantomData::<X>,
            _y: PhantomData::<Y>,
        }
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

impl<X, Y, P: CoordinateLike<X> + IntensityMeasurement> FromIterator<(P, f64)>
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

impl<X, Y> FromIterator<(f64, f64, f32)> for SimpleFeature<X, Y> {
    fn from_iter<T: IntoIterator<Item = (f64, f64, f32)>>(iter: T) -> Self {
        let mut it = iter.into_iter();
        if let Some((x, time, intensity)) = it.next() {
            let mut this = Self::empty(x);
            this.push_raw(x, time, intensity);
            this.extend(it);
            this
        } else {
            Self::empty(0.0)
        }
    }
}

#[derive(Debug)]
pub struct IterMut<'a, X, Y> {
    yiter: std::slice::IterMut<'a, f64>,
    ziter: std::slice::IterMut<'a, f32>,
    value: f64,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<'a, X, Y> Iterator for IterMut<'a, X, Y> {
    type Item = (&'a mut f64, &'a mut f64, &'a mut f32);

    fn next(&mut self) -> Option<Self::Item> {
        let x = std::ptr::addr_of_mut!(self.value);
        let y = self.yiter.next();
        let z = self.ziter.next();

        match (y, z) {
            (Some(y), Some(z)) => Some((unsafe { &mut *x }, y, z)),
            _ => None,
        }
    }
}

impl<X, Y> FeatureLikeMut<X, Y> for SimpleFeature<X, Y> {
    fn iter_mut(&mut self) -> impl Iterator<Item = (&mut f64, &mut f64, &mut f32)> {
        IterMut {
            yiter: self.y.iter_mut(),
            ziter: self.z.iter_mut(),
            value: self.label,
            _x: PhantomData::<X>,
            _y: PhantomData::<Y>,
        }
    }

    fn push<T: CoordinateLike<X> + IntensityMeasurement>(&mut self, pt: &T, time: f64) {
        self.push(pt, time);
    }

    fn push_raw(&mut self, x: f64, y: f64, z: f32) {
        self.push_raw(x, y, z);
    }
}

/// A non-owning version of [`SimpleFeature`]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct SimpleFeatureView<'a, X, Y> {
    pub label: f64,
    y: &'a [f64],
    z: &'a [f32],
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}
impl<'a, X, Y> CoArrayOps for SimpleFeatureView<'a, X, Y> {}

impl<'a, X, Y> SimpleFeatureView<'a, X, Y> {
    pub fn new(label: f64, y: &'a [f64], z: &'a [f32]) -> Self {
        Self {
            label,
            y,
            z,
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    pub fn empty(label: f64) -> Self {
        Self {
            label,
            y: EMPTY_Y,
            z: EMPTY_Z,
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    fn apex_y(&self) -> Option<f64> {
        self.apex_of(self.y, self.z)
    }

    pub fn to_owned(&self) -> SimpleFeature<X, Y> {
        SimpleFeature::new(self.label, self.y.to_vec(), self.z.to_vec())
    }

    pub fn into_inner(self) -> (f64, &'a [f64], &'a [f32]) {
        (self.label, self.y, self.z)
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
}

impl<'a, X, Y> PartialEq for SimpleFeatureView<'a, X, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.label == other.label
            && self.y == other.y
            && self.z == other.z
            && self._x == other._x
            && self._y == other._y
    }
}

impl<'a, X, Y> PartialOrd for SimpleFeatureView<'a, X, Y> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match self.label.partial_cmp(&other.label) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self.y.partial_cmp(other.y) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        self.z.partial_cmp(other.z)
    }
}

impl<'a, X, Y> CoordinateLike<X> for SimpleFeatureView<'a, X, Y> {
    fn coordinate(&self) -> f64 {
        self.label
    }
}

impl<'a, X, Y> TimeInterval<Y> for SimpleFeatureView<'a, X, Y> {
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

impl<'a, X, Y> IntensityMeasurement for SimpleFeatureView<'a, X, Y> {
    fn intensity(&self) -> f32 {
        self.z.iter().sum()
    }
}

impl<'a, X, Y> FeatureLike<X, Y> for SimpleFeatureView<'a, X, Y> {
    fn len(&self) -> usize {
        self.y.len()
    }

    fn iter(&self) -> impl Iterator<Item = (f64, f64, f32)> {
        self.y
            .iter()
            .zip(self.z.iter())
            .map(|(y, z)| (self.label, *y, *z))
    }
}

impl<'a, X, Y> SplittableFeatureLike<'a, X, Y> for SimpleFeature<X, Y> {
    type ViewType = SimpleFeatureView<'a, X, Y>;

    fn split_at(&'a self, index: usize) -> (Self::ViewType, Self::ViewType) {
        if !self.is_empty() {
            let before = Self::ViewType::new(self.label, &self.y[..index], &self.z[..index]);
            let after = Self::ViewType::new(self.label, &self.y[index..], &self.z[index..]);
            (before, after)
        } else {
            let before = Self::ViewType::new(self.label, EMPTY_Y, EMPTY_Z);
            let after = Self::ViewType::new(self.label, EMPTY_Y, EMPTY_Z);
            (before, after)
        }
    }

    fn split_at_time(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        if let Some(point) = self.find_time(point).0 {
            self.split_at(point)
        } else {
            let before = Self::ViewType::new(self.label, EMPTY_Y, EMPTY_Z);
            let after = Self::ViewType::new(self.label, EMPTY_Y, EMPTY_Z);
            (before, after)
        }
    }

    fn slice<I: RangeBounds<usize> + Clone>(&'a self, bounds: I) -> Self::ViewType {
        let start = match bounds.start_bound() {
            Bound::Included(i) => *i,
            Bound::Excluded(i) => *i,
            Bound::Unbounded => 0,
        };
        let end = match bounds.end_bound() {
            Bound::Included(i) => *i + 1,
            Bound::Excluded(i) => *i,
            Bound::Unbounded => self.len(),
        };
        let y = &self.y[start..end];
        let z = &self.z[start..end];
        Self::ViewType::new(self.label, y, z)
    }
}

impl<'a, X, Y> SplittableFeatureLike<'a, X, Y> for SimpleFeatureView<'a, X, Y> {
    type ViewType = Self;

    fn split_at_time(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        if let Some(point) = self.find_time(point).0 {
            self.split_at(point)
        } else {
            let before = Self::ViewType::new(self.label, EMPTY_Y, EMPTY_Z);
            let after = Self::ViewType::new(self.label, EMPTY_Y, EMPTY_Z);
            (before, after)
        }
    }

    fn slice<I: RangeBounds<usize> + Clone>(&'a self, bounds: I) -> Self::ViewType {
        let start = bounds.start_bound();
        let end = bounds.end_bound();
        match (start, end) {
            (Bound::Included(i), Bound::Included(j)) => {
                Self::ViewType::new(self.label, &self.y[*i..=*j], &self.z[*i..=*j])
            }
            (Bound::Included(i), Bound::Excluded(j)) => {
                Self::ViewType::new(self.label, &self.y[*i..*j], &self.z[*i..*j])
            }
            (Bound::Included(i), Bound::Unbounded) => {
                Self::ViewType::new(self.label, &self.y[*i..], &self.z[*i..])
            }
            (Bound::Excluded(i), Bound::Included(j)) => {
                Self::ViewType::new(self.label, &self.y[*i..=*j], &self.z[*i..=*j])
            }
            (Bound::Excluded(i), Bound::Excluded(j)) => {
                Self::ViewType::new(self.label, &self.y[*i..*j], &self.z[*i..*j])
            }
            (Bound::Excluded(i), Bound::Unbounded) => {
                Self::ViewType::new(self.label, &self.y[*i..], &self.z[*i..])
            }
            (Bound::Unbounded, Bound::Included(j)) => {
                Self::ViewType::new(self.label, &self.y[..=*j], &self.z[..=*j])
            }
            (Bound::Unbounded, Bound::Excluded(j)) => {
                Self::ViewType::new(self.label, &self.y[..*j], &self.z[..*j])
            }
            (Bound::Unbounded, Bound::Unbounded) => Self::ViewType::new(self.label, self.y, self.z),
        }
    }

    fn split_at(&'a self, index: usize) -> (Self::ViewType, Self::ViewType) {
        if !self.is_empty() {
            let before = Self::ViewType::new(self.label, &self.y[..index], &self.z[..index]);
            let after = Self::ViewType::new(self.label, &self.y[index..], &self.z[index..]);
            (before, after)
        } else {
            let before = Self::ViewType::new(self.label, EMPTY_Y, EMPTY_Z);
            let after = Self::ViewType::new(self.label, EMPTY_Y, EMPTY_Z);
            (before, after)
        }
    }
}

impl<'a, X, Y> TimeArray<Y> for SimpleFeatureView<'a, X, Y> {
    fn time_view(&self) -> &[f64] {
        self.y
    }

    fn intensity_view(&self) -> &[f32] {
        &self.z
    }
}
