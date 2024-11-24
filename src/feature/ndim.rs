use core::slice;
use std::marker::PhantomData;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{coordinate::CoordinateLike, IntensityMeasurement};
use crate::{IntensityMeasurementMut, Mass, MZ};

use super::traits::{CoArrayOps, FeatureLikeNDLike, FeatureLikeNDLikeMut};
use super::util::NonNan;
use super::{TimeArray, TimeInterval};

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NDFeature<T, Y> {
    dimensions: Vec<Vec<f64>>,
    time: Vec<f64>,
    intensity: Vec<f32>,
    _t: PhantomData<T>,
    _y: PhantomData<Y>,
}

impl<T, Y> Default for NDFeature<T, Y> {
    fn default() -> Self {
        Self::new(Vec::new(), Vec::new(), Vec::new())
    }
}

impl<T, Y> NDFeature<T, Y> {
    pub fn new(dimensions: Vec<Vec<f64>>, time: Vec<f64>, intensity: Vec<f32>) -> Self {
        Self {
            dimensions,
            time,
            intensity,
            _t: PhantomData,
            _y: PhantomData,
        }
    }

    #[allow(unused)]
    pub fn coordinate_dim(&self, i: usize) -> f64 {
        self.weighted_average(&self.dimensions[i], &self.intensity)
    }

    fn coordinate_y(&self) -> f64 {
        self.weighted_average(&self.time, &self.intensity)
    }

    pub unsafe fn push_raw_unchecked(&mut self, point: NDPoint<T, Y>) {
        self.time.push(point.time);
        self.intensity.push(point.intensity);
        if self.dimensions.is_empty() {
            for _ in 0..point.dimensions.len() {
                self.dimensions.push(Vec::with_capacity(2));
            }
        }
        for (d, p) in self.dimensions.iter_mut().zip(point.dimensions.into_iter()) {
            d.push(p)
        }
    }

    pub(crate) fn sort_by_time(&mut self) {
        let n = self.len();
        let mut mask: Vec<_> = (0..n).collect();
        mask.sort_by_key(|i| NonNan::new(self.time[*i]));

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

                    self.time.swap(current_idx, next_idx);
                    self.intensity.swap(current_idx, next_idx);
                    for dim in self.dimensions.iter_mut() {
                        dim.swap(current_idx, next_idx);
                    }
                    current_idx = next_idx;
                }
            }
        }
    }
}

impl<T, Y> CoArrayOps for NDFeature<T, Y> {}

impl<T, Y> TimeInterval<Y> for NDFeature<T, Y> {
    fn start_time(&self) -> Option<f64> {
        self.time.first().copied()
    }

    fn end_time(&self) -> Option<f64> {
        self.time.last().copied()
    }

    fn apex_time(&self) -> Option<f64> {
        self.apex_of(&self.time, &self.intensity)
    }

    fn area(&self) -> f32 {
        self.trapezoid_integrate(&self.time, &self.intensity)
    }

    fn iter_time(&self) -> impl Iterator<Item = f64> {
        self.time.iter().copied()
    }
}

impl<T, Y> TimeArray<Y> for NDFeature<T, Y> {
    fn time_view(&self) -> &[f64] {
        &self.time
    }

    fn intensity_view(&self) -> &[f32] {
        &self.intensity
    }
}

impl<T, Y> IntensityMeasurement for NDFeature<T, Y> {
    fn intensity(&self) -> f32 {
        self.intensity.iter().sum()
    }
}

impl<T, Y> FeatureLikeNDLike<T, Y> for NDFeature<T, Y> {
    type Point = NDPoint<T, Y>;

    fn len(&self) -> usize {
        self.time.len()
    }

    fn iter(&self) -> impl Iterator<Item = Self::Point> {
        NDIter::new(
            self.dimensions.iter().map(|d| d.iter()).collect(),
            self.time.iter(),
            self.intensity.iter(),
        )
    }
}

macro_rules! impl_dim_feat {
    ($pt:ident, $dim:ty, $idx:literal) => {
        impl $crate::coordinate::CoordinateLike<$dim>
            for $pt<($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                self.weighted_average(&self.dimensions[$idx], &self.intensity)
            }
        }

        impl $crate::coordinate::CoordinateLike<$crate::IonMobility>
            for $pt<($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                self.weighted_average(&self.dimensions[$idx + 1], &self.intensity)
            }
        }

        impl $crate::coordinate::CoordinateLike<$dim>
            for $pt<($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                self.weighted_average(&self.dimensions[$idx], &self.intensity)
            }
        }

        impl $crate::coordinate::CoordinateLike<$crate::Time>
            for $pt<($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                self.weighted_average(&self.dimensions[$idx + 1], &self.intensity)
            }
        }
    };
}

impl<T, Y> PartialEq for NDFeature<T, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.dimensions == other.dimensions
            && self.time == other.time
            && self.intensity == other.intensity
            && self._t == other._t
            && self._y == other._y
    }
}

impl<T, Y> PartialOrd for NDFeature<T, Y> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        for (sd, od) in self.dimensions.iter().zip(other.dimensions.iter()) {
            let comp = self
                .weighted_average(sd, &self.intensity)
                .total_cmp(&other.weighted_average(&od, &other.intensity));
            if !comp.is_eq() {
                return Some(comp);
            }
        }
        match self.coordinate_y().total_cmp(&other.coordinate_y()) {
            std::cmp::Ordering::Equal => {}
            ord => return Some(ord),
        }
        match self.intensity.partial_cmp(&other.intensity) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self._t.partial_cmp(&other._t) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        self._y.partial_cmp(&other._y)
    }
}

impl_dim_feat!(NDFeature, MZ, 0);
impl_dim_feat!(NDFeature, Mass, 0);

#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NDPoint<T, Y> {
    pub dimensions: Vec<f64>,
    time: f64,
    intensity: f32,
    _t: PhantomData<T>,
    _y: PhantomData<Y>,
}

impl<T, Y> NDPoint<T, Y> {
    pub fn new(dims: Vec<f64>, time: f64, intensity: f32) -> Self {
        Self {
            dimensions: dims,
            time,
            intensity,
            _t: PhantomData,
            _y: PhantomData,
        }
    }
}

impl<T, Y> PartialEq for NDPoint<T, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.dimensions == other.dimensions
            && self.time == other.time
            && self.intensity == other.intensity
            && self._t == other._t
            && self._y == other._y
    }
}

impl<T, Y> PartialOrd for NDPoint<T, Y> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.dimensions.partial_cmp(&other.dimensions) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self.time.partial_cmp(&other.time) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self.intensity.partial_cmp(&other.intensity) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self._t.partial_cmp(&other._t) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        self._y.partial_cmp(&other._y)
    }
}

impl<T, Y> CoordinateLike<Y> for NDPoint<T, Y> {
    fn coordinate(&self) -> f64 {
        self.time
    }
}

impl<T, Y> IntensityMeasurement for NDPoint<T, Y> {
    fn intensity(&self) -> f32 {
        self.intensity
    }
}

pub struct NDIter<'a, T, Y> {
    dims_iter: Vec<slice::Iter<'a, f64>>,
    time_iter: slice::Iter<'a, f64>,
    intensity_iter: slice::Iter<'a, f32>,
    _t: PhantomData<T>,
    _y: PhantomData<Y>,
}

impl<'a, T, Y> NDIter<'a, T, Y> {
    pub fn new(
        dims_iter: Vec<slice::Iter<'a, f64>>,
        time_iter: slice::Iter<'a, f64>,
        intensity_iter: slice::Iter<'a, f32>,
    ) -> Self {
        Self {
            dims_iter,
            time_iter,
            intensity_iter,
            _t: PhantomData,
            _y: PhantomData,
        }
    }

    fn read_next(&mut self) -> Option<NDPoint<T, Y>> {
        let time = *self.time_iter.next()?;
        let intensity = *self.intensity_iter.next()?;
        let dims = self
            .dims_iter
            .iter_mut()
            .map(|i| *i.next().unwrap())
            .collect();
        Some(NDPoint::new(dims, time, intensity))
    }
}

impl<'a, T, Y> Iterator for NDIter<'a, T, Y> {
    type Item = NDPoint<T, Y>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_next()
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        let time = *self.time_iter.nth(n)?;
        let intensity = *self.intensity_iter.nth(n)?;
        let dims = self
            .dims_iter
            .iter_mut()
            .map(|i| *i.nth(n).unwrap())
            .collect();
        Some(NDPoint::new(dims, time, intensity))
    }
}

macro_rules! impl_dim_pt {
    ($pt:ident, $dim:ty, $idx:literal) => {
        impl $crate::coordinate::CoordinateLike<$dim>
            for $pt<($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                self.dimensions[$idx]
            }
        }

        impl $crate::coordinate::CoordinateLike<$crate::IonMobility>
            for $pt<($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                self.dimensions[$idx + 1]
            }
        }

        impl $crate::coordinate::CoordinateLike<$dim>
            for $pt<($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                self.dimensions[$idx]
            }
        }

        impl $crate::coordinate::CoordinateLike<$crate::Time>
            for $pt<($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                self.dimensions[$idx + 1]
            }
        }
    };
}

impl_dim_pt!(NDPoint, MZ, 0);
impl_dim_pt!(NDPoint, Mass, 0);

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(Serialize))]
pub struct NDPointMutRef<'a, T, Y> {
    pub dimensions: Vec<&'a mut f64>,
    time: f64,
    intensity: &'a mut f32,
    _t: PhantomData<T>,
    _y: PhantomData<Y>,
}

impl<'a, T, Y> NDPointMutRef<'a, T, Y> {
    pub fn new(dimensions: Vec<&'a mut f64>, time: f64, intensity: &'a mut f32) -> Self {
        Self {
            dimensions,
            time,
            intensity,
            _t: PhantomData,
            _y: PhantomData,
        }
    }
}

impl<'a, T, Y> From<NDPointMutRef<'a, T, Y>> for NDPoint<T, Y> {
    fn from(value: NDPointMutRef<'a, T, Y>) -> Self {
        Self::new(value.dimensions.into_iter().map(|d| *d).collect(), value.time, *value.intensity)
    }
}

macro_rules! impl_dim_pt_ref {
    ($pt:ident, $dim:ty, $idx:literal) => {
        impl<'a> $crate::coordinate::CoordinateLike<$dim>
            for $pt<'a, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                *self.dimensions[$idx]
            }
        }

        impl<'a> $crate::coordinate::CoordinateLikeMut<$dim>
            for $pt<'a, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate_mut(&mut self) -> &mut f64 {
                self.dimensions[$idx]
            }
        }

        impl<'a> $crate::coordinate::CoordinateLike<$crate::IonMobility>
            for $pt<'a, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                *self.dimensions[$idx + 1]
            }
        }

        impl<'a> $crate::coordinate::CoordinateLike<$dim>
            for $pt<'a, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                *self.dimensions[$idx]
            }
        }

        impl<'a> $crate::coordinate::CoordinateLikeMut<$dim>
            for $pt<'a, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate_mut(&mut self) -> &mut f64 {
                self.dimensions[$idx]
            }
        }

        impl<'a> $crate::coordinate::CoordinateLike<$crate::Time>
            for $pt<'a, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                *self.dimensions[$idx + 1]
            }
        }
    };
}

impl_dim_pt_ref!(NDPointMutRef, MZ, 0);
impl_dim_pt_ref!(NDPointMutRef, Mass, 0);

impl<'a, T, Y> PartialEq for NDPointMutRef<'a, T, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.dimensions == other.dimensions
            && self.time == other.time
            && self.intensity == other.intensity
            && self._t == other._t
            && self._y == other._y
    }
}

impl<'a, T, Y> PartialOrd for NDPointMutRef<'a, T, Y> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.dimensions.partial_cmp(&other.dimensions) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self.time.partial_cmp(&other.time) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self.intensity.partial_cmp(&other.intensity) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self._t.partial_cmp(&other._t) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        self._y.partial_cmp(&other._y)
    }
}

impl<'a, T, Y> CoordinateLike<Y> for NDPointMutRef<'a, T, Y> {
    fn coordinate(&self) -> f64 {
        self.time
    }
}

impl<'a, T, Y> IntensityMeasurement for NDPointMutRef<'a, T, Y> {
    fn intensity(&self) -> f32 {
        *self.intensity
    }
}

impl<'a, T, Y> IntensityMeasurementMut for NDPointMutRef<'a, T, Y> {
    fn intensity_mut(&mut self) -> &mut f32 {
        &mut self.intensity
    }
}

impl<'a, T, Y> FeatureLikeNDLikeMut<'a, T, Y> for NDFeature<T, Y>
where
    NDFeature<T, Y>: 'a,
{
    type PointMutRef = NDPointMutRef<'a, T, Y>;

    fn iter_mut(&'a mut self) -> impl Iterator<Item = Self::PointMutRef> {
        NDIterMut::new(
            self.dimensions.iter_mut().map(|d| d.iter_mut()).collect(),
            self.time.iter(),
            self.intensity.iter_mut(),
        )
    }

    fn push<X: Into<Self::Point>>(&mut self, pt: X, time: f64) {
        let mut pt: Self::Point = pt.into();
        pt.time = time;
        self.push_raw(pt);
    }

    fn push_raw(&mut self, point: Self::Point) {
        let needs_sort = !self.is_empty() && point.time < self.time.last().copied().unwrap();
        unsafe { self.push_raw_unchecked(point) };
        if needs_sort {
            self.sort_by_time();
        }
    }
}

impl<T, Y> FromIterator<NDPoint<T, Y>> for NDFeature<T, Y> {
    fn from_iter<I: IntoIterator<Item = NDPoint<T, Y>>>(iter: I) -> Self {
        let mut this = Self::default();
        for pt in iter {
            this.push_raw(pt);
        }
        this
    }
}

impl<T, Y> Extend<NDPoint<T, Y>> for NDFeature<T, Y> {
    fn extend<I: IntoIterator<Item = NDPoint<T, Y>>>(&mut self, iter: I) {
        for pt in iter {
            self.push_raw(pt);
        }
    }
}

pub struct NDIterMut<'a, T, Y> {
    dims_iter: Vec<slice::IterMut<'a, f64>>,
    time_iter: slice::Iter<'a, f64>,
    intensity_iter: slice::IterMut<'a, f32>,
    _t: PhantomData<T>,
    _y: PhantomData<Y>,
}

impl<'a, T, Y> NDIterMut<'a, T, Y> {
    pub fn new(
        dims_iter: Vec<slice::IterMut<'a, f64>>,
        time_iter: slice::Iter<'a, f64>,
        intensity_iter: slice::IterMut<'a, f32>,
    ) -> Self {
        Self {
            dims_iter,
            time_iter,
            intensity_iter,
            _t: PhantomData,
            _y: PhantomData,
        }
    }

    fn read_next(&mut self) -> Option<NDPointMutRef<'a, T, Y>> {
        let time = *self.time_iter.next()?;
        let intensity = self.intensity_iter.next()?;
        let dims = self
            .dims_iter
            .iter_mut()
            .map(|i| i.next().unwrap())
            .collect();
        Some(NDPointMutRef::new(dims, time, intensity))
    }
}

impl<'a, T, Y> Iterator for NDIterMut<'a, T, Y> {
    type Item = NDPointMutRef<'a, T, Y>;

    fn next(&mut self) -> Option<Self::Item> {
        self.read_next()
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        let time = *self.time_iter.nth(n)?;
        let intensity = self.intensity_iter.nth(n)?;
        let dims = self
            .dims_iter
            .iter_mut()
            .map(|i| i.nth(n).unwrap())
            .collect();
        Some(NDPointMutRef::new(dims, time, intensity))
    }
}

#[cfg(test)]
mod test {
    use crate::{
        coordinate::{IonMobilityLocated, TimeLocated},
        IonMobility, MZLocated, Time,
    };

    use super::*;

    #[test]
    fn test_create() {
        let pt = NDPoint::<(MZ, IonMobility), Time>::new(vec![50.0, 0.5], 12.6, 1000.0);
        assert_eq!(pt.time(), 12.6);
        assert_eq!(pt.mz(), 50.0);
        assert_eq!(pt.ion_mobility(), 0.5);
    }

    #[test]
    fn test_create_feature() {
        let points: Vec<NDPoint<(MZ, IonMobility), Time>> = vec![
            NDPoint::new(vec![50.0, 0.5], 12.6, 1000.0),
            NDPoint::new(vec![50.01, 0.49], 12.7, 1200.0),
            NDPoint::new(vec![50.0, 0.51], 12.9, 800.0),
        ];

        let f: NDFeature<_, _> = points.into_iter().collect();
        assert_eq!(f.len(), 3);
        assert!(
            (f.ion_mobility() - 0.498666).abs() < 1e-3,
            "ion_mobility was {}",
            f.ion_mobility()
        );
        assert!((f.mz() - 50.004).abs() < 1e-3, "mz was {}", f.mz());
        assert!(
            (f.apex_time().unwrap() - 12.7).abs() < 1e-3,
            "apex_time was {}",
            f.apex_time().unwrap()
        );
    }
}
