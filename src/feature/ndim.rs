use std::{
    iter::FusedIterator,
    marker::PhantomData,
    ops::{Index, IndexMut},
};

use crate::{coordinate::CoordinateLike, IntensityMeasurement, KnownChargeMut};
use crate::{
    feature::{AsPeakIter, BuildFromPeak},
    peak::{IonMobilityAwareCentroidPeak, IonMobilityAwareDeconvolutedPeak},
    prelude::IonMobilityLocated,
    CoordinateLikeMut, IntensityMeasurementMut, IonMobility, KnownCharge, MZLocated, Mass,
    MassLocated, Time, MZ,
};
#[cfg(feature = "serde")]
use {
    serde::{Deserialize, Serialize},
    serde_big_array::BigArray,
};

use super::{
    charged::Charged,
    traits::{CoArrayOps, NDFeatureLike, NDFeatureLikeMut},
    FeatureLikeMut,
};
use super::{util::NonNan, FeatureLike};
use super::{TimeArray, TimeInterval};

macro_rules! impl_dim_pt {
    ($pt:ident, $dim:ty, $idx:literal) => {
        impl $crate::coordinate::CoordinateLike<$dim>
            for $pt<2, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                self.dimensions[$idx]
            }
        }

        impl $crate::coordinate::CoordinateLike<$crate::IonMobility>
            for $pt<2, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                self.dimensions[$idx + 1]
            }
        }

        impl $crate::coordinate::CoordinateLike<$dim>
            for $pt<2, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                self.dimensions[$idx]
            }
        }

        impl $crate::coordinate::CoordinateLike<$crate::Time>
            for $pt<2, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                self.dimensions[$idx + 1]
            }
        }

        //

        impl $crate::coordinate::CoordinateLikeMut<$dim>
            for $pt<2, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate_mut(&mut self) -> &mut f64 {
                &mut self.dimensions[$idx]
            }
        }

        impl $crate::coordinate::CoordinateLikeMut<$crate::IonMobility>
            for $pt<2, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate_mut(&mut self) -> &mut f64 {
                &mut self.dimensions[$idx + 1]
            }
        }

        impl $crate::coordinate::CoordinateLikeMut<$dim>
            for $pt<2, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate_mut(&mut self) -> &mut f64 {
                &mut self.dimensions[$idx]
            }
        }

        impl $crate::coordinate::CoordinateLikeMut<$crate::Time>
            for $pt<2, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate_mut(&mut self) -> &mut f64 {
                &mut self.dimensions[$idx + 1]
            }
        }
    };
}

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NDPoint<const N: usize, T, Y> {
    #[cfg_attr(feature = "serde", serde(with = "BigArray"))]
    pub dimensions: [f64; N],
    pub time: f64,
    pub intensity: f32,
    #[cfg_attr(feature = "serde", serde(skip))]
    _t: PhantomData<T>,
    #[cfg_attr(feature = "serde", serde(skip))]
    _y: PhantomData<Y>,
}

impl<const N: usize, T, Y> Index<usize> for NDPoint<N, T, Y> {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        &self.dimensions[index]
    }
}

impl<const N: usize, T, Y> IndexMut<usize> for NDPoint<N, T, Y> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.dimensions[index]
    }
}

impl<const N: usize, T, Y> NDPoint<N, T, Y> {
    pub fn new(dimensions: [f64; N], time: f64, intensity: f32) -> Self {
        Self {
            dimensions,
            time,
            intensity,
            _t: PhantomData,
            _y: PhantomData,
        }
    }

    pub fn get(&self, index: usize) -> Option<f64> {
        self.dimensions.get(index).copied()
    }
}

impl<const N: usize, T, Y> IntensityMeasurement for NDPoint<N, T, Y> {
    fn intensity(&self) -> f32 {
        self.intensity
    }
}

impl<const N: usize, T, Y> IntensityMeasurementMut for NDPoint<N, T, Y> {
    fn intensity_mut(&mut self) -> &mut f32 {
        &mut self.intensity
    }
}

impl<const N: usize, T, Y> PartialEq for NDPoint<N, T, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.dimensions == other.dimensions
            && self.time == other.time
            && self.intensity == other.intensity
            && self._t == other._t
            && self._y == other._y
    }
}

impl<const N: usize, T, Y> PartialOrd for NDPoint<N, T, Y> {
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

impl<const N: usize, T, Y> CoordinateLike<Y> for NDPoint<N, T, Y> {
    fn coordinate(&self) -> f64 {
        self.time
    }
}

impl<const N: usize, T, Y> CoordinateLikeMut<Y> for NDPoint<N, T, Y> {
    fn coordinate_mut(&mut self) -> &mut f64 {
        &mut self.time
    }
}

impl<const N: usize, T, Y> Default for NDPoint<N, T, Y> {
    fn default() -> Self {
        Self {
            dimensions: [0.0; N],
            time: Default::default(),
            intensity: Default::default(),
            _t: Default::default(),
            _y: Default::default(),
        }
    }
}

impl_dim_pt!(NDPoint, MZ, 0);
impl_dim_pt!(NDPoint, Mass, 0);

#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(Serialize))]
pub struct NDPointMutRef<'a, const N: usize, T, Y> {
    #[cfg_attr(feature = "serde", serde(with = "BigArray"))]
    pub dimensions: [&'a mut f64; N],
    time: f64,
    intensity: &'a mut f32,
    #[cfg_attr(feature = "serde", serde(skip))]
    _t: PhantomData<T>,
    #[cfg_attr(feature = "serde", serde(skip))]
    _y: PhantomData<Y>,
}

impl<'a, const N: usize, T, Y> Index<usize> for NDPointMutRef<'a, N, T, Y> {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        self.dimensions[index]
    }
}

impl<'a, const N: usize, T, Y> IndexMut<usize> for NDPointMutRef<'a, N, T, Y> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.dimensions[index]
    }
}

impl<'a, const N: usize, T, Y> NDPointMutRef<'a, N, T, Y> {
    pub fn new(dimensions: [&'a mut f64; N], time: f64, intensity: &'a mut f32) -> Self {
        Self {
            dimensions,
            time,
            intensity,
            _t: PhantomData,
            _y: PhantomData,
        }
    }
}

impl<'a, const N: usize, T, Y> From<NDPointMutRef<'a, N, T, Y>> for NDPoint<N, T, Y> {
    fn from(value: NDPointMutRef<'a, N, T, Y>) -> Self {
        Self::new(value.dimensions.map(|d| *d), value.time, *value.intensity)
    }
}

impl<'a, const N: usize, T, Y> PartialEq for NDPointMutRef<'a, N, T, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.dimensions == other.dimensions
            && self.time == other.time
            && self.intensity == other.intensity
            && self._t == other._t
            && self._y == other._y
    }
}

impl<'a, const N: usize, T, Y> PartialOrd for NDPointMutRef<'a, N, T, Y> {
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

impl<'a, const N: usize, T, Y> CoordinateLike<Y> for NDPointMutRef<'a, N, T, Y> {
    fn coordinate(&self) -> f64 {
        self.time
    }
}

impl<'a, const N: usize, T, Y> CoordinateLikeMut<Y> for NDPointMutRef<'a, N, T, Y> {
    fn coordinate_mut(&mut self) -> &mut f64 {
        &mut self.time
    }
}

impl<'a, const N: usize, T, Y> IntensityMeasurement for NDPointMutRef<'a, N, T, Y> {
    fn intensity(&self) -> f32 {
        *self.intensity
    }
}

impl<'a, const N: usize, T, Y> IntensityMeasurementMut for NDPointMutRef<'a, N, T, Y> {
    fn intensity_mut(&mut self) -> &mut f32 {
        &mut self.intensity
    }
}

macro_rules! impl_dim_pt_ref {
    ($pt:ident, $dim:ty, $idx:literal) => {
        impl<'a> $crate::coordinate::CoordinateLike<$dim>
            for $pt<'a, 2, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                *self.dimensions[$idx]
            }
        }

        impl<'a> $crate::coordinate::CoordinateLikeMut<$dim>
            for $pt<'a, 2, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate_mut(&mut self) -> &mut f64 {
                self.dimensions[$idx]
            }
        }

        impl<'a> $crate::coordinate::CoordinateLike<$crate::IonMobility>
            for $pt<'a, 2, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                *self.dimensions[$idx + 1]
            }
        }

        impl<'a> $crate::coordinate::CoordinateLike<$dim>
            for $pt<'a, 2, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                *self.dimensions[$idx]
            }
        }

        impl<'a> $crate::coordinate::CoordinateLikeMut<$dim>
            for $pt<'a, 2, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate_mut(&mut self) -> &mut f64 {
                self.dimensions[$idx]
            }
        }

        impl<'a> $crate::coordinate::CoordinateLike<$crate::Time>
            for $pt<'a, 2, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                *self.dimensions[$idx + 1]
            }
        }
    };
}

impl_dim_pt_ref!(NDPointMutRef, MZ, 0);
impl_dim_pt_ref!(NDPointMutRef, Mass, 0);

#[derive(Debug)]
pub struct NDIter<'a, const N: usize, T, Y> {
    dim_iters: [core::slice::Iter<'a, f64>; N],
    time_iter: core::slice::Iter<'a, f64>,
    intensity_iter: core::slice::Iter<'a, f32>,
    _t: PhantomData<(T, Y)>,
}

impl<'a, const N: usize, T, Y> DoubleEndedIterator for NDIter<'a, N, T, Y> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if let Some(time) = self.time_iter.next_back().copied() {
            let intensity = self.intensity_iter.next_back().copied().unwrap();
            let dimensions = self
                .dim_iters
                .each_mut()
                .map(|i| i.next_back().copied().unwrap());
            Some(NDPoint::new(dimensions, time, intensity))
        } else {
            None
        }
    }
}

impl<'a, const N: usize, T, Y> FusedIterator for NDIter<'a, N, T, Y> {}

impl<'a, const N: usize, T, Y> ExactSizeIterator for NDIter<'a, N, T, Y> {
    fn len(&self) -> usize {
        self.time_iter.len()
    }
}

impl<'a, const N: usize, T, Y> Iterator for NDIter<'a, N, T, Y> {
    type Item = NDPoint<N, T, Y>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(time) = self.time_iter.next().copied() {
            let intensity = self.intensity_iter.next().copied().unwrap();
            let dimensions = self
                .dim_iters
                .each_mut()
                .map(|i| i.next().copied().unwrap());
            Some(NDPoint::new(dimensions, time, intensity))
        } else {
            None
        }
    }
}

impl<'a, const N: usize, T, Y> NDIter<'a, N, T, Y> {
    pub fn new(
        dim_iters: [core::slice::Iter<'a, f64>; N],
        time_iter: core::slice::Iter<'a, f64>,
        intensity_iter: core::slice::Iter<'a, f32>,
    ) -> Self {
        Self {
            dim_iters,
            time_iter,
            intensity_iter,
            _t: PhantomData,
        }
    }
}

#[derive(Debug)]
pub struct NDIterMut<'a, const N: usize, T, Y> {
    dim_iters: [core::slice::IterMut<'a, f64>; N],
    time_iter: core::slice::Iter<'a, f64>,
    intensity_iter: core::slice::IterMut<'a, f32>,
    _t: PhantomData<(T, Y)>,
}

impl<'a, const N: usize, T, Y> DoubleEndedIterator for NDIterMut<'a, N, T, Y> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if let Some(time) = self.time_iter.next_back().copied() {
            let intensity = self.intensity_iter.next_back().unwrap();
            let dimensions = self.dim_iters.each_mut().map(|i| i.next_back().unwrap());
            Some(Self::Item::new(dimensions, time, intensity))
        } else {
            None
        }
    }
}

impl<'a, const N: usize, T, Y> FusedIterator for NDIterMut<'a, N, T, Y> {}

impl<'a, const N: usize, T, Y> ExactSizeIterator for NDIterMut<'a, N, T, Y> {
    fn len(&self) -> usize {
        self.time_iter.len()
    }
}

impl<'a, const N: usize, T, Y> Iterator for NDIterMut<'a, N, T, Y> {
    type Item = NDPointMutRef<'a, N, T, Y>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(time) = self.time_iter.next().copied() {
            let intensity = self.intensity_iter.next().unwrap();
            let dimensions = self.dim_iters.each_mut().map(|i| i.next().unwrap());
            Some(Self::Item::new(dimensions, time, intensity))
        } else {
            None
        }
    }
}

impl<'a, const N: usize, T, Y> NDIterMut<'a, N, T, Y> {
    pub fn new(
        dim_iters: [core::slice::IterMut<'a, f64>; N],
        time_iter: core::slice::Iter<'a, f64>,
        intensity_iter: core::slice::IterMut<'a, f32>,
    ) -> Self {
        Self {
            dim_iters,
            time_iter,
            intensity_iter,
            _t: PhantomData,
        }
    }
}

macro_rules! impl_dim_feat {
    ($pt:ident, $dim:ty, $idx:literal) => {
        impl $crate::coordinate::CoordinateLike<$dim>
            for $pt<2, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                self.weighted_average(&self.dimensions[$idx], &self.intensity)
            }
        }

        impl $crate::coordinate::CoordinateLike<$crate::IonMobility>
            for $pt<2, ($dim, $crate::coordinate::IonMobility), $crate::coordinate::Time>
        {
            fn coordinate(&self) -> f64 {
                self.weighted_average(&self.dimensions[$idx + 1], &self.intensity)
            }
        }

        impl $crate::coordinate::CoordinateLike<$dim>
            for $pt<2, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                self.weighted_average(&self.dimensions[$idx], &self.intensity)
            }
        }

        impl $crate::coordinate::CoordinateLike<$crate::Time>
            for $pt<2, ($dim, $crate::coordinate::Time), $crate::coordinate::IonMobility>
        {
            fn coordinate(&self) -> f64 {
                self.weighted_average(&self.dimensions[$idx + 1], &self.intensity)
            }
        }
    };
}

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NDFeature<const N: usize, T, Y> {
    #[cfg_attr(feature = "serde", serde(with = "BigArray"))]
    dimensions: [Vec<f64>; N],
    time: Vec<f64>,
    intensity: Vec<f32>,
    #[cfg_attr(feature = "serde", serde(skip))]
    _t: PhantomData<T>,
    #[cfg_attr(feature = "serde", serde(skip))]
    _y: PhantomData<Y>,
}

impl<const N: usize, T, Y> NDFeatureLikeMut<T, Y> for NDFeature<N, T, Y> {
    type PointMutRef<'a>
        = NDPointMutRef<'a, N, T, Y>
    where
        Self: 'a;

    fn iter_mut(&mut self) -> impl Iterator<Item = Self::PointMutRef<'_>> {
        let dim_iter = self.dimensions.each_mut().map(|i| i.iter_mut());
        let time_iter = self.time.iter();
        let intensity_iter = self.intensity.iter_mut();
        NDIterMut::new(dim_iter, time_iter, intensity_iter)
    }

    fn push<X: Into<Self::Point>>(&mut self, pt: X, time: f64) {
        let mut pt = pt.into();
        pt.time = time;
        self.push_raw(pt);
    }

    fn push_raw(&mut self, point: Self::Point) {
        self.push(point);
    }
}

impl<const N: usize, T, Y> PartialOrd for NDFeature<N, T, Y> {
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

impl<const N: usize, T, Y> PartialEq for NDFeature<N, T, Y> {
    fn eq(&self, other: &Self) -> bool {
        self.dimensions == other.dimensions
            && self.time == other.time
            && self.intensity == other.intensity
            && self._t == other._t
            && self._y == other._y
    }
}

impl_dim_feat!(NDFeature, MZ, 0);
impl_dim_feat!(NDFeature, Mass, 0);

impl<const N: usize, T, Y> IntensityMeasurement for NDFeature<N, T, Y> {
    fn intensity(&self) -> f32 {
        self.intensity.iter().sum()
    }
}

impl<const N: usize, T, Y> NDFeatureLike<T, Y> for NDFeature<N, T, Y> {
    type Point = NDPoint<N, T, Y>;

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> impl Iterator<Item = Self::Point> {
        let dim_iter = self.dimensions.each_ref().map(|i| i.iter());
        let time_iter = self.time.iter();
        let intensity_iter = self.intensity.iter();
        NDIter::new(dim_iter, time_iter, intensity_iter)
    }

    fn coordinate(&self) -> Self::Point {
        let dims = self
            .dimensions
            .each_ref()
            .map(|d| self.weighted_average(d.as_slice(), &self.intensity));
        let time = self.apex_time().unwrap_or_default();
        let intensity = self.intensity();
        Self::Point::new(dims, time, intensity)
    }
}

impl<const N: usize, T, Y> NDFeature<N, T, Y> {
    pub fn new(dimensions: [Vec<f64>; N], time: Vec<f64>, intensity: Vec<f32>) -> Self {
        Self {
            dimensions,
            time,
            intensity,
            _t: PhantomData,
            _y: PhantomData,
        }
    }

    pub fn into_inner(self) -> ([Vec<f64>; N], Vec<f64>, Vec<f32>) {
        (self.dimensions, self.time, self.intensity)
    }

    fn iter(&self) -> NDIter<'_, N, T, Y> {
        let dim_iter = self.dimensions.each_ref().map(|i| i.iter());
        let time_iter = self.time.iter();
        let intensity_iter = self.intensity.iter();
        NDIter::new(dim_iter, time_iter, intensity_iter)
    }

    pub fn as_inner(&self) -> ([&[f64]; N], &[f64], &[f32]) {
        (
            self.dimensions.each_ref().map(|i| i.as_slice()),
            &self.time,
            &self.intensity,
        )
    }

    pub fn len(&self) -> usize {
        self.time.len()
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

    pub fn push(&mut self, point: NDPoint<N, T, Y>) {
        let needs_sort = if let Some(t) = self.time.last().copied() {
            point.time < t
        } else {
            false
        };

        self.time.push(point.time);
        self.intensity.push(point.intensity);
        for (dim, val) in self.dimensions.iter_mut().zip(point.dimensions.into_iter()) {
            dim.push(val)
        }

        if needs_sort {
            self.sort_by_time();
        }
    }
}

impl<const N: usize, T, Y> CoArrayOps for NDFeature<N, T, Y> {}

impl<const N: usize, T, Y> Default for NDFeature<N, T, Y> {
    fn default() -> Self {
        Self {
            dimensions: core::array::from_fn(|_| Vec::new()),
            time: Default::default(),
            intensity: Default::default(),
            _t: Default::default(),
            _y: Default::default(),
        }
    }
}

impl<const N: usize, T, Y> TimeInterval<Y> for NDFeature<N, T, Y> {
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

impl<const N: usize, T, Y> TimeArray<Y> for NDFeature<N, T, Y> {
    fn time_view(&self) -> &[f64] {
        &self.time
    }

    fn intensity_view(&self) -> &[f32] {
        &self.intensity
    }
}

impl<const N: usize, T, Y> FromIterator<NDPoint<N, T, Y>> for NDFeature<N, T, Y> {
    fn from_iter<I: IntoIterator<Item = NDPoint<N, T, Y>>>(iter: I) -> Self {
        let mut this = Self::default();
        for pt in iter {
            this.push_raw(pt);
        }
        this
    }
}

impl<const N: usize, T, Y> Extend<NDPoint<N, T, Y>> for NDFeature<N, T, Y> {
    fn extend<I: IntoIterator<Item = NDPoint<N, T, Y>>>(&mut self, iter: I) {
        for pt in iter {
            self.push_raw(pt);
        }
    }
}

// -------- Specific supporting specializations: m/z + ion mobility ---------------

pub struct IMMZPeakIter<'a> {
    iter: NDIter<'a, 2, (MZ, IonMobility), Time>,
}

impl<'a> FusedIterator for IMMZPeakIter<'a> {}

impl<'a> ExactSizeIterator for IMMZPeakIter<'a> {
    fn len(&self) -> usize {
        <NDIter<'a, 2, (MZ, IonMobility), Time> as ExactSizeIterator>::len(&self.iter)
    }
}

impl<'a> Iterator for IMMZPeakIter<'a> {
    type Item = (IonMobilityAwareCentroidPeak, f64);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|pt| {
            (
                IonMobilityAwareCentroidPeak::new(pt.mz(), pt.ion_mobility(), pt.intensity(), 0),
                pt.time,
            )
        })
    }
}

impl<'a> IMMZPeakIter<'a> {
    pub fn new(iter: NDIter<'a, 2, (MZ, IonMobility), Time>) -> Self {
        Self { iter }
    }
}

impl AsPeakIter for NDFeature<2, (MZ, IonMobility), Time> {
    type Peak = IonMobilityAwareCentroidPeak;

    type Iter<'a>
        = IMMZPeakIter<'a>
    where
        Self: 'a;

    fn iter_peaks(&self) -> Self::Iter<'_> {
        IMMZPeakIter::new(self.iter())
    }
}

impl BuildFromPeak<IonMobilityAwareCentroidPeak> for NDFeature<2, (MZ, IonMobility), Time> {
    fn push_peak(&mut self, value: &IonMobilityAwareCentroidPeak, time: f64) {
        self.push_raw(NDPoint::new(
            [value.mz, value.ion_mobility],
            time,
            value.intensity,
        ))
    }
}

impl From<NDFeature<2, (MZ, IonMobility), Time>>
    for NDFeatureAdapter<MZ, (MZ, IonMobility), Time, NDFeature<2, (MZ, IonMobility), Time>>
{
    fn from(value: NDFeature<2, (MZ, IonMobility), Time>) -> Self {
        NDFeatureAdapter::new(value)
    }
}

// -------- Specific supporting specializations: mass + ion mobility ---------------

use super::charged::ChargedFeatureWrapper;

pub struct IMMassPeakIter<'a> {
    iter: NDIter<'a, 2, (Mass, IonMobility), Time>,
    charge: i32,
}

impl<'a> Iterator for IMMassPeakIter<'a> {
    type Item = (IonMobilityAwareDeconvolutedPeak, f64);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|pt| {
            (
                IonMobilityAwareDeconvolutedPeak::new(
                    pt.neutral_mass(),
                    pt.ion_mobility(),
                    self.charge,
                    pt.intensity(),
                    0,
                ),
                pt.time,
            )
        })
    }
}

impl<'a> IMMassPeakIter<'a> {
    pub fn new(iter: NDIter<'a, 2, (Mass, IonMobility), Time>, charge: i32) -> Self {
        Self { iter, charge }
    }
}

impl AsPeakIter
    for ChargedFeatureWrapper<(Mass, IonMobility), Time, NDFeature<2, (Mass, IonMobility), Time>>
{
    type Peak = IonMobilityAwareDeconvolutedPeak;

    type Iter<'a>
        = IMMassPeakIter<'a>
    where
        Self: 'a;

    fn iter_peaks(&self) -> Self::Iter<'_> {
        let (inner, charge) = self.as_inner();
        IMMassPeakIter::new(inner.iter(), charge)
    }
}

impl BuildFromPeak<IonMobilityAwareDeconvolutedPeak>
    for ChargedFeatureWrapper<(Mass, IonMobility), Time, NDFeature<2, (Mass, IonMobility), Time>>
{
    fn push_peak(&mut self, value: &IonMobilityAwareDeconvolutedPeak, time: f64) {
        self.push_raw(Charged::new(
            NDPoint::new(
                [value.neutral_mass(), value.ion_mobility()],
                time,
                value.intensity(),
            ),
            value.charge(),
        ));
    }
}


impl From<ChargedFeatureWrapper<(Mass, IonMobility), Time, NDFeature<2, (Mass, IonMobility), Time>>>
    for NDFeatureAdapter<Mass, (Mass, IonMobility), Time, ChargedFeatureWrapper<(Mass, IonMobility), Time, NDFeature<2, (Mass, IonMobility), Time>>>
{
    fn from(value: ChargedFeatureWrapper<(Mass, IonMobility), Time, NDFeature<2, (Mass, IonMobility), Time>>) -> Self {
        NDFeatureAdapter::new(value)
    }
}


/// An adapter from [`NDFeatureLike`] to [`FeatureLike`].
///
/// # Note
/// Assumes dimension `T0` is the 0th dimension in `T`, but this cannot be
/// enforced as Rust lacks variadic generics.
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NDFeatureAdapter<
    T0,
    T,
    Y,
    F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + PartialOrd, //+ CoordinateLike<T0>,
> {
    inner: F,
    _x: PhantomData<(T, Y, T0)>,
}

impl<
        T0,
        T,
        Y,
        F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<T0> + KnownCharge,
    > KnownCharge for NDFeatureAdapter<T0, T, Y, F>
{
    fn charge(&self) -> i32 {
        self.inner.charge()
    }
}

impl<
        T0,
        T,
        Y,
        F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<T0> + KnownChargeMut,
    > KnownChargeMut for NDFeatureAdapter<T0, T, Y, F>
{
    fn charge_mut(&mut self) -> &mut i32 {
        self.inner.charge_mut()
    }
}

impl<T0, T, Y, F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<T0>>
    PartialEq for NDFeatureAdapter<T0, T, Y, F>
{
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<T0, T, Y, F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<T0>>
    PartialOrd for NDFeatureAdapter<T0, T, Y, F>
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.inner.partial_cmp(&other.inner)
    }
}

impl<Q, T0, T, Y, F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<T0>>
    CoordinateLike<Q> for NDFeatureAdapter<T0, T, Y, F> where F: CoordinateLike<Q>
{
    fn coordinate(&self) -> f64 {
        <F as CoordinateLike<Q>>::coordinate(&self.inner)
    }
}

impl<T0, T, Y, F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<T0>>
    IntensityMeasurement for NDFeatureAdapter<T0, T, Y, F>
{
    fn intensity(&self) -> f32 {
        self.inner.intensity()
    }
}

impl<T0, T, Y, F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<T0>>
    FeatureLike<T0, Y> for NDFeatureAdapter<T0, T, Y, F>
{
    fn len(&self) -> usize {
        self.inner.len()
    }

    fn iter(&self) -> impl Iterator<Item = (f64, f64, f32)> {
        self.inner.iter().map(|pt| {
            (
                pt[0],
                <F::Point as CoordinateLike<Y>>::coordinate(&pt),
                pt.intensity(),
            )
        })
    }
}

impl<T0, T, Y, F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<T0>>
    TimeArray<Y> for NDFeatureAdapter<T0, T, Y, F>
{
    fn time_view(&self) -> &[f64] {
        <F as TimeArray<Y>>::time_view(&self.inner)
    }

    fn intensity_view(&self) -> &[f32] {
        <F as TimeArray<Y>>::intensity_view(&self.inner)
    }
}

impl<T0, T, Y, F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<T0>>
    TimeInterval<Y> for NDFeatureAdapter<T0, T, Y, F>
{
    fn start_time(&self) -> Option<f64> {
        <F as TimeInterval<Y>>::start_time(&self.inner)
    }

    fn end_time(&self) -> Option<f64> {
        <F as TimeInterval<Y>>::end_time(&self.inner)
    }

    fn apex_time(&self) -> Option<f64> {
        <F as TimeInterval<Y>>::apex_time(&self.inner)
    }

    fn area(&self) -> f32 {
        <F as TimeInterval<Y>>::area(&self.inner)
    }

    fn as_range(&self) -> crate::CoordinateRange<Y> {
        <F as TimeInterval<Y>>::as_range(&self.inner)
    }

    fn spans(&self, time: f64) -> bool {
        <F as TimeInterval<Y>>::spans(&self.inner, time)
    }

    fn iter_time(&self) -> impl Iterator<Item = f64> {
        <F as TimeInterval<Y>>::iter_time(&self.inner)
    }

    fn find_time(&self, time: f64) -> (Option<usize>, f64) {
        <F as TimeInterval<Y>>::find_time(&self.inner, time)
    }
}

impl<X, T, Y, F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<X>>
    NDFeatureAdapter<X, T, Y, F>
{
    pub fn new(inner: F) -> Self {
        Self {
            inner,
            _x: PhantomData,
        }
    }

    pub fn as_inner(&self) -> &F {
        &self.inner
    }

    pub fn into_inner(self) -> F {
        self.inner
    }

    pub fn as_mut(&mut self) -> &mut F {
        &mut self.inner
    }
}

impl<
        T0,
        T,
        Y,
        F: NDFeatureLike<T, Y> + TimeInterval<Y> + TimeArray<Y> + CoordinateLike<T0> + AsPeakIter,
    > AsPeakIter for NDFeatureAdapter<T0, T, Y, F>
{
    type Peak = F::Peak;

    type Iter<'a>
        = F::Iter<'a>
    where
        Self: 'a;

    fn iter_peaks(&self) -> Self::Iter<'_> {
        self.inner.iter_peaks()
    }
}

impl<
        P,
        T0,
        T,
        Y,
        F: NDFeatureLike<T, Y>
            + TimeInterval<Y>
            + TimeArray<Y>
            + CoordinateLike<T0>
            + BuildFromPeak<P>,
    > BuildFromPeak<P> for NDFeatureAdapter<T0, T, Y, F>
{
    fn push_peak(&mut self, value: &P, time: f64) {
        self.inner.push_peak(value, time);
    }
}

impl<const N: usize, T0, T, Y> FeatureLikeMut<T0, Y>
    for NDFeatureAdapter<T0, T, Y, NDFeature<N, T, Y>>
where
    NDFeature<N, T, Y>: CoordinateLike<T0>,
{
    fn iter_mut(&mut self) -> impl Iterator<Item = (&mut f64, &mut f64, &mut f32)> {
        super::feature::IterMut::<T0, Y>::from_parts(
            self.inner.dimensions[0].iter_mut(),
            self.inner.time.iter_mut(),
            self.inner.intensity.iter_mut(),
        )
    }

    fn push<P: CoordinateLike<T0> + IntensityMeasurement>(&mut self, pt: &P, time: f64) {
        let mut raw_pt = NDPoint::default();
        raw_pt.time = time;
        raw_pt.dimensions[0] = pt.coordinate();
        raw_pt.intensity = pt.intensity();
        self.inner.push_raw(raw_pt);
    }

    fn push_raw(&mut self, x: f64, y: f64, z: f32) {
        let mut raw_pt = NDPoint::default();
        raw_pt.time = y;
        raw_pt.dimensions[0] = x;
        raw_pt.intensity = z;
        self.inner.push_raw(raw_pt);
    }
}

impl<const N: usize, T0, T, Y> FeatureLikeMut<T0, Y>
    for NDFeatureAdapter<T0, T, Y, ChargedFeatureWrapper<T, Y, NDFeature<N, T, Y>>>
where
    ChargedFeatureWrapper<T, Y, NDFeature<N, T, Y>>: CoordinateLike<T0>,
{
    fn iter_mut(&mut self) -> impl Iterator<Item = (&mut f64, &mut f64, &mut f32)> {
        let (inner, _) = self.inner.as_mut();
        super::feature::IterMut::<T0, Y>::from_parts(
            inner.dimensions[0].iter_mut(),
            inner.time.iter_mut(),
            inner.intensity.iter_mut(),
        )
    }

    fn push<P: CoordinateLike<T0> + IntensityMeasurement>(&mut self, pt: &P, time: f64) {
        let mut raw_pt = NDPoint::default();
        raw_pt.time = time;
        raw_pt.dimensions[0] = pt.coordinate();
        raw_pt.intensity = pt.intensity();
        self.inner.as_mut().0.push_raw(raw_pt);
    }

    fn push_raw(&mut self, x: f64, y: f64, z: f32) {
        let mut raw_pt = NDPoint::default();
        raw_pt.time = y;
        raw_pt.dimensions[0] = x;
        raw_pt.intensity = z;
        self.inner.as_mut().0.push_raw(raw_pt);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        coordinate::{IonMobilityLocated, TimeLocated},
        prelude::*,
        IonMobility, MZLocated, Time, MZ,
    };

    #[test]
    fn test_create() {
        let pt = NDPoint::<2, (MZ, IonMobility), Time>::new([50.0, 0.5], 12.6, 1000.0);
        assert_eq!(pt.time(), 12.6);
        assert_eq!(pt.mz(), 50.0);
        assert_eq!(pt.ion_mobility(), 0.5);
    }

    fn test_creation_behavior(f: &NDFeature<2, (MZ, IonMobility), Time>) {
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

        assert_eq!(f.start_time(), Some(12.6));
        assert_eq!(f.end_time(), Some(12.9));
        assert_eq!(f.time_view().len(), 3);

        assert_eq!(f.iter().count(), 3);
        assert_eq!(f.intensity_view().iter().sum::<f32>(), f.intensity());
    }

    #[test]
    fn test_create_feature() {
        let points: Vec<NDPoint<2, (MZ, IonMobility), Time>> = vec![
            NDPoint::new([50.0, 0.5], 12.6, 1000.0),
            NDPoint::new([50.01, 0.49], 12.7, 1200.0),
            NDPoint::new([50.0, 0.51], 12.9, 800.0),
        ];

        let f: NDFeature<2, _, _> = points.clone().into_iter().collect();
        test_creation_behavior(&f);

        let mut f2 = NDFeature::default();
        for pt in points.clone() {
            f2.push_raw(pt);
        }
        test_creation_behavior(&f2);
        assert_eq!(f, f2);

        assert!(f <= f2);
        assert!(f >= f2);

        for (x, y) in f.iter().zip(f2.iter()) {
            assert_eq!(x, y);
            assert!(x <= y);
            assert!(x >= y);
        }

        f2 = NDFeature::default();
        for pt in points.clone().into_iter().rev() {
            f2.push_raw(pt);
        }
        test_creation_behavior(&f2);
        assert_eq!(f, f2);
    }
}
