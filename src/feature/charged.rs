use std::{cmp::Ordering, marker::PhantomData, ops::{Index, RangeBounds}};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{
    coordinate::{CoordinateLike, IonMobility, Mass, Time, MZ},
    CoordinateLikeMut, DeconvolutedPeak, IndexedCoordinate, IntensityMeasurement,
    IntensityMeasurementMut, KnownCharge, KnownChargeMut, MZLocated, MassLocated,
};

use super::traits::{FeatureLike, FeatureLikeMut, SplittableFeatureLike, TimeInterval};
use super::{
    feature::{Feature, FeatureView, Iter, IterMut},
    traits::BuildFromPeak,
    AsPeakIter, TimeArray,
};

/// A [`Feature`] with an associated `charge`, implementing the [`KnownCharge`] trait.
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
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

    fn clear(&mut self) {
        self.feature.clear();
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

    fn iter(&self) -> impl Iterator<Item = (f64, f64, f32)> {
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

    pub fn with_capacity(capacity: usize, charge: i32) -> Self {
        Self {
            feature: Feature::with_capacity(capacity),
            charge,
        }
    }

    pub fn empty(charge: i32) -> Self {
        Self::with_capacity(0, charge)
    }

    pub fn iter(&self) -> Iter<'_, X, Y> {
        self.feature.iter()
    }

    pub fn iter_mut(&mut self) -> IterMut<'_, X, Y> {
        self.feature.iter_mut()
    }

    pub fn as_view(&self) -> ChargedFeatureView<'_, X, Y> {
        ChargedFeatureView::new(self.feature.as_view(), self.charge)
    }

    pub fn into_inner(self) -> (Feature<X, Y>, i32) {
        (self.feature, self.charge)
    }

    pub fn as_inner(&self) -> (&Feature<X, Y>, i32) {
        (&self.feature, self.charge)
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

impl<Y> AsPeakIter for ChargedFeature<Mass, Y> {
    type Peak = DeconvolutedPeak;
    type Iter<'a>
        = DeconvolutedPeakIter<'a, Y>
    where
        Self: 'a;

    fn iter_peaks(&self) -> Self::Iter<'_> {
        self.iter_peaks()
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

impl<X, Y> KnownChargeMut for ChargedFeature<X, Y> {
    fn charge_mut(&mut self) -> &mut i32 {
        &mut self.charge
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

impl<'a, X, Y> TimeArray<Y> for ChargedFeature<X, Y> {
    fn time_view(&self) -> &[f64] {
        self.feature.time_view()
    }

    fn intensity_view(&self) -> &[f32] {
        self.feature.intensity_view()
    }
}

pub type DeconvolvedLCMSFeature = ChargedFeature<Mass, Time>;
pub type DeconvolvedIMSFeature = ChargedFeature<Mass, IonMobility>;

/// An iterator over a [`ChargedFeature`] which produces [`DeconvolutedPeak`]
/// instances with an associated time point.
pub struct DeconvolutedPeakIter<'a, Y> {
    charge: i32,
    point_iter: Iter<'a, Mass, Y>,
}

impl<'a, Y> DeconvolutedPeakIter<'a, Y> {
    pub fn new(source: &'a ChargedFeature<Mass, Y>) -> Self {
        Self {
            charge: source.charge(),
            point_iter: source.iter(),
        }
    }
}

impl<'a, Y> Iterator for DeconvolutedPeakIter<'a, Y> {
    type Item = (DeconvolutedPeak, f64);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((mass, time, intensity)) = self.point_iter.next() {
            Some((DeconvolutedPeak::new(mass, intensity, self.charge, 0), time))
        } else {
            None
        }
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        if let Some((mass, time, intensity)) = self.point_iter.nth(n) {
            Some((DeconvolutedPeak::new(mass, intensity, self.charge, 0), time))
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.point_iter.size_hint()
    }
}

impl<'a, Y> ExactSizeIterator for DeconvolutedPeakIter<'a, Y> {
    fn len(&self) -> usize {
        self.point_iter.len()
    }
}

impl<'a, Y> DoubleEndedIterator for DeconvolutedPeakIter<'a, Y> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if let Some((mass, time, intensity)) = self.point_iter.next_back() {
            Some((DeconvolutedPeak::new(mass, intensity, self.charge, 0), time))
        } else {
            None
        }
    }
}

impl<X, Y> IntoIterator for ChargedFeature<X, Y> {
    type Item = <Feature<X, Y> as IntoIterator>::Item;

    type IntoIter = <Feature<X, Y> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.feature.into_iter()
    }
}

impl<Y, T: MZLocated + IntensityMeasurement + KnownCharge> BuildFromPeak<T>
    for ChargedFeature<MZ, Y>
{
    fn push_peak(&mut self, value: T, time: f64) {
        self.push(&value, time);
    }
}

impl<Y, T: MassLocated + IntensityMeasurement + KnownCharge> BuildFromPeak<T>
    for ChargedFeature<Mass, Y>
{
    fn push_peak(&mut self, value: T, time: f64) {
        self.push(&value, time);
    }
}

/// A non-owning version of [`ChargedFeature`]
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
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

impl<'a, X, Y> KnownCharge for ChargedFeatureView<'a, X, Y> {
    fn charge(&self) -> i32 {
        self.charge
    }
}

impl<'a, X, Y> FeatureLike<X, Y> for ChargedFeatureView<'a, X, Y>
where
    ChargedFeatureView<'a, X, Y>: TimeInterval<Y>,
{
    fn len(&self) -> usize {
        self.feature.len()
    }

    fn iter(&self) -> impl Iterator<Item = (f64, f64, f32)> {
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

    pub fn empty(charge: i32) -> Self {
        Self::new(FeatureView::empty(), charge)
    }

    pub fn into_inner(self) -> (FeatureView<'a, X, Y>, i32) {
        (self.feature, self.charge)
    }

    pub fn as_inner(&self) -> (&FeatureView<'a, X, Y>, i32) {
        (&self.feature, self.charge)
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
        ChargedFeature::new(self.feature.to_owned(), self.charge)
    }
}

impl<'a, X, Y> SplittableFeatureLike<'a, X, Y> for ChargedFeature<X, Y> {
    type ViewType = ChargedFeatureView<'a, X, Y>;

    fn split_at_time(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        let (before, after) = self.feature.split_at_time(point);
        (
            Self::ViewType::new(before, self.charge),
            Self::ViewType::new(after, self.charge),
        )
    }

    fn slice<I: RangeBounds<usize> + Clone>(&'a self, bounds: I) -> Self::ViewType {
        let part = self.feature.slice(bounds);

        Self::ViewType::new(part, self.charge)
    }

    fn split_at(&'a self, index: usize) -> (Self::ViewType, Self::ViewType) {
        let (before, after) = self.feature.split_at(index);
        (
            Self::ViewType::new(before, self.charge),
            Self::ViewType::new(after, self.charge),
        )
    }
}

impl<'a, X, Y> SplittableFeatureLike<'a, X, Y> for ChargedFeatureView<'a, X, Y> {
    type ViewType = Self;

    fn split_at_time(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        let (before, after) = self.feature.split_at_time(point);
        (
            Self::ViewType::new(before, self.charge),
            Self::ViewType::new(after, self.charge),
        )
    }

    fn slice<I: RangeBounds<usize> + Clone>(&'a self, bounds: I) -> Self::ViewType {
        let part = self.feature.slice(bounds);
        Self::ViewType::new(part, self.charge)
    }

    fn split_at(&'a self, index: usize) -> (Self::ViewType, Self::ViewType) {
        let (before, after) = self.feature.split_at(index);
        (
            Self::ViewType::new(before, self.charge),
            Self::ViewType::new(after, self.charge),
        )
    }
}

impl<'a, X, Y> TimeArray<Y> for ChargedFeatureView<'a, X, Y> {
    fn time_view(&self) -> &[f64] {
        self.feature.time_view()
    }

    fn intensity_view(&self) -> &[f32] {
        self.feature.intensity_view()
    }
}

impl<'a, X, Y> IntoIterator for ChargedFeatureView<'a, X, Y> {
    type Item = (f64, f64, f32);

    type IntoIter = Iter<'a, X, Y>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ChargedFeatureWrapper<X, Y, F: TimeInterval<Y> + TimeArray<Y>> {
    inner: F,
    charge: i32,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<X, Y, F: TimeInterval<Y> + TimeArray<Y>> TimeArray<Y> for ChargedFeatureWrapper<X, Y, F> {
    fn time_view(&self) -> &[f64] {
        <F as TimeArray<Y>>::time_view(&self.inner)
    }

    fn intensity_view(&self) -> &[f32] {
        <F as TimeArray<Y>>::intensity_view(&self.inner)
    }
}

impl<X, Y, F: TimeInterval<Y> + TimeArray<Y>> ChargedFeatureWrapper<X, Y, F> {
    pub fn new(inner: F, charge: i32) -> Self {
        Self {
            inner,
            charge,
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    pub fn into_inner(self) -> (F, i32) {
        (self.inner, self.charge)
    }

    pub fn as_inner(&self) -> (&F, i32) {
        (&self.inner, self.charge)
    }

    pub fn as_mut(&mut self) -> (&mut F, &mut i32) {
        (&mut self.inner, &mut self.charge)
    }
}

impl<'a, X, Y, F: FeatureLike<X, Y> + TimeInterval<Y> + TimeArray<Y>>
    SplittableFeatureLike<'a, X, Y> for ChargedFeatureWrapper<X, Y, F>
where
    F: SplittableFeatureLike<'a, X, Y>,
    F::ViewType: FeatureLike<X, Y> + TimeInterval<Y> + TimeArray<Y>,
{
    type ViewType = ChargedFeatureWrapper<X, Y, F::ViewType>;

    fn split_at(&'a self, index: usize) -> (Self::ViewType, Self::ViewType) {
        let (a, b) = self.inner.split_at(index);
        (
            ChargedFeatureWrapper::new(a, self.charge),
            ChargedFeatureWrapper::new(b, self.charge),
        )
    }

    fn split_at_time(&'a self, point: f64) -> (Self::ViewType, Self::ViewType) {
        let (a, b) = self.inner.split_at_time(point);
        (
            ChargedFeatureWrapper::new(a, self.charge),
            ChargedFeatureWrapper::new(b, self.charge),
        )
    }

    fn slice<I: RangeBounds<usize> + Clone>(&'a self, bounds: I) -> Self::ViewType {
        let a = self.inner.slice(bounds);
        ChargedFeatureWrapper::new(a, self.charge)
    }
}

impl<X, Y, F: TimeInterval<Y> + TimeArray<Y>> KnownChargeMut for ChargedFeatureWrapper<X, Y, F> {
    fn charge_mut(&mut self) -> &mut i32 {
        &mut self.charge
    }
}

impl<X, Y, F: TimeInterval<Y> + TimeArray<Y>> KnownCharge for ChargedFeatureWrapper<X, Y, F> {
    fn charge(&self) -> i32 {
        self.charge
    }
}

impl<X, Y, F: IntensityMeasurement + TimeInterval<Y> + TimeArray<Y>> IntensityMeasurement
    for ChargedFeatureWrapper<X, Y, F>
{
    fn intensity(&self) -> f32 {
        self.inner.intensity()
    }
}

impl<X, Y, F: TimeInterval<Y> + TimeArray<Y>> TimeInterval<Y> for ChargedFeatureWrapper<X, Y, F> {
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

impl<X, Y, F: PartialEq<F> + TimeInterval<Y> + TimeArray<Y>> PartialEq
    for ChargedFeatureWrapper<X, Y, F>
{
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
            && self.charge == other.charge
            && self._x == other._x
            && self._y == other._y
    }
}

impl<X, Y, F: PartialOrd<F> + TimeInterval<Y> + TimeArray<Y>> PartialOrd
    for ChargedFeatureWrapper<X, Y, F>
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match self.inner.partial_cmp(&other.inner) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self.charge.partial_cmp(&other.charge) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self._x.partial_cmp(&other._x) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        self._y.partial_cmp(&other._y)
    }
}

impl<X, Y, F: CoordinateLike<X> + TimeInterval<Y> + TimeArray<Y>> CoordinateLike<X>
    for ChargedFeatureWrapper<X, Y, F>
{
    fn coordinate(&self) -> f64 {
        self.inner.coordinate()
    }
}

impl<X, Y, F: FeatureLike<X, Y> + TimeInterval<Y> + TimeArray<Y>> FeatureLike<X, Y>
    for ChargedFeatureWrapper<X, Y, F>
{
    fn len(&self) -> usize {
        <F as FeatureLike<X, Y>>::len(&self.inner)
    }

    fn iter(&self) -> impl Iterator<Item = (f64, f64, f32)> {
        <F as FeatureLike<X, Y>>::iter(&self.inner)
    }

    fn is_empty(&self) -> bool {
        <F as FeatureLike<X, Y>>::is_empty(&self.inner)
    }

    fn at(&self, index: usize) -> Option<(f64, f64, f32)> {
        <F as FeatureLike<X, Y>>::at(&self.inner, index)
    }

    fn first(&self) -> Option<(f64, f64, f32)> {
        <F as FeatureLike<X, Y>>::first(&self.inner)
    }

    fn last(&self) -> Option<(f64, f64, f32)> {
        <F as FeatureLike<X, Y>>::last(&self.inner)
    }

    fn at_time(&self, time: f64) -> Option<(f64, f64, f32)> {
        <F as FeatureLike<X, Y>>::at_time(&self.inner, time)
    }
}

impl<X, Y, F: FeatureLike<X, Y> + TimeInterval<Y> + TimeArray<Y>> FeatureLikeMut<X, Y>
    for ChargedFeatureWrapper<X, Y, F>
where
    F: FeatureLikeMut<X, Y>,
{
    fn iter_mut(&mut self) -> impl Iterator<Item = (&mut f64, &mut f64, &mut f32)> {
        <F as FeatureLikeMut<X, Y>>::iter_mut(&mut self.inner)
    }

    fn push<T: CoordinateLike<X> + IntensityMeasurement>(&mut self, pt: &T, time: f64) {
        <F as FeatureLikeMut<X, Y>>::push(&mut self.inner, pt, time)
    }

    fn push_raw(&mut self, x: f64, y: f64, z: f32) {
        <F as FeatureLikeMut<X, Y>>::push_raw(&mut self.inner, x, y, z)
    }

    fn clear(&mut self) {
        <F as FeatureLikeMut<X, Y>>::clear(&mut self.inner)
    }

    fn at_mut(&mut self, index: usize) -> Option<(&mut f64, f64, &mut f32)> {
        <F as FeatureLikeMut<X, Y>>::at_mut(&mut self.inner, index)
    }

    fn first_mut(&mut self) -> Option<(&mut f64, f64, &mut f32)> {
        <F as FeatureLikeMut<X, Y>>::first_mut(&mut self.inner)
    }

    fn last_mut(&mut self) -> Option<(&mut f64, f64, &mut f32)> {
        <F as FeatureLikeMut<X, Y>>::last_mut(&mut self.inner)
    }

    fn at_time_mut(&mut self, time: f64) -> Option<(&mut f64, f64, &mut f32)> {
        <F as FeatureLikeMut<X, Y>>::at_time_mut(&mut self.inner, time)
    }
}

use super::{NDFeatureLike, NDFeatureLikeMut};

#[derive(Debug, Default, Clone, PartialEq, PartialOrd)]
pub struct Charged<T>(pub T, pub i32);

impl<T> AsRef<T> for Charged<T> {
    fn as_ref(&self) -> &T {
        &self.0
    }
}

impl<T> Charged<T> {
    pub fn new(point: T, charge: i32) -> Self {
        Self(point, charge)
    }

    pub fn into(self) -> T {
        self.0
    }
}

impl<D, T: CoordinateLike<D>> CoordinateLike<D> for Charged<T> {
    fn coordinate(&self) -> f64 {
        self.0.coordinate()
    }
}

impl<T> KnownCharge for Charged<T> {
    fn charge(&self) -> i32 {
        self.1
    }
}

impl<T> KnownChargeMut for Charged<T> {
    fn charge_mut(&mut self) -> &mut i32 {
        &mut self.1
    }
}

impl<D, T: CoordinateLikeMut<D>> CoordinateLikeMut<D> for Charged<T> {
    fn coordinate_mut(&mut self) -> &mut f64 {
        self.0.coordinate_mut()
    }
}

impl<T> IntensityMeasurement for Charged<T>
where
    T: IntensityMeasurement,
{
    fn intensity(&self) -> f32 {
        self.0.intensity()
    }
}

impl<T> IntensityMeasurementMut for Charged<T>
where
    T: IntensityMeasurementMut,
{
    fn intensity_mut(&mut self) -> &mut f32 {
        self.0.intensity_mut()
    }
}

impl<D, T: IndexedCoordinate<D>> IndexedCoordinate<D> for Charged<T> {
    fn get_index(&self) -> crate::IndexType {
        self.0.get_index()
    }

    fn set_index(&mut self, index: crate::IndexType) {
        self.0.set_index(index);
    }
}

impl<X, Y, F: NDFeatureLike<X, Y> + TimeInterval<Y> + TimeArray<Y>> NDFeatureLike<X, Y>
    for ChargedFeatureWrapper<X, Y, F> where Charged<F::Point>: Index<usize, Output = f64>
{
    type Point = Charged<F::Point>;

    fn len(&self) -> usize {
        self.inner.len()
    }

    fn iter(&self) -> impl Iterator<Item = Self::Point> {
        self.inner.iter().map(|pt| Charged::new(pt, self.charge))
    }

    fn coordinate(&self) -> Self::Point {
        let pt = self.inner.coordinate();
        Charged::new(pt, self.charge)
    }
}

impl<X, Y, F: NDFeatureLikeMut<X, Y> + TimeInterval<Y> + TimeArray<Y>>
    NDFeatureLikeMut<X, Y> for ChargedFeatureWrapper<X, Y, F> where Charged<F::Point>: Index<usize, Output = f64>
{
    type PointMutRef<'a>
        = Charged<F::PointMutRef<'a>>
    where
        Self: 'a;

    fn iter_mut(&mut self) -> impl Iterator<Item = Self::PointMutRef<'_>> {
        self.inner
            .iter_mut()
            .map(|pt| Charged::new(pt, self.charge))
    }

    fn push<P: Into<Self::Point>>(&mut self, pt: P, time: f64) {
        let point = pt.into();
        self.inner.push(point.0, time);
    }

    fn push_raw(&mut self, point: Self::Point) {
        self.inner.push_raw(point.0);
    }
}
