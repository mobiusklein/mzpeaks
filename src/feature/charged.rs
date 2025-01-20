use std::{cmp::Ordering, ops::RangeBounds};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{
    coordinate::{CoordinateLike, IonMobility, Mass, Time, MZ}, DeconvolutedPeak, IntensityMeasurement, KnownCharge, MZLocated, MassLocated
};

use super::{feature::{Feature, FeatureView, Iter, IterMut}, traits::BuildFromPeak, PeakSeries, TimeArray};
use super::traits::{FeatureLike, FeatureLikeMut, SplittableFeatureLike, TimeInterval};

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

impl<'a, Y: 'a> PeakSeries<'a> for ChargedFeature<Mass, Y> {
    type Peak = DeconvolutedPeak;

    type Iter = DeconvolutedPeakIter<'a, Y>;

    fn iter_peaks(&'a self) -> Self::Iter {
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
                DeconvolutedPeak::new(mass, intensity, self.source.charge, 0),
                time,
            ))
        } else {
            None
        }
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        if let Some((mass, time, intensity)) = self.point_iter.nth(n) {
            Some((
                DeconvolutedPeak::new(mass, intensity, self.source.charge, 0),
                time,
            ))
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
        self.source.len()
    }
}

impl<'a, Y> DoubleEndedIterator for DeconvolutedPeakIter<'a, Y> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if let Some((mass, time, intensity)) = self.point_iter.next_back() {
            Some((
                DeconvolutedPeak::new(mass, intensity, self.source.charge, 0),
                time,
            ))
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


impl<Y, T: MZLocated + IntensityMeasurement + KnownCharge> BuildFromPeak<T> for ChargedFeature<MZ, Y> {
    fn push_peak(&mut self, value: T, time: f64) {
        self.push(&value, time);
    }
}

impl<Y, T: MassLocated + IntensityMeasurement + KnownCharge> BuildFromPeak<T> for ChargedFeature<Mass, Y> {
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