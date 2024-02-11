use core::slice;
use std::{cmp::Ordering, marker::PhantomData};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{
    coordinate::{CoordinateLike, IonMobility, Mass, Time, MZ}, CentroidPeak, CoordinateRange, DeconvolutedPeak, IntensityMeasurement, KnownCharge, MassLocated
};

pub trait TimeInterval<T> {
    fn start_time(&self) -> Option<f64>;
    fn end_time(&self) -> Option<f64>;
    fn apex_time(&self) -> Option<f64>;
    fn area(&self) -> f64;

    fn as_range(&self) -> CoordinateRange<T> {
        CoordinateRange::new(self.start_time(), self.end_time())
    }
}

#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Feature<X, Y> {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f32>,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

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

    fn empty() -> Self {
        Self {
            x: Vec::new(),
            y: Vec::new(),
            z: Vec::new(),
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    fn weighted_average(&self, x: &[f64]) -> f64 {
        let (acc, norm) = x
            .iter()
            .zip(self.z.iter())
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

    fn coordinate_x(&self) -> f64 {
        self.weighted_average(&self.x)
    }

    fn coordinate_y(&self) -> f64 {
        self.weighted_average(&self.y)
    }

    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn transpose(self) -> Feature<Y, X> {
        Feature::<Y, X> {
            x: self.y,
            y: self.x,
            z: self.z,
            _x: self._y,
            _y: self._x,
        }
    }

    fn idxmax(&self) -> Option<usize> {
        let pt = self
            .z
            .iter()
            .enumerate()
            .reduce(|(best_i, best), (current_i, current)| {
                if *current > *best {
                    (current_i, current)
                } else {
                    (best_i, best)
                }
            });

        pt.and_then(|(i, _)| Some(i))
    }

    fn apex_y(&self) -> Option<f64> {
        self.idxmax().and_then(|i| self.y.get(i).copied())
    }

    pub fn push<T: CoordinateLike<X> + IntensityMeasurement>(&mut self, pt: &T, time: f64) {
        let x = pt.coordinate();
        let z = pt.intensity();
        self.x.push(x);
        self.y.push(time);
        self.z.push(z);
    }

    pub fn push_raw(&mut self, x: f64, y: f64, z: f32) {
        self.x.push(x);
        self.y.push(y);
        self.z.push(z);
    }

    pub fn iter(&self) -> Iter<'_, X, Y> {
        Iter::new(self)
    }

    pub fn iter_mut(&mut self) -> IterMut<'_, X, Y> {
        IterMut::new(self)
    }

    fn integrate_y(&self) -> f64 {
        let mut it = self.y.iter().zip(self.z.iter());
        if let Some((first_y, first_z)) = it.next() {
            let (_y, _z, acc) =
                it.fold((first_y, first_z, 0.0), |(last_y, last_z, acc), (y, z)| {
                    let step = (last_z + z) / 2.0;
                    let dy = y - last_y;
                    (y, z, acc + (step as f64 * dy))
                });
            acc
        } else {
            0.0
        }
    }

    pub fn total_intensity(&self) -> f32 {
        self.z.iter().sum()
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
        Some(self.coordinate_y().total_cmp(&other.coordinate_y()))
    }
}

impl<Y> CoordinateLike<MZ> for Feature<MZ, Y> {
    fn coordinate(&self) -> f64 {
        self.coordinate_x()
    }
}

impl<Y> CoordinateLike<Mass> for Feature<Mass, Y> {
    fn coordinate(&self) -> f64 {
        self.coordinate_x()
    }
}

impl<X, Y> IntensityMeasurement for Feature<X, Y> {
    fn intensity(&self) -> f32 {
        self.total_intensity()
    }
}

impl<X> TimeInterval<Time> for Feature<X, Time> {
    fn apex_time(&self) -> Option<f64> {
        self.apex_y()
    }

    fn area(&self) -> f64 {
        self.integrate_y()
    }

    fn end_time(&self) -> Option<f64> {
        self.y.last().copied()
    }

    fn start_time(&self) -> Option<f64> {
        self.y.first().copied()
    }
}

impl<X> TimeInterval<IonMobility> for Feature<X, IonMobility> {
    fn apex_time(&self) -> Option<f64> {
        self.apex_y()
    }

    fn area(&self) -> f64 {
        self.integrate_y()
    }

    fn start_time(&self) -> Option<f64> {
        self.y.first().copied()
    }

    fn end_time(&self) -> Option<f64> {
        self.y.last().copied()
    }
}

impl<Y> Feature<MZ, Y> {
    pub fn iter_peaks(&self) -> MZPeakIter<'_, Y> {
        MZPeakIter::new(self)
    }
}

pub type MZLCMSFeature = Feature<MZ, Time>;
pub type MZIMSFeature = Feature<MZ, IonMobility>;

pub struct Iter<'a, X, Y> {
    source: &'a Feature<X, Y>,
    xiter: slice::Iter<'a, f64>,
    yiter: slice::Iter<'a, f64>,
    ziter: slice::Iter<'a, f32>,
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
        self.source.len()
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
            source,
            xiter: source.x.iter(),
            yiter: source.y.iter(),
            ziter: source.z.iter(),
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

    pub fn len(&self) -> usize {
        self.feature.len()
    }
}

impl<Y> ChargedFeature<Mass, Y> {
    pub fn iter_peaks(&self) -> DeconvolutedPeakIter<'_, Y> {
        DeconvolutedPeakIter::new(self)
    }
}

impl<X> TimeInterval<Time> for ChargedFeature<X, Time> {
    fn apex_time(&self) -> Option<f64> {
        self.feature.apex_time()
    }

    fn area(&self) -> f64 {
        self.feature.area()
    }

    fn start_time(&self) -> Option<f64> {
        self.feature.start_time()
    }

    fn end_time(&self) -> Option<f64> {
        self.feature.end_time()
    }
}

impl<X> TimeInterval<IonMobility> for ChargedFeature<X, IonMobility> {
    fn apex_time(&self) -> Option<f64> {
        <Feature<X, IonMobility> as TimeInterval<IonMobility>>::apex_time(&self.feature)
    }

    fn area(&self) -> f64 {
        <Feature<X, IonMobility> as TimeInterval<IonMobility>>::area(&self.feature)
    }

    fn start_time(&self) -> Option<f64> {
        <Feature<X, IonMobility> as TimeInterval<IonMobility>>::start_time(&self.feature)
    }

    fn end_time(&self) -> Option<f64> {
        <Feature<X, IonMobility> as TimeInterval<IonMobility>>::end_time(&self.feature)
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

impl<Y> CoordinateLike<MZ> for ChargedFeature<MZ, Y> {
    fn coordinate(&self) -> f64 {
        <Feature<MZ, Y> as CoordinateLike<MZ>>::coordinate(&self.feature)
    }
}

impl<Y> CoordinateLike<MZ> for ChargedFeature<Mass, Y> {
    fn coordinate(&self) -> f64 {
        let charge_carrier: f64 = 1.007276;
        let charge = self.charge as f64;
        (self.neutral_mass() + charge_carrier * charge) / charge
    }
}

impl<Y> CoordinateLike<Mass> for ChargedFeature<Mass, Y> {
    fn coordinate(&self) -> f64 {
        <Feature<Mass, Y> as CoordinateLike<Mass>>::coordinate(&self.feature)
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

#[cfg(test)]
mod test {
    use super::*;
    use crate::{CentroidPeak, DeconvolutedPeak, MZLocated};

    #[test]
    fn test_build_raw() {
        let mut x = MZLCMSFeature::empty();

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
