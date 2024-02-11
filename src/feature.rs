use core::slice;
use std::{cmp::Ordering, marker::PhantomData};

use crate::{
    coordinate::{CoordinateLike, IonMobility, Mass, Time, MZ},
    IntensityMeasurement,
};

#[derive(Debug, Default, Clone)]
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

impl<X> Feature<X, Time> {
    pub fn apex_time(&self) -> Option<f64> {
        self.apex_y()
    }

    pub fn area(&self) -> f64 {
        self.integrate_y()
    }

    pub fn start_time(&self) -> Option<f64> {
        self.y.first().copied()
    }

    pub fn end_time(&self) -> Option<f64> {
        self.y.last().copied()
    }
}

impl<X> Feature<X, IonMobility> {
    pub fn apex_time(&self) -> Option<f64> {
        self.apex_y()
    }

    pub fn area(&self) -> f64 {
        self.integrate_y()
    }

    pub fn start_time(&self) -> Option<f64> {
        self.y.first().copied()
    }

    pub fn end_time(&self) -> Option<f64> {
        self.y.last().copied()
    }
}

pub type MZLCMSFeature = Feature<MZ, Time>;
pub type MassLCMSFeature = Feature<Mass, Time>;

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


#[cfg(test)]
mod test {
    use super::*;
    use crate::{CentroidPeak, MZLocated};

    #[test]
    fn test_build() {
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

        assert_eq!(x, points.into_iter().collect());
    }
}