//! Collections of features that are ordered and searchable by a coordinate,
//! support fast access and are growable. While the main behaviors are provided
//! through the [`FeatureMapLike`] generic trait, a (generic) full implementation
//! is given by [`FeatureMap`].
//!

use crate::{feature::FeatureLike, CoordinateLike, Tolerance};
use std::{marker::PhantomData, ops::{self, Range}};


/// A two dimensional feature collection where features are sorted by the `X` dimension
/// and each feature is internally sorted by the `Y` dimension.
pub trait FeatureMapLike<X, Y, T: FeatureLike<X, Y>>: ops::Index<usize, Output = T>
{
    fn search_by(&self, query: f64) -> Result<usize, usize>;
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Implement index access
    fn get_item(&self, i: usize) -> &T;
    fn get_slice(&self, i: ops::Range<usize>) -> &[T];

    #[inline]
    fn _closest_feature(&self, query: f64, error_tolerance: Tolerance, i: usize) -> Option<usize> {
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
    /// this feature collection, or `None`.
    fn search(&self, query: f64, error_tolerance: Tolerance) -> Option<usize> {
        let lower_bound = error_tolerance.bounds(query).0;

        match self.search_by(lower_bound) {
            Ok(j) => self._closest_feature(query, error_tolerance, j),
            Err(j) => self._closest_feature(query, error_tolerance, j),
        }
    }

    #[inline]
    /// Return the feature nearest to `query` within `error_tolerance` in
    /// this feature collection, or `None`.
    fn has_feature(&self, query: f64, error_tolerance: Tolerance) -> Option<&T> {
        return match self.search(query, error_tolerance) {
            Some(j) => Some(self.get_item(j)),
            None => None,
        };
    }

    #[inline]
    /// Return the index range containing all features between `low` and `high` coordinates within
    /// `error_tolerance`.
    ///
    /// See [`FeatureMapLike::between`] for that retrieves these features.
    fn indices_between(&self, low: f64, high: f64, error_tolerance: Tolerance) -> Range<usize> {
        let lower_bound = error_tolerance.bounds(low).0;
        let upper_bound = error_tolerance.bounds(high).1;

        let n = self.len();
        if n == 0 {
            return 0..0;
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
            return 0..0;
        }

        lower_index..upper_index
    }

    #[inline]
    /// Return a slice containing all features between `low` and `high` coordinates within
    /// `error_tolerance`.
    fn between(&self, low: f64, high: f64, error_tolerance: Tolerance) -> &[T] {
        let indices = self.indices_between(low, high, error_tolerance);
        let subset = self.get_slice(indices);
        subset
    }

    #[inline]
    /// Find all feature indices which could match `query` within `error_tolerance` units.
    ///
    /// See [`FeatureMapLike::all_features_for`] for the function that retrieves these features.
    fn all_indices_for(&self, query: f64, error_tolerance: Tolerance) -> Range<usize> {
        let (lower_bound, upper_bound) = error_tolerance.bounds(query);

        let n = self.len();
        if n == 0 {
            return 0..0
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

        lower_index..upper_index + 1
    }

    #[inline]
    /// Find all features which could match `query` within `error_tolerance` units
    fn all_features_for(&self, query: f64, error_tolerance: Tolerance) -> &[T] {
        let c = self.all_indices_for(query, error_tolerance);
        return self.get_slice(c);
    }
}

/// A mutable kind of [`FeatureMapLike`] which new features can be added to.
pub trait FeatureMapLikeMut<X, Y, T: FeatureLike<X, Y>> : FeatureMapLike<X, Y, T>
{

    /// Add `feature` to the collection, maintaining sort order and feature
    /// indexing.
    fn push(&mut self, feature: T);

    /// Sort the collection, updating the feature indexing.
    fn sort(&mut self);
}


/// Represents a sorted list of mass spectral features that is a concrete implementation
/// of [`FeatureMapLike`] and [`FeatureMapLikeMut`]
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FeatureMap<X, Y, T: FeatureLike<X, Y>> {
    features: Vec<T>,
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<'a, X, Y, T: FeatureLike<X, Y>> FeatureMap<X, Y, T> {
    /// Create a new [`FeatureMap`] from an existing `Vec<T>` and sorts
    /// the newly created structure to ensure it is ordered by coordinate `X`
    pub fn new(features: Vec<T>) -> Self {
        let mut inst = Self::wrap(features);
        inst.sort();
        inst
    }

    pub fn empty() -> Self {
        Self {
            features: Vec::new(),
            _x: PhantomData,
            _y: PhantomData
        }
    }

    /// Create a new [`FeatureMap`] from an existing `Vec<T>`, but does not actively
    /// sort the collection. It is up to the caller to ensure that the provided `Vec`
    /// is sorted or that it will be sorted prior to any of its search functionality
    /// is used.
    pub fn wrap(features: Vec<T>) -> Self {
        Self {
            features,
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    pub fn len(&self) -> usize {
        self.features.len()
    }

    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }

    /// Iterate over references to features
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.features.iter()
    }

    /// Iterate over mutable reference to features
    pub fn iter_mut(&'a mut self) -> std::slice::IterMut<'a, T> {
        self.features.iter_mut()
    }

    pub fn search_by(&self, query: f64) -> Result<usize, usize> {
        self.features
            .binary_search_by(|feature| feature.coordinate().partial_cmp(&query).unwrap())
    }

    /// Extract a subset of this [`FeatureMap`] that overlap the specified `y` coordinate
    pub fn spanning(&'a self, y: f64) -> FeatureMap<X, Y, &T> {
        let subset: Vec<_> = self.iter().filter(|f| f.spans(y)).collect();
        FeatureMap::wrap(subset)
    }

    pub fn from_iter<I: Iterator<Item=T>>(iter: I, sort: bool) -> Self {
        let features = iter.collect();
        if sort {
            Self::new(features)
        } else {
            Self::wrap(features)
        }
    }

    pub fn earliest_time(&self) -> Option<f64> {
        self.features.iter().fold(Option::<f64>::None, |prev, feat| {
            match (prev, feat.start_time()) {
                (Some(prev), Some(cur)) => {
                    Some(prev.min(cur))
                },
                (None, Some(cur)) => {
                    Some(cur)
                },
                (_, _) => prev
            }
        })
    }

    pub fn latest_time(&self) -> Option<f64> {
        self.features.iter().fold(Option::<f64>::None, |prev, feat| {
            match (prev, feat.start_time()) {
                (Some(prev), Some(cur)) => {
                    Some(prev.max(cur))
                },
                (None, Some(cur)) => {
                    Some(cur)
                },
                (_, _) => prev
            }
        })
    }
}

impl<X, Y, T: FeatureLike<X, Y>> FeatureMapLike<X, Y, T> for FeatureMap<X, Y, T> {
    fn search_by(&self, query: f64) -> Result<usize, usize> {
        self.search_by(query)
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn get_item(&self, i: usize) -> &T {
        &self.features[i]
    }

    fn get_slice(&self, i: ops::Range<usize>) -> &[T] {
        &self.features[i]
    }
}

impl<X, Y, T: FeatureLike<X, Y>> FeatureMapLikeMut<X, Y, T> for FeatureMap<X, Y, T> {
    fn push(&mut self, feature: T) {
        if self.is_empty() {
            self.features.push(feature)
        } else {
            let is_tail = self.features.last().as_ref().unwrap().coordinate() <= feature.coordinate();
            self.features.push(feature);
            if !is_tail {
                self.sort();
            }
        }
    }

    fn sort(&mut self) {
        self.features.sort_by(|x, y| {
            x.partial_cmp(y).unwrap()
        })
    }
}

impl<X, Y, T: FeatureLike<X, Y>> FromIterator<T> for FeatureMap<X, Y, T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let items = iter.into_iter().collect();
        Self::new(items)
    }
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

impl<X, Y, T: FeatureLike<X, Y>> ops::Index<usize> for FeatureMap<X, Y, T> {
    type Output = T;

    fn index(&self, i: usize) -> &Self::Output {
        &(self.features[i])
    }
}

impl<X, Y, T: FeatureLike<X, Y>> ops::IndexMut<usize> for FeatureMap<X, Y, T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.features[index]
    }
}

impl_slicing!(FeatureMap<X, Y, T>, X, Y, T: FeatureLike<X, Y>);


#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct FeatureMapView<'a, X, Y, T: FeatureLike<X, Y>> {
    features: &'a [T],
    _x: PhantomData<X>,
    _y: PhantomData<Y>,
}

impl<'a, X, Y, T: FeatureLike<X, Y>> FeatureMapView<'a, X, Y, T> {
    pub fn new(features: &'a [T]) -> Self {
        Self {
            features,
            _x: PhantomData,
            _y: PhantomData,
        }
    }

    pub fn len(&self) -> usize {
        self.features.len()
    }

    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.features.iter()
    }

    fn search_by(&self, query: f64) -> Result<usize, usize> {
        self.features
            .binary_search_by(|feature| feature.coordinate().partial_cmp(&query).unwrap())
    }
}

impl<'a, X, Y, T: FeatureLike<X, Y>> ops::Index<usize> for FeatureMapView<'a, X, Y, T> {
    type Output = T;

    fn index(&self, i: usize) -> &Self::Output {
        &(self.features[i])
    }
}

impl_slicing!(FeatureMapView<'a, X, Y, T>, 'a, X, Y, T: FeatureLike<X, Y>);

impl<'a, X, Y, T: FeatureLike<X, Y>> FeatureMapLike<X, Y, T> for FeatureMapView<'a, X, Y, T> {
    fn search_by(&self, query: f64) -> Result<usize, usize> {
        self.search_by(query)
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn get_item(&self, i: usize) -> &T {
        &self.features[i]
    }

    fn get_slice(&self, i: ops::Range<usize>) -> &[T] {
        &self.features[i]
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::feature::LCMSFeature;
    use crate::test_data;
    use crate::prelude::*;

    #[test]
    fn test_sequence_behavior() {
        let source_peaks = test_data::read_peaks_from_file("./test/data/test.txt").unwrap();

        let feature_map: FeatureMap<_, _, LCMSFeature> = source_peaks.iter().map(|p| {
            let mut feature = LCMSFeature::empty();
            feature.push(&p, 1.0);
            feature
        }).collect();

        assert_eq!(feature_map.len(), 485);
        assert!((feature_map[0].mz() - 231.3888).abs() < 1e-3);

        for (i, feature) in feature_map.iter().enumerate() {
            if i > 0 {
                assert!(feature.mz() > feature_map[i - 1].mz());
            }
        }

        let part = feature_map.search(773.4414, Tolerance::Da(0.01));
        assert_eq!(part.expect("Match feature"), 300);

        let part = feature_map.all_features_for(773.4414, Tolerance::PPM(10.0));
        assert_eq!(part.len(), 1);
        // assert_eq!(part[0].index, 300);

        let part = feature_map.search(773.4414, Tolerance::Da(50.0));
        assert_eq!(part.expect("Match feature"), 300);

        let part = feature_map.all_features_for(736.637, Tolerance::PPM(10.0));
        assert_eq!(part.len(), 1);
        let part = feature_map.all_features_for(736.237, Tolerance::PPM(10.0));
        assert_eq!(part.len(), 0);

        let q = 1221.639893;
        let block = feature_map.all_features_for(q, Tolerance::Da(0.5));
        assert_eq!(block.len(), 1);

        let q = 2000.0;
        let block = feature_map.all_features_for(q, Tolerance::PPM(10.));
        assert_eq!(block.len(), 0);

        let q = -2000.0;
        let block = feature_map.all_features_for(q, Tolerance::PPM(10.));
        assert_eq!(block.len(), 0);

        let block = feature_map.between(-2000f64, 2000f64, Tolerance::PPM(10.0));
        assert_eq!(block.len(), feature_map.len());

        let block = feature_map.between(0.0, 2000f64, Tolerance::PPM(10.0));
        assert_eq!(block.len(), feature_map.len());

        let block = feature_map.between(1313.0, 1316.0, "10.0ppm".parse().unwrap());
        assert_eq!(block.len(), 3);

    }

    #[test]
    fn test_edgecases() {
        let features = FeatureMap::new(vec![
            LCMSFeature::from_iter([(500.0, 2., 1.0)])
        ]);

        let p = features.has_feature(500.0, Tolerance::Da(1.0));
        assert!(p.is_some());

        let p = features.all_features_for(500.0, Tolerance::Da(1.0));
        assert!(p.len() == 1);

        let features: FeatureMap<_, _, LCMSFeature> = FeatureMap::new(vec![]);

        let p = features.has_feature(500.0, Tolerance::Da(1.0));
        assert!(p.is_none());

        let p = features.all_features_for(500.0, Tolerance::Da(1.0));
        assert!(p.len() == 0);
    }
}