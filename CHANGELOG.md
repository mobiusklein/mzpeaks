# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.8] - 2025-03-08

### Added

- Add AVX implementation of weighted average
- Add `FeatureLikeMut::reserve` to allow memory allocation
- Add `FeatureMap::into_inner` to get access to the underlying memory allocation

### Fixed

- Fix overlap queries on `QuadTree` using containment instead of overlap

## [1.0.7] - 2025-02-11

### Added

- Add faster `FeatureMap` search and sort by caching ordering statistics

## [1.0.6] - 2025-01-25

### Added

- Add generalized `CoordinateLike` impl to `NDFeatureAdapter`

## [1.0.5] - 2025-01-25

### Added

- Add `Default` impl for `NDFeatureAdapter`

## [1.0.4] - 2025-01-25

### Changed

- Pass peaks by reference in `BuildFromPeak::push_peak` to mirror `FeatureLikeMut::push`

## [1.0.3] - 2025-01-24

### Added

- Add `KnownChargeMut` to `ChargedFeature`
- Add more `NDFeature` compatibility via `NDFeatureAdapter`
- Add `get_item_unchecked` to `FeatureMapLike`
- Add `From` impls for `NDFeatureAdapter`

### Changed

- Rewrite `NDFeature` in terms of const generic arrays, recast `PeakSeries` as combination of iterate over peaks and build from peaks

## [1.0.2] - 2025-01-20

### Added

- Add `PeakSeries` and `BuildFromPeak` traits for features
- Add `rayon` feature for `PeakSetVec` and `FeatureMap` `par_iter` methods

## [1.0.1] - 2024-12-24

### Added

- Add `contains_iter_mut` and `overlaps_iter_mut` to `IntervalTree`
- Add `clear` to `FeatureLikeMut`

## [1.0.0] - 2024-12-13

### Added

- Add approximate equality to geometry types
- Add draft of *n*-ary `Feature`-like type
- Add more serde support throughout the type system

### Changed

- `PeakCollection::iter` lifetime bound is no longer `'static`

### Fixed

- Make automatic comparator macro for peak types more consistent
- More testing of geometry types
- Fix `SplittableFeatureLike` slicing behavior
- Fix up `IntervalTree` and test coverage

### Removed

- Remove `Hash` implementation from implementation macros
- Remove `peak_index` submodule. Useless concept
- Remove `BetweenIter`, no value over slice

## [0.23.0] - 2024-11-10

### Added

- Add `first` and `last` methods, as well as mutable variants to `FeatureLike` and `FeatureLikeMut`
	- Also changes the return type for `FeatureLikeMut::at_mut` to
	  prevent mutating the time dimension.
- Add `HasProximity` trait to the coordinate collection type bounds

## [0.22.0] - 2024-10-14

### Fixed

- Revert `ExactSizeIterator` bound

## [0.21.0] - 2024-10-06

### Added

- Add `IndexedCoordinate` impl for `&mut T` where `T: IndexedCoordinate`
- Add `intensity_view` to `TimeArray` trait

### Changed

- `PeakCollection` now directly enforces `ops::Index` prepreq
- Replace `PeakSetVec` iterator types are now just aliases for `std::slice::Iter` types
- `FeatureMapLike` now directly enforces `ops::Index` prepreq

## [0.20.0] - 2024-08-29

### Added

- Add rough draft of `QuadTree`
- Add `CentroidRef` and `DeconvolutedCentroidRef` for wrapping peak points indexed independently
- Add `FusedIterator` to Feature-type iterators

### Fixed

- Fix bug in `IntervalTree.overlaps` and harmonize API with query iterators

## [0.19.0] - 2024-08-09

### Added

- Add `first` and `last` to `PeakSetVec`

## [0.18.0] - 2024-08-09

### Added

- Add marker traits to `PeakSetVec` iterators

### Changed

- Introduce loop breaking in `IntervalTree::new`

### Fixed

- Fix infinite loop bug in `IntervalTree::new`

### Removed

- Remove `Deref` implementation from `PeakSetVec`

## [0.17.0] - 2024-07-12

### Added

- Add missing `FromIterator` for `SimpleFeature`
- Add `name` to `CoordinateSystem`

## [0.16.0] - 2024-06-26

### Added

- Add `into_inner` to owning `FeatureLike` types to allows these types to expose their internal data

### Fixed

- Fix macro implementer expecting symbols to be in global scope

## [0.14.0] - 2024-06-05

### Added

- Support `Span1D` and `Span2D` on references to implementing types.

### Fixed

- Fix time search behaviors and mutation behaviors of `Feature` types

## [0.13.0] - 2024-05-24

### Changed

- Refactor coodinates and implement interval tree

[1.0.8]: https://github.com/mobiusklein/mzpeaks/compare/v1.0.7..v1.0.8
[1.0.7]: https://github.com/mobiusklein/mzpeaks/compare/v1.0.6..v1.0.7
[1.0.6]: https://github.com/mobiusklein/mzpeaks/compare/v1.0.5..v1.0.6
[1.0.5]: https://github.com/mobiusklein/mzpeaks/compare/v1.0.4..v1.0.5
[1.0.4]: https://github.com/mobiusklein/mzpeaks/compare/v1.0.3..v1.0.4
[1.0.3]: https://github.com/mobiusklein/mzpeaks/compare/v1.0.2..v1.0.3
[1.0.2]: https://github.com/mobiusklein/mzpeaks/compare/v1.0.1..v1.0.2
[1.0.1]: https://github.com/mobiusklein/mzpeaks/compare/v1.0.0..v1.0.1
[1.0.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.23.0..v1.0.0
[0.23.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.22.0..v0.23.0
[0.22.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.21.0..v0.22.0
[0.21.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.20.0..v0.21.0
[0.20.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.19.0..v0.20.0
[0.19.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.18.0..v0.19.0
[0.18.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.17.0..v0.18.0
[0.17.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.16.0..v0.17.0
[0.16.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.15.0..v0.16.0
[0.14.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.13.0..v0.14.0
[0.13.0]: https://github.com/mobiusklein/mzpeaks/compare/v0.12.0..v0.13.0

<!-- generated by git-cliff -->
