# mzpeaks
[![Latest Version](https://img.shields.io/crates/v/mzpeaks.svg)](https://crates.io/crates/mzpeaks)

`mzpeaks` implements the building blocks and machinery for representing peaks
in a mass spectrum.

It's meant to be used as a building block for other tools and does not provide
any I/O machinery for peak lists.

If you're looking for the mzPeak file format, you may want to instead look at https://github.com/mobiusklein/mzpeak_prototyping.


## Usage

```rust
use mzpeaks::{CentroidPeak, PeakSet, PeakCollection, Tolerance};

let peaks = PeakSet::new(vec![
    CentroidPeak::new(186.04, 522.0, 0),
    CentroidPeak::new(204.07, 9800.0, 1),
    CentroidPeak::new(205.07, 150.0, 2)
]);

assert_eq!(peaks.search(204.05, Tolerance::Da(0.02)).unwrap(), 1);

let peak = match peaks.has_peak(204.05, Tolerance::Da(0.02)) {
    Some(p) => p,
    None => panic!("Failed to retrieve peak!")
};

assert!((peak.mz - 204.07).abs() < 1e-6);
```

## Compile time features

`mzpeask` can be modified with

- `serde`: Adds `serde::Serialize` and `serde::Deserialize` implementations to most types in the library
- `rayon`: Adds `par_iter` and similar methods to `PeakSetVec` and `FeatureMap`


## Python Bindings

The Python bindings located in `bindings/pymzpeaks` can be built with `maturin`. They provide very basic access to
`PeakSetVec`-like and peak-like types via `CentroidPeak`,  `DeconvolutedPeak`, `PeakSet`, and `DeconvolutedPeakSet`.
These bindings are not expected to be useful, but were written as part of a learning process.