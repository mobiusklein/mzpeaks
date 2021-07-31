# mzpeaks
[![Latest Version](https://img.shields.io/crates/v/mzpeaks.svg)](https://crates.io/crates/mzpeaks)

`mzpeaks` implements the building blocks and machinery for representing peaks
in a mass spectrum.

It's meant to be used as a building block for other tools and does not provide
any I/O machinery for peak lists


## Usage

```rust
use mzpeaks::{CentroidPeak, PeakSet, PeakCollection, MassErrorType};

let peaks = PeakSet::new(vec![
    CentroidPeak::new(186.04, 522.0, 0),
    CentroidPeak::new(204.07, 9800.0, 1),
    CentroidPeak::new(205.07, 150.0, 2)
]);

assert_eq!(peaks.search(204.05, 0.02, MassErrorType::Absolute).unwrap(), 1);

let peak = match peaks.has_peak(204.05, 0.02, MassErrorType::Absolute) {
    Some(p) => p,
    None => panic!("Failed to retrieve peak!")
};

assert!((peak.mz - 204.07).abs() < 1e-6);
```