# mzpeaks
[![Latest Version](https://img.shields.io/crates/v/mzpeaks.svg)](https://crates.io/crates/mzpeaks)

`mzpeaks` implements the building blocks and machinery for representing peaks
in a mass spectrum.

It's meant to be used as a building block for other tools and does not provide
any I/O machinery for peak lists


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

## Building the Rust crate

### Prerequisites

- Rust toolchain (via `rustup`) and build tools

```bash
sudo apt-get update && sudo apt-get install -y build-essential
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

### Build and test

```bash
cd mzpeaks
cargo build --release
cargo test
```

### Optional crate features

```bash
# Enable serde support
cargo build --release --features serde-support

# Enable parallelism via rayon
cargo build --release --features rayon

# Combine features
cargo build --release --features "serde-support rayon"
```

### Build docs (optional)

```bash
cargo doc --lib --no-deps
```

## Python bindings (pymzpeaks)

The Python bindings are located in `bindings/pymzpeaks` and built with `maturin`/PyO3.

### Prerequisites

- Python 3.7+
- Python headers on Linux: `sudo apt-get install -y python3-dev`

### Option A — Develop install with maturin (recommended)

```bash
# Activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install maturin as specified in pyproject
pip install "maturin[patchelf]>=0.14,<0.15"

# Build and install into the active venv
cd mzpeaks/bindings/pymzpeaks
maturin develop --release

# Optional: enable crate features for the Python build
maturin develop --release --features serde-support
```

### Option B — Install via pip (PEP 517 build)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install mzpeaks/bindings/pymzpeaks
```

### Build a wheel (without installing)

```bash
cd mzpeaks/bindings/pymzpeaks
maturin build --release
# Wheels will be under mzpeaks/bindings/pymzpeaks/target/wheels
```

### Quick verification

```bash
python -c "import pymzpeaks as m; print('ok', m.__name__)"
```