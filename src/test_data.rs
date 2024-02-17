use std::fs;
use std::io;
use std::io::prelude::*;

use crate::peak::CentroidPeak;
use crate::peak_set::{PeakCollectionMut, PeakSet};

pub(crate) fn read_peaks_from_file(path: &str) -> io::Result<PeakSet> {
    let reader = io::BufReader::new(fs::File::open(path)?);
    let mut peak_set = PeakSet::empty();
    for (i, line) in reader.lines().enumerate() {
        let line = line.unwrap();
        let pref = line.trim();
        let chunks: Vec<&str> = pref.split(' ').collect();
        let mz = chunks[0].parse::<f64>().expect("Expected number for m/z");
        let intensity = chunks[1]
            .parse::<f32>()
            .expect("Expected number for intensity");
        let peak = CentroidPeak {
            mz,
            intensity,
            index: i as u32,
        };
        peak_set.push(peak);
    }
    Ok(peak_set)
}
