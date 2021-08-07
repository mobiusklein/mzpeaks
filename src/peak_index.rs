use std::collections::hash_map::{HashMap};
use std::marker::PhantomData;

use crate::prelude::*;


#[derive(Debug, Default)]
pub struct CoordinateKey(f64);

impl PartialEq for CoordinateKey {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).abs() < 1e-6
    }
}

impl Eq for CoordinateKey {}

impl PartialEq<f64> for CoordinateKey {
    fn eq(&self, other: &f64) -> bool {
        (self.0 - other).abs() < 1e-6
    }
}


impl std::hash::Hash for CoordinateKey {
    #[inline]
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let v: i64 = self.0.round() as i64;
        v.hash(state);
    }
}


#[derive(Default, Debug)]
pub struct PeakSliceMap<'lifespan, T: CoordinateLike<C>, C> {
    pub map: HashMap<CoordinateKey, &'lifespan [T]>,
    _1: PhantomData<C>,
}


impl<'transient, 'lifespan: 'transient, T: CoordinateLike<C>, C> PeakSliceMap<'lifespan, T, C> {
    pub fn with_capacity(capacity: usize) -> PeakSliceMap<'lifespan, T, C> {
        PeakSliceMap {
            map: HashMap::with_capacity(capacity),
            _1: PhantomData
        }
    }

    pub fn get(&'lifespan self, key: f64) -> Option<&&'lifespan [T]> {
        let k = CoordinateKey(key);
        self.map.get(&k)
    }

    pub fn insert(& mut self, key: f64, chunk: &'lifespan [T]) {
        let k = CoordinateKey(key);
        self.map.insert(k, chunk);
    }

    pub fn clear(&mut self) {
        self.map.clear()
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use crate::{CentroidPeak, MZ, MassErrorType};
    use crate::test_data;
    use std::io;

    #[test]
    fn test_slice_map() -> io::Result<()> {
        let mut peaks = test_data::read_peaks_from_file("./test/data/test.txt")?;

        let mut map = PeakSliceMap::<CentroidPeak, MZ>::default();

        let q = 1221.639893;
        let block = peaks.all_peaks_for(q, 0.4, MassErrorType::Absolute);
        assert_eq!(block.len(), 1);

        map.insert(q, block);
        let block2 = map.get(q).unwrap();
        assert_eq!(block2.len(), 1);

        let i = block2[0].get_index();
        let val = block2[0].intensity();
        peaks[i as usize].intensity *= 10.0;
        assert_ne!(val, peaks[i as usize].intensity());
        Ok(())
    }
}
