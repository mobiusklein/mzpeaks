//! A type system implementation of a coordinate system that attempts to deal with the different dimensions
//! an observation may be placed in simultaneously.
mod dim;
mod range;
mod bbox;
mod ivtree;

pub use dim::{
    CoordinateLike, CoordinateLikeMut, CoordinateSystem,
    IndexType, IndexedCoordinate,
    IonMobility, IonMobilityLocated,
    MZ, MZLocated,
    Mass, MassLocated,
    Time, TimeLocated,
};

pub use range::{CoordinateRange, CoordinateRangeParseError, SimpleInterval, Span1D};
pub use bbox::{BoundingBox, CoordinateBox, Span2D};
pub use ivtree::{IntervalTree, IntervalTreeNode, PreorderIter};

#[cfg(test)]
mod test {
    use super::*;
    use crate::DeconvolutedPeak;

    fn check_coord<C: CoordinateSystem, P: CoordinateLike<C>>(peak: &P, cs: &C) -> f64 {
        cs.coordinate(peak)
    }

    #[test]
    fn test_coordinate_system() {
        let mut peak = DeconvolutedPeak::new(204.09, 300.0, 2, 0);
        let mass = check_coord(&peak, &Mass());
        let mz = check_coord(&peak, &MZ());

        assert_eq!(peak.neutral_mass(), mass);
        assert_eq!(peak.mz(), mz);

        *Mass().coordinate_mut(&mut peak) = 9001.0;
    }
}
