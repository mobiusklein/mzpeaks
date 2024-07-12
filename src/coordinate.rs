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

    #[test]
    fn test_coordinate_system() {
        let mut peak = DeconvolutedPeak::new(204.09, 300.0, 2, 0);
        let mass = <Mass as CoordinateSystem>::coordinate(&peak);
        let mz = MZ::coordinate(&peak);

        assert_eq!(peak.neutral_mass(), mass);
        assert_eq!(peak.mz(), mz);

        *Mass::coordinate_mut(&mut peak) = 9001.0;
    }
}
