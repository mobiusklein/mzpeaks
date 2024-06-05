//! A feature is a two dimensional mass spectrum concept for a measure over some time unit.
//! It represents something that is located at a constrained but varying coordinate system `X`
//! over a sequentially ordered dimension `Y` with an abundance measure at each time point.
//!

mod charged;
mod feature;
mod simple;
mod traits;
mod util;

pub use traits::{FeatureLike, FeatureLikeMut, SplittableFeatureLike, TimeInterval};

pub use charged::{
    ChargedFeature, ChargedFeatureView, DeconvolutedPeakIter, DeconvolvedIMSFeature,
    DeconvolvedLCMSFeature,
};
pub use feature::{
    Feature, FeatureView, IMSFeature, IntoIter, Iter, IterMut, LCMSFeature, MZPeakIter,
};
pub use simple::{SimpleFeature, SimpleFeatureView};

#[cfg(test)]
mod test {
    use super::*;
    use crate::prelude::*;
    use crate::{CentroidPeak, DeconvolutedPeak, MZLocated, Time, MZ};

    #[test]
    fn test_build_raw() {
        let mut x = LCMSFeature::empty();

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

        if let Some((pt, t)) = x.iter_peaks().next() {
            let (rpt, rt) = &points[0];
            assert_eq!(pt, rpt);
            assert_eq!(t, *rt);
        }

        assert_eq!(x, points.into_iter().collect());

        let (i, e) = x.find_time(0.3);
        assert_eq!(i, Some(2));
        assert_eq!(e, 0.0);

        let (i, e) = x.find_time(0.5);
        assert_eq!(i, Some(2));
        assert_eq!(e, 0.2);

        let (b, a) = x.split_at_time(0.2);
        assert_eq!(b.len(), 1);
        assert_eq!(a.len(), 2);
    }

    #[test]
    fn test_build_charged() {
        let mut x = DeconvolvedLCMSFeature::empty(1);

        let points = vec![
            (DeconvolutedPeak::new(203.08, 3432.1, 1, 0), 0.1),
            (DeconvolutedPeak::new(203.07, 7251.9, 1, 0), 0.2),
            (DeconvolutedPeak::new(203.08, 5261.7, 1, 0), 0.3),
        ];

        x.extend(points.iter().cloned());

        assert_eq!(x.len(), 3);

        let y: f32 = points.iter().map(|p| p.0.intensity()).sum();
        assert!((x.intensity() - y).abs() < 1e-6);

        let mass = 203.07545212;
        assert!((x.neutral_mass() - mass).abs() < 1e-6);

        let area = 1159.879980;
        assert!((x.area() - area).abs() < 1e-6);

        assert_eq!(x.iter().len(), 3);

        if let Some((pt, t)) = x.iter_peaks().next() {
            let (rpt, rt) = &points[0];
            assert_eq!(pt, rpt);
            assert_eq!(t, *rt);
        }

        assert_eq!(x, points.into_iter().collect());
    }

    fn make_large_feature() -> Feature<MZ, Time> {
        let points = [
            (1075.227783203125, 127.814165801867, 50647.84),
            (1075.2315673828125, 127.8677471344, 235570.69),
            (1075.228515625, 127.921503501067, 143257.8),
            (1075.2294921875, 127.975274319467, 123955.15),
            (1075.225341796875, 128.136723199467, 79086.79),
            (1075.2261962890625, 128.1903371096, 221846.45),
            (1075.2293701171875, 128.220671521333, 76135.8),
            (1075.228271484375, 128.3669475152, 372242.47),
            (1075.2281494140625, 128.416172280267, 225322.0),
            (1075.22900390625, 128.4700598232, 41362.75),
            (1075.2318115234375, 128.5238657408, 131210.38),
            (1075.225830078125, 128.5777638256, 146266.61),
            (1075.2281494140625, 128.631568028533, 111895.05),
            (1075.2293701171875, 128.6850992784, 152425.03),
            (1075.2303466796875, 128.738557833067, 215884.19),
            (1075.227294921875, 128.792116063733, 125404.945),
            (1075.22705078125, 128.8457370032, 119783.586),
            (1075.2257080078125, 128.895031093867, 244697.55),
            (1075.2279052734375, 128.948732996267, 119553.375),
            (1075.228515625, 129.002493416267, 164899.05),
            (1075.227783203125, 129.056439435733, 49560.582),
            (1075.22705078125, 129.110237865867, 262365.28),
            (1075.226318359375, 129.163855989867, 263023.44),
            (1075.22900390625, 129.2174601552, 130926.49),
            (1075.2274169921875, 129.271160236, 285440.4),
            (1075.2249755859375, 129.3246664664, 142812.7),
            (1075.22998046875, 129.373630130933, 161386.63),
            (1075.2294921875, 129.4274522376, 63290.746),
            (1075.2293701171875, 129.4810838, 132746.83),
        ];
        Feature::from_iter(points.iter().copied())
    }

    #[test]
    fn test_behaviors() {
        let feature = make_large_feature();
        let err = feature.mz() - 1075.229;
        assert!(err.abs() < 1e-3, "Error = {err}");

        let (idx, err) = feature.find_time(128.470);

        assert_eq!(idx, Some(9));

        let e = 128.4700598232 - 128.470;
        assert!((e - err).abs() < 1e-3);

        let (idx, _) = feature.find_time(0.0);
        assert_eq!(idx, Some(0));

        let (idx, _) = feature.find_time(1000.0);
        assert_eq!(idx, Some(feature.len() - 1));

        let (idx, _) = LCMSFeature::empty().find_time(128.470);
        assert_eq!(idx, None);
    }
}
