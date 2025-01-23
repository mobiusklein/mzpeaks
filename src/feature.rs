//! A feature is a two dimensional mass spectrum concept for a measure over some time unit.
//! It represents something that is located at a constrained but varying coordinate system `X`
//! over a sequentially ordered dimension `Y` with an abundance measure at each time point.
//!

mod charged;
mod feature;
mod ndim;
mod simple;
mod traits;
mod util;

pub use traits::{
    AsPeakIter, BuildFromPeak, FeatureLike, FeatureLikeMut, NDFeatureLike,
    NDFeatureLikeMut, PeakSeries, SplittableFeatureLike, TimeArray, TimeInterval,
};

pub use charged::{
    Charged as ChargedPointWrapper, ChargedFeature, ChargedFeatureView, ChargedFeatureWrapper,
    DeconvolutedPeakIter, DeconvolvedIMSFeature, DeconvolvedLCMSFeature,
};

pub use feature::{
    Feature, FeatureView, IMSFeature, IntoIter, Iter, IterMut, LCMSFeature, MZPeakIter,
};

pub use simple::{SimpleFeature, SimpleFeatureView};

pub use ndim::{NDFeature, NDIter, NDIterMut, NDPoint, NDPointMutRef, NDFeatureAdapter};

#[cfg(test)]
mod test {
    use super::*;
    use crate::{prelude::*, Mass};
    use crate::{CentroidPeak, DeconvolutedPeak, MZLocated, Time, MZ};

    #[test]
    fn test_build_raw() {
        let mut x = LCMSFeature::empty();

        let points = vec![
            (CentroidPeak::new(204.08, 3432.1, 0), 0.1),
            (CentroidPeak::new(204.07, 7251.9, 0), 0.2),
            (CentroidPeak::new(204.08, 5261.7, 0), 0.3),
        ];

        let (idx, _) = LCMSFeature::empty().find_time(128.470);
        assert_eq!(idx, None);

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

        let (x, y, z) = x.into_inner();
        let x = LCMSFeature::new(x, y, z);

        let mut y = LCMSFeature::with_capacity(3);

        for (p, t) in x.iter_peaks().rev() {
            y.push(&p, t);
        }
        assert_eq!(y.apex_time(), Some(0.2));
        assert_eq!(y.as_view().to_owned(), y);
    }

    #[test]
    fn test_empty_behavior() {
        macro_rules! do_test {
            ($x:expr) => {
                assert!($x.apex_time().is_none());
                assert!($x.is_empty());

                let (a, b) = $x.split_at(5);
                assert!(a.is_empty());
                assert!(b.is_empty());

                let (a, b) = $x.split_at_time(50.0);
                assert!(a.is_empty());
                assert!(b.is_empty());

                assert_eq!($x.area(), 0.0);
                assert_eq!($x.intensity(), 0.0);

                assert_eq!($x.iter().len(), 0);

                let mut i = 0;
                for _ in $x {
                    i += 1;
                }
                assert_eq!(i, 0);
            };
        }

        let mut x = LCMSFeature::empty();
        assert_eq!(x.iter_mut().len(), 0);
        assert_eq!(x.iter_peaks().len(), 0);
        do_test!(x);
        let mut x: ChargedFeature<crate::Mass, Time> = ChargedFeature::empty(2);

        assert_eq!(x.iter_mut().len(), 0);
        assert_eq!(x.iter_peaks().len(), 0);
        do_test!(x);

        let x: FeatureView<'_, MZ, Time> = FeatureView::empty();
        do_test!(x);

        let x: ChargedFeatureView<'_, MZ, Time> = ChargedFeatureView::empty(2);
        do_test!(x);
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

    fn behaviors_to_test_mut<
        'a,
        F: FeatureLike<MZ, Time> + FeatureLikeMut<MZ, Time> + SplittableFeatureLike<'a, MZ, Time>,
    >(
        feature: &'a mut F,
    ) {
        let y = feature.intensity();
        feature.iter_mut().for_each(|(_, _, z)| {
            *z *= 2.0;
        });
        let y2 = feature.intensity();
        assert_eq!(y * 2.0, y2);
        for i in 0..feature.len() {
            let pt = feature.at_mut(i).unwrap();
            *pt.2 /= 2.0;
        }
        let y3 = feature.intensity();
        assert_eq!(y, y3);
        let times: Vec<_> = feature.iter_time().collect();
        for t in times.iter().copied() {
            let pt = feature.at_time_mut(t).unwrap();
            *pt.2 *= 2.0;
        }
        let y4 = feature.intensity();
        assert_eq!(y * 2.0, y4);

        let y5: f32 = times
            .iter()
            .copied()
            .flat_map(|t| feature.at_time(t))
            .map(|pt| pt.2)
            .sum();
        assert_eq!(y * 2.0, y5);

        let y6: f32 = (0..feature.len())
            .flat_map(|i| feature.at(i))
            .map(|pt| pt.2)
            .sum();
        assert_eq!(y * 2.0, y6);

        let n = feature.len();
        let (a, b, c) = feature.last().unwrap();
        feature.push_raw(a, b, c);
        assert_eq!(
            feature.len(),
            n,
            "Appending to the trailing time point should not add a new entry"
        );

        let (a, b, c) = feature.first().unwrap();
        feature.push_raw(a, b, c);
        assert_eq!(
            feature.len(),
            n + 1,
            "Appending to the a new time point should add a new entry"
        );
    }

    fn behaviors_to_test<
        'a,
        F: FeatureLike<MZ, Time> + SplittableFeatureLike<'a, MZ, Time> + TimeArray<Time>,
    >(
        feature: &'a F,
    ) {
        let err = feature.mz() - 1075.229;
        assert!(err.abs() < 1e-3, "Error = {err}");

        let (idx, err) = feature.find_time(128.470);

        assert_eq!(idx, Some(9));
        assert!(feature.spans(128.47));

        assert_eq!(
            feature.as_range(),
            crate::coordinate::CoordinateRange::new(feature.start_time(), feature.end_time())
        );

        let e = 128.4700598232 - 128.470;
        assert!((e - err).abs() < 1e-3);

        let (idx, _) = feature.find_time(0.0);
        assert_eq!(idx, Some(0));

        let (idx, _) = feature.find_time(1000.0);
        assert_eq!(idx, Some(feature.len() - 1));

        let part = feature.slice(0..5);
        assert_eq!(part.end_time(), feature.at(4).map(|(_, b, _)| b));
        let part = feature.slice(..5);
        assert_eq!(part.end_time(), feature.at(4).map(|(_, b, _)| b));

        let part = feature.slice(..=5);
        assert_eq!(part.end_time(), feature.at(5).map(|(_, b, _)| b));

        let part = feature.slice(0..=5);
        assert_eq!(part.end_time(), feature.at(5).map(|(_, b, _)| b));

        let part = feature.slice(..);
        assert_eq!(part.end_time(), feature.last().map(|(_, b, _)| b));
        assert_eq!(part.start_time(), feature.first().map(|(_, b, _)| b));

        let part = feature.slice(2..);
        assert_eq!(part.end_time(), feature.last().map(|(_, b, _)| b));

        let mask: Vec<_> = feature.iter_time().map(|t| t >= 128.4700598232).collect();
        let parts = feature.split_mask(&mask);
        assert_eq!(parts.len(), 1);
        assert_eq!(parts[0].end_time(), feature.end_time());

        let parts = feature.split_when(|(_, y, _), _| y <= 128.4700598232);
        assert_eq!(parts.len(), 1);
        assert_eq!(parts[0].end_time(), feature.end_time());

        let (part_a, part_b) = feature.split_at_time(128.4700598232);
        assert_eq!(part_b.end_time(), feature.end_time());
        assert_eq!(part_a.start_time(), feature.start_time());

        let (part_a, part_b) = feature.split_at(9);
        assert_eq!(part_b.end_time(), feature.end_time());
        assert_eq!(part_a.start_time(), feature.start_time());

        assert_eq!(feature.time_view().len(), feature.len());
        assert_eq!(feature.intensity_view().len(), feature.len());
    }

    #[test]
    fn test_behaviors() {
        let mut feature = make_large_feature();
        behaviors_to_test(&feature);
        behaviors_to_test(&feature.as_view());
        behaviors_to_test_mut(&mut feature);

        let (_pt, t) = feature.iter_peaks().nth(10).unwrap();
        let e = 128.4700598232 - t;
        assert!(e.abs() < 1e-3);
        assert!(feature.iter_mut().next_back().is_some());
        assert!(feature.iter_peaks().next_back().is_some());

        let mut feature = ChargedFeature::new(make_large_feature(), 2);
        behaviors_to_test(&feature);
        behaviors_to_test(&feature.as_view());
        behaviors_to_test_mut(&mut feature);
        assert!(feature.iter_mut().next_back().is_some());
        assert_eq!(feature.charge(), 2);
        let mut feature2: ChargedFeature<Mass, Time> = ChargedFeature::empty(2);
        feature2.extend(feature.into_iter());
        let (_pt, t) = feature2.iter_peaks().nth(9).unwrap();
        let e = 128.4700598232 - t;
        assert!(e.abs() < 1e-3);
        assert!(feature2.iter_mut().next_back().is_some());
        assert!(feature2.iter_peaks().next_back().is_some());

        assert_eq!(feature2.as_inner().1, feature2.charge());
        assert_eq!(feature2.as_view(), feature2.as_view());
        let (f, z) = feature2.clone().into_inner();
        assert_eq!(z, 2);
        assert_eq!(f, *feature2.as_ref())
    }

    #[test]
    fn test_behaviors_simple() {
        let feature = make_large_feature();
        let x = feature.coordinate_x();
        let (_, y, z) = feature.into_inner();

        let mut feature: SimpleFeature<MZ, Time> = SimpleFeature::new(x, y, z);
        let f2 = feature.clone();
        assert_eq!(feature, f2);
        assert!(feature >= f2);
        assert!(feature <= f2);

        assert_eq!(feature.as_view(), f2.as_view());
        assert!(feature.as_view() >= f2.as_view());
        assert!(feature.as_view() <= f2.as_view());

        behaviors_to_test(&feature);
        behaviors_to_test(&feature.as_view());
        behaviors_to_test_mut(&mut feature);

        assert_eq!(feature.apex_time(), f2.as_view().apex_time());

        let feature: SimpleFeature<MZ, Time> = f2.into_iter().collect();
        behaviors_to_test(&feature);
    }
}
