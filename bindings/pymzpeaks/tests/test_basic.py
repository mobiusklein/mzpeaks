import math
import pytest


def test_import():
    import pymzpeaks as m
    assert m.__name__ == "pymzpeaks"


def test_centroid_peak_and_set_roundtrip():
    import pymzpeaks as m

    p0 = m.CentroidPeak(186.04, 522.0, 0)
    p1 = m.CentroidPeak(204.07, 9800.0, 1)
    p2 = m.CentroidPeak(205.07, 150.0, 2)

    ps = m.PeakSet([p0, p1, p2])

    assert len(ps) == 3

    # __getitem__ path and reference object
    r1 = ps[1]
    assert math.isclose(r1.mz, 204.07, rel_tol=0.0, abs_tol=1e-6)
    assert r1.index == 1
    assert math.isclose(r1.intensity, 9800.0, rel_tol=0.0, abs_tol=1e-6)

    # has_peak returns a reference into the set
    q = 204.05
    tol = m.Tolerance("0.02Da")
    pref = ps.has_peak(q, tol)
    assert pref
    assert math.isclose(pref.mz, 204.07, rel_tol=0.0, abs_tol=1e-6)

    # between returns a list of references
    block = ps.between(200.0, 206.0, 10.0)  # 10 PPM by float -> default PPM
    assert isinstance(block, list)
    assert len(block) == 2

    # all_peaks_for with PPM
    match = ps.all_peaks_for(204.07, 10.0)
    assert match
    assert isinstance(match, list)
    assert all(hasattr(p, "mz") for p in match)

    match = ps.all_peaks_for(204.05, 10.0)
    assert not match
    assert isinstance(match, list)
    assert all(hasattr(p, "mz") for p in match)


def test_deconvoluted_peak_set():
    import pymzpeaks as m

    p0 = m.DeconvolutedPeak(1000.0, 100.0, 2, 0)
    p1 = m.DeconvolutedPeak(1500.0, 50.0, 3, 1)
    ds = m.DeconvolutedPeakSet([p0, p1])

    assert len(ds) == 2
    assert isinstance(ds[0], m.DeconvolutedPeak)

    # Query API
    assert ds.has_peak(1000.0, 10.0) is not None
    allp = ds.all_peaks_for(1500.0, m.Tolerance("10ppm"))
    assert isinstance(allp, list)


def test_centroid_and_deconvoluted_repr_and_ordering():
    import pymzpeaks as m

    p1 = m.CentroidPeak(204.07, 5000.0, 1)
    p2 = m.CentroidPeak(205.07, 5000.0, 2)
    assert "CentroidPeak(" in repr(p1)
    assert p1 < p2
    assert p2 > p1
    assert p1 == p1
    assert p1 != p2

    d1 = m.DeconvolutedPeak(799.359964027, 5000.0, 2, 1)
    d2 = m.DeconvolutedPeak(1000.0, 100.0, 2, 2)
    assert "DeconvolutedPeak(" in repr(d1)
    assert d1 != d2


def test_deconvoluted_peak_mz_calculation():
    import pymzpeaks as m

    neutral_mass = 799.359964027
    charge = 2
    x = m.DeconvolutedPeak(neutral_mass, 5000.0, charge, 1)
    charge_carrier = 1.007276
    expected_mz = (neutral_mass + charge_carrier * charge) / charge
    assert math.isclose(x.mz, expected_mz, rel_tol=0.0, abs_tol=1e-6)


def test_peakset_sequence_behavior_from_file():
    import pymzpeaks as m

    from pathlib import Path
    repo_root = Path(__file__).resolve().parents[3]
    data_path = repo_root / "test" / "data" / "test.txt"
    peaks = []
    with open(data_path, "rt") as fh:
        for i, line in enumerate(fh):
            mz_s, inten_s = line.strip().split()
            peaks.append(m.CentroidPeak(float(mz_s), float(inten_s), i))

    ps = m.PeakSet(peaks)

    assert len(ps) == len(peaks)
    assert math.isclose(ps[0].mz, 231.3888, rel_tol=0.0, abs_tol=1e-3)

    # has_peak and all_peaks_for
    pref = ps.has_peak(773.4414, m.Tolerance("0.01Da"))
    assert pref is not None
    assert pref.index == 300

    block = ps.all_peaks_for(736.637, m.Tolerance("10ppm"))
    assert isinstance(block, list)
    assert len(block) == 1
    block = ps.all_peaks_for(736.237, m.Tolerance("10ppm"))
    assert len(block) == 0

    q = 1221.639893
    block = ps.all_peaks_for(q, m.Tolerance("0.5Da"))
    assert len(block) == 1

    q = 2000.0
    block = ps.all_peaks_for(q, 10.0)
    assert len(block) == 0

    # between
    block = ps.between(-2000.0, 2000.0, 10.0)
    assert len(block) == len(ps)
    block = ps.between(0.0, 2000.0, 10.0)
    assert len(block) == len(ps)
    block = ps.between(1313.0, 1316.0, m.Tolerance("10ppm"))
    assert len(block) == 3


def test_peakset_edgecases():
    import pymzpeaks as m

    ps = m.PeakSet([m.CentroidPeak(500.0, 2.0, 0)])
    assert ps.has_peak(500.0, m.Tolerance("1.0Da")) is not None
    assert len(ps.all_peaks_for(500.0, m.Tolerance("1.0Da"))) == 1

    ps_empty = m.PeakSet([])
    assert ps_empty.has_peak(500.0, m.Tolerance("1.0Da")) is None
    assert len(ps_empty.all_peaks_for(500.0, m.Tolerance("1.0Da"))) == 0


def test_tolerance_use_and_multiplication_in_queries():
    import pymzpeaks as m

    # Base 1 ppm at 1e6 has width 1; scaled by 10 gives width 10
    tol = m.Tolerance("1ppm")
    tol10 = tol * 10.0

    p_lo = m.CentroidPeak(1_000_000.0 - 9.0, 100.0, 0)
    p_hi = m.CentroidPeak(1_000_000.0 + 9.0, 200.0, 1)
    ps = m.PeakSet([p_lo, p_hi])

    pref = ps.has_peak(1_000_000.0, tol10)
    assert pref is not None


def test_peakref_setters_resort_and_mutate():
    import pymzpeaks as m

    p0 = m.CentroidPeak(100.0, 10.0, 0)
    p1 = m.CentroidPeak(200.0, 20.0, 1)
    p2 = m.CentroidPeak(300.0, 30.0, 2)
    ps = m.PeakSet([p0, p1, p2])

    # Grab middle ref and move it to the end by mz
    r1 = ps[1]
    r1.mz = 350.0
    # After resort, the largest mz should be at the last index
    assert math.isclose(ps[2].mz, 350.0, rel_tol=0.0, abs_tol=1e-6)

    # Mutate intensity through the ref and read it back
    r1.intensity = 999.0
    assert math.isclose(ps[2].intensity, 999.0, rel_tol=0.0, abs_tol=1e-6)

