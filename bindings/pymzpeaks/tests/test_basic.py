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

