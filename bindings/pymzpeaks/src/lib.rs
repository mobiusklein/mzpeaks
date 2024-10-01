use std::ops::{Deref, DerefMut, Index, IndexMut};
use std::str::FromStr;
use std::sync::{Arc, RwLock, RwLockReadGuard, RwLockWriteGuard};

use pyo3::exceptions::{PyException, PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::pyclass::CompareOp;
use pyo3::types::{PyFloat, PyList, PyLong, PySlice, PyString};

use mzpeaks::coordinate::{MZLocated, Mass, MZ};
use mzpeaks::peak_set::{PeakCollection, PeakSetVec};
use mzpeaks::Tolerance;
use mzpeaks::{prelude::*, IndexType};
use mzpeaks::{CentroidPeak, DeconvolutedPeak};

#[derive(FromPyObject)]
pub enum FloatOrTolerance {
    Float(f64),
    Tolerance(PyTolerance),
}

impl From<FloatOrTolerance> for PyTolerance {
    fn from(value: FloatOrTolerance) -> PyTolerance {
        match value {
            FloatOrTolerance::Float(f) => PyTolerance(Tolerance::PPM(f)),
            FloatOrTolerance::Tolerance(t) => t,
        }
    }
}

#[pyclass(module = "pymzpeaks", name = "Tolerance")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PyTolerance(Tolerance);

#[pymethods]
impl PyTolerance {
    #[new]
    fn new(value: Bound<PyAny>) -> PyResult<Self> {
        if value.is_instance_of::<PyFloat>() {
            Ok(Self(Tolerance::PPM(value.extract::<f64>()?)))
        } else if value.is_instance_of::<PyString>() {
            match value.extract::<&str>()?.parse::<Tolerance>() {
                Ok(tol) => Ok(Self(tol)),
                Err(err) => Err(PyValueError::new_err(err.to_string())),
            }
        } else {
            Err(PyTypeError::new_err(format!(
                "Could not convert {} into Tolerance instance",
                value.get_type().getattr("__name__")?.extract::<&str>()?
            )))
        }
    }

    fn __mul__(&self, value: f64) -> Self {
        Self(self.0 * value)
    }
}

impl FromStr for PyTolerance {
    type Err = <Tolerance as FromStr>::Err;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match <Tolerance as FromStr>::from_str(s) {
            Ok(tol) => Ok(Self(tol)),
            Err(err) => Err(err),
        }
    }
}

impl ToString for PyTolerance {
    fn to_string(&self) -> String {
        <Tolerance as ToString>::to_string(&self.0)
    }
}

impl Default for PyTolerance {
    fn default() -> Self {
        Self(Tolerance::PPM(20.0))
    }
}

#[pyclass(module = "pymzpeaks", name = "CentroidPeak")]
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct PyCentroidPeak(CentroidPeak);

#[pymethods]
impl PyCentroidPeak {
    #[new]
    pub fn new(mz: f64, intensity: f32, index: u32) -> Self {
        (CentroidPeak {
            mz,
            intensity,
            index,
        })
        .into()
    }

    #[getter]
    fn mz(&self) -> f64 {
        self.0.mz
    }

    #[getter]
    fn intensity(&self) -> f32 {
        self.0.intensity
    }

    #[getter]
    fn index(&self) -> u32 {
        return self.0.index;
    }

    fn __repr__(&self) -> String {
        format!(
            "CentroidPeak({:0.4}, {:0.4}, {})",
            self.mz(),
            self.intensity(),
            self.index()
        )
    }

    fn __richcmp__(&self, other: PyRef<PyCentroidPeak>, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(self.0 == other.0),
            CompareOp::Ne => Ok(self.0 != other.0),
            CompareOp::Le => Ok(self.0 <= other.0),
            CompareOp::Ge => Ok(self.0 >= other.0),
            CompareOp::Lt => Ok(self.0 < other.0),
            CompareOp::Gt => Ok(self.0 > other.0),
        }
    }
}

impl From<CentroidPeak> for PyCentroidPeak {
    fn from(value: CentroidPeak) -> Self {
        Self(value)
    }
}

impl CoordinateLike<MZ> for PyCentroidPeak {
    fn coordinate(&self) -> f64 {
        self.mz()
    }
}

impl IndexedCoordinate<MZ> for PyCentroidPeak {
    fn get_index(&self) -> mzpeaks::IndexType {
        self.0.index
    }

    fn set_index(&mut self, index: mzpeaks::IndexType) {
        self.0.index = index
    }
}

impl IntensityMeasurement for PyCentroidPeak {
    #[inline]
    fn intensity(&self) -> f32 {
        self.0.intensity
    }
}

impl IntensityMeasurementMut for PyCentroidPeak {
    fn intensity_mut(&mut self) -> &mut f32 {
        self.0.intensity_mut()
    }
}

#[pyclass(module = "pymzpeaks", name = "DeconvolutedPeak")]
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct PyDeconvolutedPeak(DeconvolutedPeak);

#[pymethods]
impl PyDeconvolutedPeak {
    #[new]
    pub fn new(neutral_mass: f64, intensity: f32, charge: i32, index: u32) -> Self {
        (DeconvolutedPeak {
            neutral_mass,
            intensity,
            index,
            charge,
        })
        .into()
    }

    #[getter]
    fn neutral_mass(&self) -> f64 {
        self.0.neutral_mass
    }

    #[getter]
    fn mz(&self) -> f64 {
        MZLocated::mz(&self.0)
    }

    #[getter]
    fn intensity(&self) -> f32 {
        self.0.intensity
    }

    #[getter]
    fn index(&self) -> u32 {
        return self.0.index;
    }

    #[getter]
    fn charge(&self) -> i32 {
        self.0.charge
    }

    fn __repr__(&self) -> String {
        format!(
            "DeconvolutedPeak({:0.4}, {:0.4}, {}, {})",
            self.neutral_mass(),
            self.intensity(),
            self.charge(),
            self.index()
        )
    }

    fn __richcmp__(&self, other: PyRef<PyDeconvolutedPeak>, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(self.0 == other.0),
            CompareOp::Ne => Ok(self.0 != other.0),
            CompareOp::Le => Ok(self.0 <= other.0),
            CompareOp::Ge => Ok(self.0 >= other.0),
            CompareOp::Lt => Ok(self.0 < other.0),
            CompareOp::Gt => Ok(self.0 > other.0),
        }
    }
}

impl From<DeconvolutedPeak> for PyDeconvolutedPeak {
    fn from(value: DeconvolutedPeak) -> Self {
        Self(value)
    }
}

impl CoordinateLike<MZ> for PyDeconvolutedPeak {
    fn coordinate(&self) -> f64 {
        self.0.mz()
    }
}

impl CoordinateLike<Mass> for PyDeconvolutedPeak {
    fn coordinate(&self) -> f64 {
        self.0.neutral_mass
    }
}

impl IndexedCoordinate<Mass> for PyDeconvolutedPeak {
    fn get_index(&self) -> mzpeaks::IndexType {
        self.0.index
    }

    fn set_index(&mut self, index: mzpeaks::IndexType) {
        self.0.index = index
    }
}

impl IntensityMeasurement for PyDeconvolutedPeak {
    fn intensity(&self) -> f32 {
        self.0.intensity
    }
}

impl IntensityMeasurementMut for PyDeconvolutedPeak {
    fn intensity_mut(&mut self) -> &mut f32 {
        self.0.intensity_mut()
    }
}

impl KnownCharge for PyDeconvolutedPeak {
    fn charge(&self) -> i32 {
        self.0.charge
    }
}

impl KnownChargeMut for PyDeconvolutedPeak {
    fn charge_mut(&mut self) -> &mut i32 {
        self.0.charge_mut()
    }
}

#[derive(Debug, Clone)]
pub enum CentroidPeakOrRef {
    Peak(PyCentroidPeak),
    Ref {
        index: usize,
        peak_set: Arc<RwLock<PeakSetVec<PyCentroidPeak, MZ>>>,
    },
}

impl From<CentroidPeak> for CentroidPeakOrRef {
    fn from(value: CentroidPeak) -> Self {
        Self::Peak(value.into())
    }
}

impl Default for CentroidPeakOrRef {
    fn default() -> Self {
        Self::Peak(PyCentroidPeak::new(0.0, 0.0, 0))
    }
}

#[pyclass]
#[derive(Debug, Clone)]
pub struct PyPeakRef(CentroidPeakOrRef);

#[pymethods]
impl PyPeakRef {
    #[new]
    fn new(mz: f64, intensity: f32) -> Self {
        Self(CentroidPeakOrRef::Peak(PyCentroidPeak::new(mz, intensity, 0)))
    }

    fn copy(&self) -> Self {
        match &self.0 {
            CentroidPeakOrRef::Peak(py_centroid_peak) => Self(CentroidPeakOrRef::Peak(py_centroid_peak.clone())),
            CentroidPeakOrRef::Ref { index, peak_set } => {
                let ps = peak_set.read().expect("Failed to acquire reading lock on peak set");
                let p = ps.get_item(*index);
                Self(CentroidPeakOrRef::Peak(p.clone()))
            },
        }
    }

    #[getter]
    fn mz(&self) -> PyResult<f64> {
        match &self.0 {
            CentroidPeakOrRef::Peak(peak) => Ok(peak.mz()),
            CentroidPeakOrRef::Ref { index, peak_set } => match peak_set.as_ref().read() {
                Ok(peaks) => Ok(peaks.deref().index(*index).mz()),
                Err(err) => Err(PyValueError::new_err(format!(
                    "Failed to get read access to peak set: {}",
                    err
                ))),
            },
        }
    }

    #[setter]
    fn mz_setter(&mut self, val: f64) -> PyResult<()> {
        match &mut self.0 {
            CentroidPeakOrRef::Peak(peak) => {
                peak.0.mz = val;
                Ok(())
            }
            CentroidPeakOrRef::Ref { index, peak_set } => match peak_set.as_ref().write() {
                Ok(mut peaks) => {
                    peaks.deref_mut().index_mut(*index).0.mz = val;
                    peaks.sort();
                    Ok(())
                }
                Err(err) => Err(PyValueError::new_err(format!(
                    "Failed to get read access to peak set: {}",
                    err
                ))),
            },
        }
    }

    #[getter]
    fn intensity(&self) -> PyResult<f32> {
        match &self.0 {
            CentroidPeakOrRef::Peak(peak) => Ok(peak.intensity()),
            CentroidPeakOrRef::Ref { index, peak_set } => match peak_set.as_ref().read() {
                Ok(peaks) => Ok(peaks.deref().index(*index).intensity()),
                Err(err) => Err(PyValueError::new_err(format!(
                    "Failed to get read access to peak set: {}",
                    err
                ))),
            },
        }
    }

    #[setter]
    fn intensity_setter(&mut self, val: f32) -> PyResult<()> {
        match &mut self.0 {
            CentroidPeakOrRef::Peak(peak) => {
                peak.0.intensity = val;
                Ok(())
            }
            CentroidPeakOrRef::Ref { index, peak_set } => match peak_set.as_ref().write() {
                Ok(mut peaks) => {
                    peaks.deref_mut().index_mut(*index).0.intensity = val;
                    Ok(())
                }
                Err(err) => Err(PyValueError::new_err(format!(
                    "Failed to get read access to peak set: {}",
                    err
                ))),
            },
        }
    }

    #[getter]
    fn index(&self) -> PyResult<u32> {
        match &self.0 {
            CentroidPeakOrRef::Peak(peak) => Ok(peak.get_index()),
            CentroidPeakOrRef::Ref { index, peak_set } => match peak_set.as_ref().read() {
                Ok(peaks) => Ok(peaks.deref().index(*index).get_index()),
                Err(err) => Err(PyValueError::new_err(format!(
                    "Failed to get read access to peak set: {}",
                    err
                ))),
            },
        }
    }

    #[setter]
    fn index_setter(&mut self, val: u32) -> PyResult<()> {
        match &mut self.0 {
            CentroidPeakOrRef::Peak(peak) => {
                peak.0.set_index(val);
                Ok(())
            }
            CentroidPeakOrRef::Ref { index, peak_set } => match peak_set.as_ref().write() {
                Ok(mut peaks) => {
                    peaks.deref_mut().index_mut(*index).0.set_index(val);
                    Ok(())
                }
                Err(err) => Err(PyValueError::new_err(format!(
                    "Failed to get read access to peak set: {}",
                    err
                ))),
            },
        }
    }
}

macro_rules! unpack {
    ($e:expr, $v:literal) => {
        $e.unwrap_or_else(|_| panic!("Failed to acquire read lock on peak {}", $v))
    };
    ($e:expr) => {
        $e.unwrap_or_else(|_| panic!("Failed to acquire read lock on peak"))
    };
}

impl PartialEq for PyPeakRef {
    fn eq(&self, other: &Self) -> bool {
        if (unpack!(self.mz(), 1) - unpack!(other.mz(), 2)).abs() > 1e-3 {
            false
        } else {
            (unpack!(self.intensity(), 1) - unpack!(other.intensity(), 2)).abs() > 1e-3
        }
    }
}

impl PartialOrd for PyPeakRef {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match unpack!(self.mz(), 1).partial_cmp(&unpack!(other.mz(), 2)) {
            Some(val) => match val {
                std::cmp::Ordering::Equal => {
                    unpack!(self.intensity(), 1).partial_cmp(&unpack!(other.intensity(), 2))
                }
                _ => Some(val),
            },
            None => None,
        }
    }
}

impl CoordinateLike<MZ> for PyPeakRef {
    fn coordinate(&self) -> f64 {
        unpack!(self.mz())
    }
}

impl IntensityMeasurement for PyPeakRef {
    fn intensity(&self) -> f32 {
        unpack!(self.intensity())
    }
}

impl IndexedCoordinate<MZ> for PyPeakRef {
    fn get_index(&self) -> IndexType {
        unpack!(self.index())
    }

    fn set_index(&mut self, index: IndexType) {
        self.index_setter(index)
            .unwrap_or_else(|_| panic!("Failed to acquire write lock on peak"))
    }
}

#[pyclass(module = "pymzpeaks", sequence, name = "PeakSet")]
#[derive(Debug, Clone)]
pub struct PyPeakSet(Arc<RwLock<PeakSetVec<PyCentroidPeak, MZ>>>);

impl PyPeakSet {
    pub fn read_peak_set(&self) -> PyResult<RwLockReadGuard<'_, PeakSetVec<PyCentroidPeak, MZ>>> {
        self.0.read().map_err(|e| {
            PyException::new_err(format!("Failed to acquire lock for reading peak set: {e}"))
        })
    }

    pub fn write_peak_set(
        &self,
    ) -> Result<RwLockWriteGuard<'_, PeakSetVec<PyCentroidPeak, MZ>>, PyErr> {
        self.0.write().map_err(|e| {
            PyException::new_err(format!("Failed to acquire lock for writing peak set: {e}"))
        })
    }
}

#[pymethods]
impl PyPeakSet {
    #[new]
    pub fn py_new(peaks: Vec<PyCentroidPeak>) -> Self {
        PyPeakSet(Arc::new(RwLock::new(PeakSetVec::new(peaks))))
    }

    #[pyo3(signature=(query, error_tolerance=FloatOrTolerance::Float(10.0)))]
    pub fn has_peak(
        &self,
        query: f64,
        error_tolerance: FloatOrTolerance,
    ) -> PyResult<Option<PyPeakRef>> {
        match self.0.read() {
            Ok(ps) => {
                if let Some(p) = ps.has_peak(
                    query,
                    <FloatOrTolerance as Into<PyTolerance>>::into(error_tolerance).0,
                ) {
                    let i = p.index();
                    Ok(Some(PyPeakRef(CentroidPeakOrRef::Ref {
                        index: i as usize,
                        peak_set: self.0.clone(),
                    })))
                } else {
                    Ok(None)
                }
            }
            Err(e) => Err(PyException::new_err(format!(
                "Failed to acquire lock for reading peak set: {e}"
            ))),
        }
    }

    #[pyo3(signature=(low, high, error_tolerance=FloatOrTolerance::Float(10.0)))]
    pub fn between(
        &self,
        low: f64,
        high: f64,
        error_tolerance: FloatOrTolerance,
    ) -> PyResult<Py<PyList>> {
        match self.0.read() {
            Ok(ps) => {
                let peaks = ps.between(
                    low,
                    high,
                    <FloatOrTolerance as Into<PyTolerance>>::into(error_tolerance).0,
                );
                Python::with_gil(|py| {
                    let pl = PyList::empty_bound(py);
                    peaks
                        .into_iter()
                        .map(|p| -> PyResult<()> {
                            let p = CentroidPeakOrRef::Ref {
                                index: p.index() as usize,
                                peak_set: self.0.clone(),
                            };
                            let py_p = PyPeakRef(p);
                            pl.append(py_p.into_py(py))?;
                            Ok(())
                        })
                        .collect::<PyResult<()>>()?;
                    Ok(pl.into())
                })
            }
            Err(e) => Err(PyException::new_err(format!(
                "Failed to acquire lock for reading peak set: {e}"
            ))),
        }
    }

    #[pyo3(signature=(query, error_tolerance=FloatOrTolerance::Float(10.0)))]
    pub fn all_peaks_for(
        &self,
        query: f64,
        error_tolerance: FloatOrTolerance,
    ) -> PyResult<Py<PyList>> {
        match self.0.read() {
            Ok(ps) => {
                let peaks = ps.all_peaks_for(
                    query,
                    <FloatOrTolerance as Into<PyTolerance>>::into(error_tolerance).0,
                );
                Python::with_gil(|py| {
                    let pl = PyList::empty_bound(py);
                    peaks
                        .into_iter()
                        .map(|p| -> PyResult<()> {
                            let p = CentroidPeakOrRef::Ref {
                                index: p.index() as usize,
                                peak_set: self.0.clone(),
                            };
                            let py_p = PyPeakRef(p);
                            pl.append(py_p.into_py(py))?;
                            Ok(())
                        })
                        .collect::<PyResult<()>>()?;
                    Ok(pl.into())
                })
            }
            Err(e) => Err(PyException::new_err(format!(
                "Failed to acquire lock for reading peak set: {e}"
            ))),
        }
    }

    fn __getitem__(&self, i: Bound<PyAny>) -> PyResult<PyPeakRef> {
        if i.is_instance_of::<PySlice>() {
            Err(PyTypeError::new_err("Could not select indices by slice"))
        } else if i.is_instance_of::<PyLong>() {
            let ps = self.0.read().map_err(|e| {
                PyException::new_err(format!("Failed to acquire lock for reading peak set: {e}"))
            })?;
            let i: usize = i.extract()?;
            if i >= ps.len() {
                Err(PyIndexError::new_err(i))
            } else {
                let p = PyPeakRef(CentroidPeakOrRef::Ref {
                    index: i,
                    peak_set: self.0.clone(),
                });
                Ok(p)
            }
        } else {
            Err(PyTypeError::new_err("Could not select indices from input"))
        }
    }

    fn __len__(&self) -> PyResult<usize> {
        Ok(self.read_peak_set()?.len())
    }

    fn __repr__(&self) -> String {
        match self.__len__() {
            Ok(i) => format!("PyPeakSet({i} peaks)"),
            Err(e) => format!("PyPeakSet(error={e})"),
        }
    }
}

#[pyclass(module = "pymzpeaks", sequence, name = "PeakSet")]
#[derive(Debug, Clone)]
pub struct PyDeconvolutedPeakSet(PeakSetVec<PyDeconvolutedPeak, Mass>);

#[pymethods]
impl PyDeconvolutedPeakSet {
    #[new]
    pub fn py_new(peaks: Vec<PyDeconvolutedPeak>) -> Self {
        PyDeconvolutedPeakSet(PeakSetVec::new(peaks))
    }

    #[pyo3(signature=(query, error_tolerance=FloatOrTolerance::Float(10.0)))]
    pub fn has_peak(
        &self,
        query: f64,
        error_tolerance: FloatOrTolerance,
    ) -> Option<PyDeconvolutedPeak> {
        if let Some(peak) = self.0.has_peak(
            query,
            <FloatOrTolerance as Into<PyTolerance>>::into(error_tolerance).0,
        ) {
            Some(peak.clone())
        } else {
            None
        }
    }

    #[pyo3(signature=(query, error_tolerance=FloatOrTolerance::Float(10.0)))]
    pub fn all_peaks_for(
        &self,
        query: f64,
        error_tolerance: FloatOrTolerance,
    ) -> PyResult<Py<PyList>> {
        let peaks = self.0.all_peaks_for(
            query,
            <FloatOrTolerance as Into<PyTolerance>>::into(error_tolerance).0,
        );
        Python::with_gil(|py| {
            let pl = PyList::empty_bound(py);
            peaks
                .into_iter()
                .map(|p| -> PyResult<()> {
                    let py_p = (p.clone()).into_py(py);
                    pl.append(py_p)
                })
                .collect::<PyResult<()>>()?;
            Ok(pl.into())
        })
    }

    fn __getitem__(&self, i: Bound<PyAny>) -> PyResult<PyDeconvolutedPeak> {
        if i.is_instance_of::<PySlice>() {
            Err(PyTypeError::new_err("Could not select indices by slice"))
        } else if i.is_instance_of::<PyLong>() {
            let i: usize = i.extract()?;
            if i >= self.0.len() {
                Err(PyIndexError::new_err(i))
            } else {
                let p = self.0.get_item(i);
                Ok(p.clone())
            }
        } else {
            Err(PyTypeError::new_err("Could not select indices from input"))
        }
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    fn __repr__(&self) -> String {
        format!("PyDeconvolutedPeakSet({} peaks)", self.0.len())
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn pymzpeaks(_py: Python, m: Bound<PyModule>) -> PyResult<()> {
    m.add_class::<PyTolerance>()?;
    m.add_class::<PyCentroidPeak>()?;
    m.add_class::<PyPeakSet>()?;
    m.add_class::<PyDeconvolutedPeak>()?;
    m.add_class::<PyDeconvolutedPeakSet>()?;
    Ok(())
}
