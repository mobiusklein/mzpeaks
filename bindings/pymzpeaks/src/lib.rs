use std::ops::{Index, IndexMut, Deref};
use std::sync::{Arc, RwLock};

use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::pyclass::CompareOp;
use pyo3::types::{PySlice, PyLong, PyFloat, PyString, PyList};
use pyo3::prelude::*;

use mzpeaks::{
    CoordinateLike, IndexedCoordinate, IntensityMeasurement, CentroidPeak,
    KnownCharge, DeconvolutedCentroidLike, DeconvolutedPeak, MassLocated
};
use mzpeaks::coordinate::{MZ, Mass, MZLocated};
use mzpeaks::Tolerance;
use mzpeaks::peak_set::{PeakCollection, PeakSetVec};


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


#[pyclass(module="pymzpeaks", name="Tolerance")]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PyTolerance(Tolerance);

#[pymethods]
impl PyTolerance {
    #[new]
    fn new(value: &PyAny) -> PyResult<Self> {
        if value.is_instance_of::<PyFloat>()? {
            Ok(Self(Tolerance::PPM(value.extract::<f64>()?)))
        }
        else if value.is_instance_of::<PyString>()? {
            match value.extract::<&str>()?.parse::<Tolerance>() {
                Ok(tol) => Ok(Self(tol)),
                Err(err) => {
                    Err(PyValueError::new_err(err.to_string()))
                }
            }
        } else {
            Err(PyTypeError::new_err(format!("Could not convert {} into Tolerance instance", value.get_type().getattr("__name__")?.extract::<&str>()?)))
        }
    }
}

impl std::str::FromStr for PyTolerance {
    type Err = <Tolerance as std::str::FromStr>::Err;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match <Tolerance as std::str::FromStr>::from_str(s) {
            Ok(tol) => Ok(Self(tol)),
            Err(err) => Err(err)
        }
    }
}

impl ToString for PyTolerance {
    fn to_string(&self) -> String {
        <Tolerance as ToString>::to_string( &self.0 )
    }
}

impl Default for PyTolerance {
    fn default() -> Self {
        Self(Tolerance::PPM(20.0))
    }
}


#[pyclass(module="pymzpeaks", name="CentroidPeak")]
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct PyCentroidPeak(CentroidPeak);


#[pymethods]
impl PyCentroidPeak {
    #[new]
    pub fn new(mz: f64, intensity: f32, index: u32) -> Self {
        (CentroidPeak {mz, intensity, index}).into()
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
        return self.0.index
    }

    fn __repr__(&self) -> String {
        format!("CentroidPeak({:0.4}, {:0.4}, {})", self.mz(), self.intensity(), self.index())
    }

    fn __richcmp__(&self, other: PyRef<PyCentroidPeak>, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(self.0 == other.0),
            CompareOp::Ne => Ok(self.0 != other.0),
            CompareOp::Le => Ok(self.0 <= other.0),
            CompareOp::Ge => Ok(self.0 >= other.0),
            CompareOp::Lt => Ok(self.0 < other.0),
            CompareOp::Gt => Ok(self.0 > other.0)
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
    fn intensity(&self) -> f32 {
        self.0.intensity
    }
}


#[pyclass(module="pymzpeaks", name="DeconvolutedPeak")]
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct PyDeconvolutedPeak(DeconvolutedPeak);

#[pymethods]
impl PyDeconvolutedPeak {
    #[new]
    pub fn new(neutral_mass: f64, intensity: f32, charge: i32, index: u32) -> Self {
        (DeconvolutedPeak {neutral_mass, intensity, index, charge}).into()
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
        return self.0.index
    }

    #[getter]
    fn charge(&self) -> i32 {
        self.0.charge
    }

    fn __repr__(&self) -> String {
        format!(
            "DeconvolutedPeak({:0.4}, {:0.4}, {}, {})",
            self.neutral_mass(), self.intensity(), self.charge(), self.index()
        )
    }

    fn __richcmp__(&self, other: PyRef<PyDeconvolutedPeak>, op: CompareOp) -> PyResult<bool> {
        match op {
            CompareOp::Eq => Ok(self.0 == other.0),
            CompareOp::Ne => Ok(self.0 != other.0),
            CompareOp::Le => Ok(self.0 <= other.0),
            CompareOp::Ge => Ok(self.0 >= other.0),
            CompareOp::Lt => Ok(self.0 < other.0),
            CompareOp::Gt => Ok(self.0 > other.0)
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

impl KnownCharge for PyDeconvolutedPeak {
    fn charge(&self) -> i32 {
        self.0.charge
    }
}


#[pyclass]
#[derive(Debug, Clone)]
pub struct PyPeakRef(usize, Arc<RwLock<PeakSetVec<PyCentroidPeak, MZ>>>);

#[pymethods]
impl PyPeakRef {
    #[getter]
    fn mz(&self) -> PyResult<f64> {
        match self.1.as_ref().read() {
            Ok(peaks) => {
                Ok(peaks.deref().index(self.0).mz())
            },
            Err(err) => {
                Err(PyValueError::new_err(format!("Failed to get read access to peak set: {}", err)))
            }
        }
    }

    #[setter]
    fn mz_setter(&self, val: f64) -> PyResult<()>{
        match self.1.as_ref().write() {
            Ok(mut peaks) => {
                peaks.index_mut(self.0).0.mz = val;
                Ok(())
            },
            Err(err) => {
                Err(PyValueError::new_err(format!("Failed to get read access to peak set: {}", err)))
            }
        }
    }

    #[getter]
    fn intensity(&self) -> PyResult<f32> {
        match self.1.as_ref().read() {
            Ok(peaks) => {
                Ok(peaks.deref().index(self.0).intensity())
            },
            Err(err) => {
                Err(PyValueError::new_err(format!("Failed to get read access to peak set: {}", err)))
            }
        }
    }

    #[setter]
    fn intensity_setter(&self, val: f32) -> PyResult<()>{
        match self.1.as_ref().write() {
            Ok(mut peaks) => {
                peaks.index_mut(self.0).0.intensity = val;
                Ok(())
            },
            Err(err) => {
                Err(PyValueError::new_err(format!("Failed to get read access to peak set: {}", err)))
            }
        }
    }

    #[getter]
    fn index(&self) -> PyResult<u32> {
        match self.1.as_ref().read() {
            Ok(peaks) => {
                Ok(peaks.deref().index(self.0).index())
            },
            Err(err) => {
                Err(PyValueError::new_err(format!("Failed to get read access to peak set: {}", err)))
            }
        }
    }

    #[setter]
    fn index_setter(&self, val: u32) -> PyResult<()>{
        match self.1.as_ref().write() {
            Ok(mut peaks) => {
                peaks.index_mut(self.0).0.index = val;
                Ok(())
            },
            Err(err) => {
                Err(PyValueError::new_err(format!("Failed to get read access to peak set: {}", err)))
            }
        }
    }
}


#[pyclass(module="pymzpeaks", sequence, name="PeakSet")]
#[derive(Debug, Clone)]
pub struct PyPeakSet(PeakSetVec<PyCentroidPeak, MZ>);

#[pymethods]
impl PyPeakSet {

    #[new]
    pub fn py_new(peaks: Vec<PyCentroidPeak>) -> Self {
        PyPeakSet(PeakSetVec::new(peaks))
    }

    #[pyo3(signature=(query, error_tolerance=FloatOrTolerance::Float(10.0)))]
    pub fn has_peak(&self, query: f64, error_tolerance: FloatOrTolerance) -> Option<PyCentroidPeak> {
        if let Some(peak) = self.0.has_peak(query, <FloatOrTolerance as Into<PyTolerance>>::into(error_tolerance).0) {
            Some(peak.clone())
        } else {
            None
        }
    }

    #[pyo3(signature=(query, error_tolerance=FloatOrTolerance::Float(10.0)))]
    pub fn all_peaks_for(&self, query: f64, error_tolerance: FloatOrTolerance) -> PyResult<Py<PyList>> {
        let peaks = self.0.all_peaks_for(query, <FloatOrTolerance as Into<PyTolerance>>::into(error_tolerance).0);
        Python::with_gil(|py| {
            let pl = PyList::empty(py);
            peaks.into_iter().map(|p| -> PyResult<()> {
                let py_p = (p.clone()).into_py(py);
                pl.append(py_p)
            }).collect::<PyResult<()>>()?;
            Ok(pl.into())
        })
    }

    fn __getitem__(&self,  i: &PyAny) -> PyResult<PyCentroidPeak> {
        if i.is_instance_of::<PySlice>()? {
            Err(PyTypeError::new_err("Could not select indices by slice"))
        } else if i.is_instance_of::<PyLong>()? {
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
        format!("PyPeakSet({} peaks)", self.0.len())
    }
}


#[pyclass(module="pymzpeaks", sequence, name="PeakSet")]
#[derive(Debug, Clone)]
pub struct PyDeconvolutedPeakSet(PeakSetVec<PyDeconvolutedPeak, Mass>);

#[pymethods]
impl PyDeconvolutedPeakSet {

    #[new]
    pub fn py_new(peaks: Vec<PyDeconvolutedPeak>) -> Self {
        PyDeconvolutedPeakSet(PeakSetVec::new(peaks))
    }

    #[pyo3(signature=(query, error_tolerance=FloatOrTolerance::Float(10.0)))]
    pub fn has_peak(&self, query: f64, error_tolerance: FloatOrTolerance) -> Option<PyDeconvolutedPeak> {
        if let Some(peak) = self.0.has_peak(query, <FloatOrTolerance as Into<PyTolerance>>::into(error_tolerance).0) {
            Some(peak.clone())
        } else {
            None
        }
    }

    #[pyo3(signature=(query, error_tolerance=FloatOrTolerance::Float(10.0)))]
    pub fn all_peaks_for(&self, query: f64, error_tolerance: FloatOrTolerance) -> PyResult<Py<PyList>> {
        let peaks = self.0.all_peaks_for(query, <FloatOrTolerance as Into<PyTolerance>>::into(error_tolerance).0);
        Python::with_gil(|py| {
            let pl = PyList::empty(py);
            peaks.into_iter().map(|p| -> PyResult<()> {
                let py_p = (p.clone()).into_py(py);
                pl.append(py_p)
            }).collect::<PyResult<()>>()?;
            Ok(pl.into())
        })
    }

    fn __getitem__(&self,  i: &PyAny) -> PyResult<PyDeconvolutedPeak> {
        if i.is_instance_of::<PySlice>()? {
            Err(PyTypeError::new_err("Could not select indices by slice"))
        } else if i.is_instance_of::<PyLong>()? {
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
fn pymzpeaks(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyTolerance>()?;
    m.add_class::<PyCentroidPeak>()?;
    m.add_class::<PyPeakSet>()?;
    m.add_class::<PyDeconvolutedPeak>()?;
    m.add_class::<PyDeconvolutedPeakSet>()?;
    Ok(())
}