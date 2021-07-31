#[derive(Debug, Clone, Copy)]
/// Represents a unit of mass accuracy, and provides the means
/// for measuring them.
pub enum MassErrorType {
    /// An absolute mass error is measured in Daltons: `(x - y)`
    Absolute,
    /// PPM, "Parts Per Million", is relative to the masses being compared: `(x - y) / y`
    PPM,
}

impl MassErrorType {
    #[inline]
    /// The minimum value which would match `query` with `tolerance` units of error
    pub fn lower_bound(&self, query: f64, tolerance: f64) -> f64 {
        match self {
            Self::Absolute => query - tolerance,
            Self::PPM => query - (query * tolerance * 1e-6),
        }
    }

    #[inline]
    /// Compute the actual error between `query` and `alt`
    pub fn call(&self, query: f64, alt: f64) -> f64 {
        match self {
            Self::Absolute => query - alt,
            Self::PPM => (query - alt) / alt,
        }
    }

    #[inline]
    /// The maximum value which would match `query` with `tolerance` units of error
    pub fn upper_bound(&self, query: f64, tolerance: f64) -> f64 {
        self.lower_bound(query, -tolerance)
    }
}
