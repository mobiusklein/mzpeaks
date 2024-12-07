//! Measure accuracy error tolerance ranges
use std::{ops::{self, RangeInclusive}, fmt::Display, error::Error, str::FromStr};

#[cfg(feature = "serde")]
use serde::{Serialize, Deserialize};


/// A failure to parse a mass error tolerance quantity from a string
#[derive(Debug, PartialEq, Eq)]
pub enum ToleranceParsingError {
    /// The unit isn't empty, but not recognized
    UnknownUnit,
    /// The magnitude of the error tolerated couldn't be determined
    InvalidMagnitude
}

impl Display for ToleranceParsingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&format!("{:?}", self))
    }
}

impl Error for ToleranceParsingError {}

impl FromStr for Tolerance {
    type Err = ToleranceParsingError;

    /// Parse a string of the form "<magnitude:f64><unit:da|ppm?>"
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let n = s.len();
        if n <= 2 {
            return Err(ToleranceParsingError::InvalidMagnitude)
        }
        let s = s.to_lowercase();
        if s.ends_with("da") {
            if let Ok(magnitude) = s[0..n-2].parse::<f64>() {
                Ok(Self::Da(magnitude))
            } else {
                Err(ToleranceParsingError::InvalidMagnitude)
            }
        } else if s.ends_with("ppm") {
            if let Ok(magnitude) = s[0..n-3].parse::<f64>() {
                Ok(Self::PPM(magnitude))
            } else {
                Err(ToleranceParsingError::InvalidMagnitude)
            }
        } else {
            Err(ToleranceParsingError::UnknownUnit)
        }
    }
}


#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Tolerance {
    PPM(f64),
    Da(f64)
}

impl ToString for Tolerance {
    fn to_string(&self) -> String {
        match self {
            Self::Da(tol) => format!("{}Da", tol),
            Self::PPM(tol) => format!("{}PPM", tol),
        }
    }
}

impl Tolerance {

    /// The interval around `query` which is within this `Tolerance`
    /// instance's range.
    pub fn bounds(&self, query: f64) -> (f64, f64) {
        match self {
            Tolerance::PPM(tol) => {
                let width = query * *tol / 1e6;
                (query - width, query + width)
            }
            Tolerance::Da(tol) => {
                (query - *tol, query + *tol)
            }
        }
    }

    /// Compute the error between the two masses, in the appropriate units
    pub fn call(&self, query: f64, reference: f64) -> f64 {
        match self {
            Self::PPM(_tol) => {
                (query - reference) / reference * 1e6
            },
            Self::Da(_tol) => {
                query - reference
            }
        }
    }

    /// Return the numeric value of the error threshold in its units
    pub fn tol(&self) -> f64 {
        match self {
            Self::PPM(tol) => *tol,
            Self::Da(tol) => *tol
        }
    }

    /// Check if `query` is within the tolerated error interval around `reference`
    pub fn test(&self, query: f64, reference: f64) -> bool {
        let (lower_bound, upper_bound) = self.bounds(reference);
        query >= lower_bound && query <= upper_bound
    }

    /// Format the error between two masses with the appropriate units
    pub fn format_error(&self, query: f64, reference: f64) -> String {
        match self {
            Self::PPM(_tol) => {
                let magnitude = (query - reference) / reference * 1e6;
                format!("{}PPM", magnitude).to_string()
            },
            Self::Da(_tol) => {
                let magnitude = query - reference;
                format!("{}Da", magnitude).to_string()
            }
        }
    }

    pub fn as_range(&self, query: f64) -> RangeInclusive<f64> {
        let (low, hi) = self.bounds(query);
        RangeInclusive::new(low, hi)
    }

    pub fn as_interval(&self, query: f64) -> crate::coordinate::SimpleInterval<f64> {
        let bounds = self.bounds(query);
        bounds.into()
    }
}

/// Tolerance objects can by scaled up or down by a floating point value
impl ops::Mul<f64> for Tolerance {
    type Output = Tolerance;

    fn mul(self, rhs: f64) -> Self::Output {
        match self {
            Self::Da(val) => Self::Da(rhs * val),
            Self::PPM(val) => Self::PPM(rhs * val)
        }
    }
}

impl From<f64> for Tolerance {
    fn from(value: f64) -> Self {
        Self::PPM(value)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse() -> Result<(), ToleranceParsingError> {
        let mut tol: Tolerance = "10ppm".parse()?;
        assert_eq!(Tolerance::PPM(10.0), tol);

        assert_eq!(tol.to_string(), "10PPM");

        tol = "0.02Da".parse()?;
        assert_eq!(Tolerance::Da(0.02), tol);

        assert_eq!(tol.to_string(), "0.02Da");

        if let Err(e) = "0.02 da".parse::<Tolerance>() {
            match e {
                ToleranceParsingError::InvalidMagnitude => {},
                _ => panic!("Wrong error: {e}"),
            }
        }
        if let Err(e) = "1".parse::<Tolerance>() {
            match e {
                ToleranceParsingError::InvalidMagnitude => {},
                _ => panic!("Wrong error: {e}"),
            }
        }
        if let Err(e) = "0.02".parse::<Tolerance>() {
            match e {
                ToleranceParsingError::UnknownUnit => {},
                _ => panic!("Wrong error: {e}"),
            }
        }
        Ok(())
    }

    #[test]
    fn test_range() {
        let iv = Tolerance::Da(0.2).as_interval(50.0);
        assert_eq!(iv.start, 49.8);
        assert_eq!(iv.end, 50.2);

        let iv = (Tolerance::Da(0.2) * 10.0).as_interval(50.0);
        assert_eq!(iv.start, 48.0);
        assert_eq!(iv.end, 52.0);

        let iv = Tolerance::PPM(1.0).as_range(1e6);
        assert_eq!(*iv.start(), 1e6 - 1.0);
        assert_eq!(*iv.end(), 1e6 + 1.0);

        let iv = (Tolerance::PPM(1.0) * 10.0).as_range(1e6);
        assert_eq!(*iv.start(), 1e6 - 10.0);
        assert_eq!(*iv.end(), 1e6 + 10.0);
    }

    #[test]
    fn test_test() {
        assert!(Tolerance::PPM(10.0).test(100.0, 100.0));
        assert!(!Tolerance::PPM(10.0).test(100.0, 102.0));

        assert_eq!(
            Tolerance::PPM(10.0).format_error(100.0, 100.0),
            "0PPM"
        );
    }
}