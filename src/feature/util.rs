use std::cmp::Ordering;

#[derive(PartialEq)]
pub(crate) struct NonNan(f64);

impl NonNan {
    pub fn new(val: f64) -> Option<NonNan> {
        if val.is_nan() {
            None
        } else {
            Some(NonNan(val))
        }
    }
}

impl Eq for NonNan {}

impl PartialOrd for NonNan {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for NonNan {
    fn cmp(&self, other: &NonNan) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

pub const EMPTY_X: &[f64] = &[];
pub const EMPTY_Y: &[f64] = &[];
pub const EMPTY_Z: &[f32] = &[];
