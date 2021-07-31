#[derive(Default, Debug, Clone, Copy)]
pub struct MZ {}
#[derive(Default, Debug, Clone, Copy)]
pub struct Mass {}
#[derive(Default, Debug, Clone, Copy)]
pub struct Time {}
#[derive(Default, Debug, Clone, Copy)]
pub struct IonMobility {}

pub enum CoordinateDimension {
    MZ(MZ),
    Mass(Mass),
    Time(Time),
    IonMobility(IonMobility),
}

pub trait CoordinateLike<T>: PartialOrd {
    fn coordinate(&self) -> f64;
}

pub type IndexType = u32;

pub trait IndexedCoordinate<T>: CoordinateLike<T> {
    fn get_index(&self) -> IndexType;
    fn set_index(&mut self, index: IndexType);
}
