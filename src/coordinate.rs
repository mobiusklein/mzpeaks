#[derive(Default, Debug, Clone, Copy)]
pub struct MZ {}

impl MZ {
    pub fn coordinate<T: CoordinateLike<MZ>>(inst: &T) -> f64 {
        CoordinateLike::<MZ>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy)]
pub struct Mass {}

impl Mass {
    pub fn coordinate<T: CoordinateLike<Mass>>(inst: &T) -> f64 {
        CoordinateLike::<Mass>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy)]
pub struct Time {}
impl Time {
    pub fn coordinate<T: CoordinateLike<Time>>(inst: &T) -> f64 {
        CoordinateLike::<Time>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy)]
pub struct IonMobility {}
impl IonMobility {
    pub fn coordinate<T: CoordinateLike<IonMobility>>(inst: &T) -> f64 {
        CoordinateLike::<IonMobility>::coordinate(inst)
    }
}

pub enum CoordinateDimension {
    MZ(MZ),
    Mass(Mass),
    Time(Time),
    IonMobility(IonMobility),
}

pub trait CoordinateLike<T>: PartialOrd {
    fn coordinate(&self) -> f64;
}

pub trait MassLocated: CoordinateLike<Mass> {
    fn neutral_mass(&self) -> f64 {
        CoordinateLike::<Mass>::coordinate(self)
    }
}

pub trait MZLocated: CoordinateLike<MZ> {
    fn mz(&self) -> f64 {
        CoordinateLike::<MZ>::coordinate(self)
    }
}

impl<T: CoordinateLike<Mass>> MassLocated for T {}

pub type IndexType = u32;

pub trait IndexedCoordinate<T>: CoordinateLike<T> {
    fn get_index(&self) -> IndexType;
    fn set_index(&mut self, index: IndexType);
}
