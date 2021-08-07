#[derive(Default, Debug, Clone, Copy)]
/// The Mass To Charge Ratio (m/z) coordinate system
pub struct MZ {}

impl MZ {
    /// Access the m/z of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<MZ>>(inst: &T) -> f64 {
        CoordinateLike::<MZ>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy)]
/// The Mass coordinate system
pub struct Mass {}

impl Mass {
    /// Access the neutral mass of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<Mass>>(inst: &T) -> f64 {
        CoordinateLike::<Mass>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy)]
/// The Event Time coordinate system
pub struct Time {}
impl Time {
    /// Access the elapsed time of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<Time>>(inst: &T) -> f64 {
        CoordinateLike::<Time>::coordinate(inst)
    }
}

#[derive(Default, Debug, Clone, Copy)]
/// The Ion Mobility Time coordinate system
pub struct IonMobility {}
impl IonMobility {
    /// Access the ion mobility time unit of the coordinate type
    #[inline]
    pub fn coordinate<T: CoordinateLike<IonMobility>>(inst: &T) -> f64 {
        CoordinateLike::<IonMobility>::coordinate(inst)
    }
}

#[derive(Clone, Copy, Debug)]
pub enum CoordinateDimension {
    MZ(MZ),
    Mass(Mass),
    Time(Time),
    IonMobility(IonMobility),
}


/// Denote a type has a coordinate value on coordinate system `T`
pub trait CoordinateLike<T>: PartialOrd {

    /// The trait method for accessing the coordinate of the object on coordinate
    /// system `T`
    fn coordinate(&self) -> f64;
}

/// A named coordinate system membership for neutral mass
pub trait MassLocated: CoordinateLike<Mass> {
    #[inline]
    fn neutral_mass(&self) -> f64 {
        CoordinateLike::<Mass>::coordinate(self)
    }
}

/// A named coordinate system membership for m/z
pub trait MZLocated: CoordinateLike<MZ> {
    #[inline]
    fn mz(&self) -> f64 {
        CoordinateLike::<MZ>::coordinate(self)
    }
}

impl<T: CoordinateLike<Mass>> MassLocated for T {}
impl<T: CoordinateLike<MZ>> MZLocated for T {}

pub type IndexType = u32;


/// Indicate that an object may be indexed by coordinate system `T`
pub trait IndexedCoordinate<T>: CoordinateLike<T> {
    fn get_index(&self) -> IndexType;
    fn set_index(&mut self, index: IndexType);
}
