#![allow(unused)]
use std::{marker::PhantomData, ops::Sub};
use num_traits::Num;

use super::range::{CoordinateRange, SimpleInterval, Span1D};

/** An inclusive interval over two dimensions
*/
pub trait Span2D {
    type DimType1: PartialOrd + Copy;
    type DimType2: PartialOrd + Copy;

    fn start(&self) -> (Self::DimType1, Self::DimType2);
    fn end(&self) -> (Self::DimType1, Self::DimType2);

    fn contains(&self, i: (&Self::DimType1, &Self::DimType2)) -> bool {
        let (x, y) = i;
        let (sx, sy) = self.start();
        let (ex, ey) = self.end();
        SimpleInterval::new(sx, ex).contains(x) && SimpleInterval::new(sy, ey).contains(y)
    }

    fn contains_interval(&self, other: impl Span2D<DimType1 = Self::DimType1, DimType2 = Self::DimType2>) -> bool {
        self.start() <= other.start() && self.end() >= other.end()
    }

    fn overlaps<T: Span2D<DimType1 = Self::DimType1, DimType2 = Self::DimType2>>(
        &self,
        interval: &T,
    ) -> bool {
        self.end() >= interval.start() && interval.end() >= self.start()
    }
}

impl<T: Span2D> Span2D for &T {
    type DimType1 = T::DimType1;

    type DimType2 = T::DimType2;

    fn start(&self) -> (Self::DimType1, Self::DimType2) {
        (*self).start()
    }

    fn end(&self) -> (Self::DimType1, Self::DimType2) {
        (*self).end()
    }
}

/// A basic [`Span2D`] implementation
#[derive(Debug, Default, Clone, Copy, PartialEq, PartialOrd)]
pub struct BoundingBox<V1: PartialOrd + Num, V2: PartialOrd + Num> {
    pub start: (V1, V2),
    pub end: (V1, V2),
}

impl<V1: PartialOrd + Copy + Num, V2: PartialOrd + Copy + Num> Span2D for BoundingBox<V1, V2> {
    type DimType1 = V1;

    type DimType2 = V2;

    fn start(&self) -> (Self::DimType1, Self::DimType2) {
        self.start
    }

    fn end(&self) -> (Self::DimType1, Self::DimType2) {
        self.end
    }
}

impl<V1: PartialOrd + Copy + Num, V2: PartialOrd + Copy + Num> BoundingBox<V1, V2> {
    pub fn new(start: (V1, V2), end: (V1, V2)) -> Self {
        Self { start, end }
    }

    pub fn centroid(&self) -> (V1, V2) {
        let dim0 = self.end.0 - self.start.0;
        let dim1 = self.end.1 - self.start.1;
        (dim0, dim1)
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, PartialOrd)]
pub struct CoordinateBox<X, Y> {
    dim1: CoordinateRange<X>,
    dim2: CoordinateRange<Y>,
}

impl<X, Y> Span2D for CoordinateBox<X, Y> {
    type DimType1 = <CoordinateRange<X> as Span1D>::DimType;

    type DimType2 = <CoordinateRange<Y> as Span1D>::DimType;

    fn start(&self) -> (Self::DimType1, Self::DimType2) {
        (self.dim1.start(), self.dim2.start())
    }

    fn end(&self) -> (Self::DimType1, Self::DimType2) {
        (self.dim1.end(), self.dim2.end())
    }
}

impl<X, Y> CoordinateBox<X, Y> {
    pub fn new(dim1: CoordinateRange<X>, dim2: CoordinateRange<Y>) -> Self {
        Self { dim1, dim2 }
    }
}

#[derive(Debug, Clone, Default)]
pub struct QuadTreeNode<
    V1: PartialOrd + Copy,
    V2: PartialOrd + Copy,
    T: Span2D<DimType1 = V1, DimType2 = V2>,
> {
    pub start: (V1, V2),
    pub end: (V1, V2),
    pub center: (V1, V2),
    pub level: u32,
    pub members: Vec<T>,
    pub parent: Option<usize>,
    pub children: Option<[usize; 4]>,
    _t: PhantomData<T>,
}

impl<V1: PartialOrd + Copy + Num, V2: PartialOrd + Copy + Num, T: Span2D<DimType1 = V1, DimType2 = V2>> From<BoundingBox<V1, V2>> for QuadTreeNode<V1, V2, T> {
    fn from(value: BoundingBox<V1, V2>) -> Self {
        Self::empty(value)
    }
}

impl<V1: PartialOrd + Copy + Num, V2: PartialOrd + Copy + Num, T: Span2D<DimType1 = V1, DimType2 = V2>>
    QuadTreeNode<V1, V2, T>
{
    pub fn empty(bbox: BoundingBox<V1, V2>) -> Self {
        Self {
            start: bbox.start(),
            end: bbox.end(),
            center: bbox.centroid(),
            level: 0,
            members: Vec::new(),
            parent: None,
            children: None,
            _t: PhantomData
        }
    }

    pub fn is_leaf(&self) -> bool {
        !self.children.is_some()
    }

    pub fn split(&mut self) -> [Self; 4] {
        let mut upper_left: Self =
            BoundingBox::new((self.start.0, self.start.1), (self.center.0, self.center.1)).into();
        let mut lower_right: Self =
            BoundingBox::new((self.center.0, self.center.1), (self.end.0, self.end.1)).into();
        let mut lower_left: Self =
            BoundingBox::new((self.start.0, self.center.1), (self.center.0, self.end.1)).into();
        let mut upper_right: Self =
            BoundingBox::new((self.center.0, self.start.1), (self.end.0, self.center.1)).into();

        upper_left.level = self.level + 1;
        upper_right.level = self.level + 1;
        lower_left.level = self.level + 1;
        lower_right.level = self.level + 1;

        let mut leaves = [upper_left, upper_right, lower_left, lower_right];

        let mut kept = Vec::new();
        for point in self.members.drain(..) {
            let mut container: Option<usize> = None;
            for (i, child) in leaves.iter().enumerate() {
                if child.contains_interval(&point) {
                    container = Some(i);
                    break;
                }
            }
            match container {
                Some(0) | Some(1) | Some(2) | Some(3) => {
                    leaves.get_mut(container.unwrap()).unwrap().members.push(point);
                },
                _ => kept.push(point),
            }
        }
        self.members = kept;
        leaves
    }
}

impl<V1: PartialOrd + Copy, V2: PartialOrd + Copy, T: Span2D<DimType1 = V1, DimType2 = V2>> Span2D
    for QuadTreeNode<V1, V2, T>
{
    type DimType1 = V1;

    type DimType2 = V2;

    fn start(&self) -> (Self::DimType1, Self::DimType2) {
        self.start
    }

    fn end(&self) -> (Self::DimType1, Self::DimType2) {
        self.end
    }
}

#[derive(Debug, Clone, Default)]
pub struct QuadTree<
    V1: PartialOrd + Copy,
    V2: PartialOrd + Copy,
    T: Span2D<DimType1 = V1, DimType2 = V2>,
> {
    pub nodes: Vec<QuadTreeNode<V1, V2, T>>,
}
