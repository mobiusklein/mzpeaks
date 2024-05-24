#![allow(unused)]
use std::marker::PhantomData;

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

    fn overlaps<T: Span2D<DimType1 = Self::DimType1, DimType2 = Self::DimType2>>(
        &self,
        interval: &T,
    ) -> bool {
        self.end() >= interval.start() && interval.end() >= self.start()
    }
}

/// A basic [`Span2D`] implementation
#[derive(Debug, Default, Clone, Copy, PartialEq, PartialOrd)]
pub struct BoundingBox<V1: PartialOrd, V2: PartialOrd> {
    pub start: (V1, V2),
    pub end: (V1, V2),
}

impl<V1: PartialOrd + Copy, V2: PartialOrd + Copy> Span2D for BoundingBox<V1, V2> {
    type DimType1 = V1;

    type DimType2 = V2;

    fn start(&self) -> (Self::DimType1, Self::DimType2) {
        self.start
    }

    fn end(&self) -> (Self::DimType1, Self::DimType2) {
        self.end
    }
}

impl<V1: PartialOrd, V2: PartialOrd> BoundingBox<V1, V2> {
    pub fn new(start: (V1, V2), end: (V1, V2)) -> Self {
        Self { start, end }
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
pub struct QuadTreeNode<V1: PartialOrd + Copy, V2: PartialOrd + Copy, T: Span2D<DimType1 = V1, DimType2 = V2>> {
    pub start: (V1, V2),
    pub end: (V1, V2),
    pub center: (V1, V2),
    pub level: u32,
    pub members: Vec<T>,
    pub parent: Option<usize>,
    pub children: [Option<usize>; 4],
    _t: PhantomData<T>
}

impl<V1: PartialOrd + Copy, V2: PartialOrd + Copy, T: Span2D<DimType1 = V1, DimType2 = V2>> QuadTreeNode<V1, V2, T> {
    pub fn is_leaf(&self) -> bool {
        !self.children.iter().any(|s| s.is_some())
    }
}

impl<V1: PartialOrd + Copy, V2: PartialOrd + Copy, T: Span2D<DimType1 = V1, DimType2 = V2>> Span2D for QuadTreeNode<V1, V2, T> {
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
pub struct QuadTree<V1: PartialOrd + Copy, V2: PartialOrd + Copy, T: Span2D<DimType1 = V1, DimType2 = V2>> {
    pub nodes: Vec<QuadTreeNode<V1, V2, T>>
}