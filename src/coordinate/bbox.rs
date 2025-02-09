#![allow(unused)]
use num_traits::Num;
use std::{collections::VecDeque, iter::FusedIterator, marker::PhantomData, ops::Sub};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use super::{
    range::{CoordinateRange, SimpleInterval, Span1D},
    HasProximity,
};

/** An inclusive interval over two dimensions
*/
pub trait Span2D {
    type DimType1: HasProximity;
    type DimType2: HasProximity;

    fn start(&self) -> (Self::DimType1, Self::DimType2);
    fn end(&self) -> (Self::DimType1, Self::DimType2);

    fn is_close(
        &self,
        other: impl Span2D<DimType1 = Self::DimType1, DimType2 = Self::DimType2>,
    ) -> bool {
        let s1 = self.start();
        let s2 = other.start();
        if !(s1.0.is_close(&s2.0) && s1.1.is_close(&s2.1)) {
            return false;
        }
        let s1 = self.end();
        let s2 = other.end();
        if !(s1.0.is_close(&s2.0) && s1.1.is_close(&s2.1)) {
            return false;
        }
        true
    }

    fn contains(&self, i: (&Self::DimType1, &Self::DimType2)) -> bool {
        let (x, y) = i;
        let (sx, sy) = self.start();
        let (ex, ey) = self.end();
        SimpleInterval::new(sx, ex).contains(x) && SimpleInterval::new(sy, ey).contains(y)
    }

    fn as_bounding_box(&self) -> BoundingBox<Self::DimType1, Self::DimType2> {
        BoundingBox::new(self.start(), self.end())
    }

    fn contains_interval(
        &self,
        other: impl Span2D<DimType1 = Self::DimType1, DimType2 = Self::DimType2>,
    ) -> bool {
        (self.start() <= other.start() && self.end() >= other.end()) || self.is_close(other)
    }

    fn overlaps<T: Span2D<DimType1 = Self::DimType1, DimType2 = Self::DimType2>>(
        &self,
        interval: &T,
    ) -> bool {
        (self.end() >= interval.start() && interval.end() >= self.start())
            || self.is_close(interval)
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
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BoundingBox<V1: HasProximity, V2: HasProximity> {
    pub start: (V1, V2),
    pub end: (V1, V2),
}

impl<V1: HasProximity, V2: HasProximity> Span2D for BoundingBox<V1, V2> {
    type DimType1 = V1;

    type DimType2 = V2;

    fn start(&self) -> (Self::DimType1, Self::DimType2) {
        self.start
    }

    fn end(&self) -> (Self::DimType1, Self::DimType2) {
        self.end
    }
}

impl<V1: HasProximity, V2: HasProximity> BoundingBox<V1, V2> {
    pub fn new(start: (V1, V2), end: (V1, V2)) -> Self {
        Self { start, end }
    }

    pub fn combine_over<
        'a,
        T: Span2D<DimType1 = V1, DimType2 = V2> + 'a,
        I: Iterator<Item = &'a T>,
    >(
        mut iter: I,
    ) -> Option<Self> {
        if let Some(first) = iter.next() {
            let bbox = first.as_bounding_box();
            let bbox = iter.fold(bbox, |acc, item| acc.combine(&item));
            Some(bbox)
        } else {
            None
        }
    }

    pub fn combine(&self, other: &impl Span2D<DimType1 = V1, DimType2 = V2>) -> Self {
        let mut bbox = *self;
        let start = other.start();
        if start.0 < bbox.start.0 {
            bbox.start.0 = start.0;
        }
        if start.1 < bbox.start.1 {
            bbox.start.1 = start.1;
        }
        let end = other.end();
        if end.0 > bbox.end.0 {
            bbox.end.0 = end.0;
        }
        if end.1 > bbox.end.1 {
            bbox.end.1 = end.1;
        }
        bbox
    }

    pub fn centroid(&self) -> (V1, V2)
    where
        V1: Num,
        V2: Num,
    {
        let dim0 = (self.end.0 - self.start.0) / (V1::one() + V1::one()) + self.start.0;
        let dim1 = (self.end.1 - self.start.1) / (V2::one() + V2::one()) + self.start.1;
        (dim0, dim1)
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
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
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct QuadTreeNode<V1: HasProximity, V2: HasProximity, T: Span2D<DimType1 = V1, DimType2 = V2>>
{
    pub start: (V1, V2),
    pub end: (V1, V2),
    pub center: (V1, V2),
    pub level: u32,
    pub members: Vec<T>,
    pub parent: Option<usize>,
    pub children: Option<[usize; 4]>,
    _t: PhantomData<T>,
}

impl<V1: HasProximity, V2: HasProximity, T: Span2D<DimType1 = V1, DimType2 = V2>> IntoIterator
    for QuadTreeNode<V1, V2, T>
{
    type Item = T;

    type IntoIter = std::vec::IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.members.into_iter()
    }
}

impl<'a, V1: HasProximity, V2: HasProximity, T: Span2D<DimType1 = V1, DimType2 = V2>> IntoIterator
    for &'a QuadTreeNode<V1, V2, T>
{
    type Item = &'a T;

    type IntoIter = std::slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.members.iter()
    }
}

impl<V1: HasProximity + Num, V2: HasProximity + Num, T: Span2D<DimType1 = V1, DimType2 = V2>>
    From<BoundingBox<V1, V2>> for QuadTreeNode<V1, V2, T>
{
    fn from(value: BoundingBox<V1, V2>) -> Self {
        Self::empty(value)
    }
}

impl<V1: HasProximity + Num, V2: HasProximity + Num, T: Span2D<DimType1 = V1, DimType2 = V2>>
    QuadTreeNode<V1, V2, T>
{
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.members.iter()
    }

    pub fn empty(bbox: BoundingBox<V1, V2>) -> Self {
        Self {
            start: bbox.start(),
            end: bbox.end(),
            center: bbox.centroid(),
            level: 0,
            members: Vec::new(),
            parent: None,
            children: None,
            _t: PhantomData,
        }
    }

    pub fn is_leaf(&self) -> bool {
        self.children.is_none()
    }

    fn update_bounds(&mut self, children: &[BoundingBox<V1, V2>]) -> bool {
        let bbox = if !self.members.is_empty() {
            let bbox = self.members.first().unwrap().as_bounding_box();
            let bbox = self
                .members
                .iter()
                .skip(1)
                .fold(bbox, |bbox: BoundingBox<V1, V2>, member| {
                    bbox.combine(member)
                });
            bbox.combine(&self.as_bounding_box())
        } else {
            self.as_bounding_box()
        };

        let bbox = children
            .iter()
            .fold(bbox, |bbox: BoundingBox<V1, V2>, node| bbox.combine(node));
        let change_start = if self.start != bbox.start {
            self.start = bbox.start;
            true
        } else {
            false
        };
        let change_end = if self.end != bbox.end {
            self.end = bbox.end;
            true
        } else {
            false
        };
        change_start || change_end
    }

    fn split(&mut self) -> [Self; 4] {
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
                    leaves
                        .get_mut(container.unwrap())
                        .unwrap()
                        .members
                        .push(point);
                }
                _ => kept.push(point),
            }
        }
        self.members = kept;
        leaves
    }
}

impl<V1: HasProximity, V2: HasProximity, T: Span2D<DimType1 = V1, DimType2 = V2>> Span2D
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

const SPLIT_THRESHOLD: usize = 8;

#[derive(Debug, Clone, Default)]
pub struct QuadTree<
    V1: HasProximity + Num,
    V2: HasProximity + Num,
    T: Span2D<DimType1 = V1, DimType2 = V2>,
> {
    pub nodes: Vec<QuadTreeNode<V1, V2, T>>,
}

impl<V1: HasProximity + Num, V2: HasProximity + Num, T: Span2D<DimType1 = V1, DimType2 = V2>>
    FromIterator<T> for QuadTree<V1, V2, T>
{
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let items = iter.into_iter().collect();
        Self::build(items)
    }
}

type ContainsBBoxIter<'a, V1, V2, T> = QueryIter<
        'a,
        V1,
        V2,
        T,
        BoundingBox<V1, V2>,
        ContainsPredicate<V1, V2, T, BoundingBox<V1, V2>>,
    >;

impl<
        'a,
        V1: HasProximity + Num,
        V2: HasProximity + Num,
        T: Span2D<DimType1 = V1, DimType2 = V2>,
    > QuadTree<V1, V2, T>
{
    pub fn empty() -> Self {
        Self::new(Vec::new())
    }

    pub fn build(members: Vec<T>) -> Self {
        if members.is_empty() {
            Self::new(Vec::new())
        } else {
            let root_box = BoundingBox::combine_over(members.iter()).unwrap();
            let root = QuadTreeNode::empty(root_box);
            let mut this = Self::new(vec![root]);

            for member in members {
                this.insert(member);
            }

            this
        }
    }

    pub fn new(nodes: Vec<QuadTreeNode<V1, V2, T>>) -> Self {
        Self { nodes }
    }

    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    fn root(&self) -> &QuadTreeNode<V1, V2, T> {
        &self.nodes[0]
    }

    pub fn contains_point(&self, value: (&V1, &V2)) -> Vec<&T> {
        let mut results: Vec<&T> = Vec::new();

        if self.nodes.is_empty() {
            return results;
        }
        if !self.nodes[0].contains(value) {
            return results;
        }

        let mut queue = VecDeque::new();
        queue.push_back(&self.nodes[0]);

        while !queue.is_empty() {
            let node = queue.pop_front().unwrap();
            for i in node.members.iter() {
                if i.contains(value) {
                    results.push(i);
                }
            }
            queue.extend(node.children.into_iter().flatten().map(|i| &self.nodes[i]));
        }

        results
    }

    pub fn contains_iter(
        &'a self,
        value: (V1, V2),
    ) -> ContainsBBoxIter<'a, V1, V2, T> {
        let v = BoundingBox::new(value, value);
        QueryIter::new(self, v)
    }

    pub fn insert(&mut self, item: T) {
        if self.is_empty() {
            let mut node = QuadTreeNode::from(item.as_bounding_box());
            node.members.push(item);
            self.nodes.push(node);
            return;
        }

        let deepest_node = self
            .nodes_overlaps(&item)
            .into_iter()
            .map(|(i, node)| i)
            .max().unwrap_or(0);

        let (needs_splitting, mut children_idxes) = {
            let node = self.nodes.get_mut(deepest_node).unwrap();
            node.members.push(item);
            (
                node.members.len() > SPLIT_THRESHOLD && node.is_leaf(),
                node.children,
            )
        };

        if needs_splitting {
            let mut children = self.nodes.get_mut(deepest_node).unwrap().split();
            let mut idxes = [0; 4];
            for (j, mut child) in children.into_iter().enumerate() {
                child.parent = Some(deepest_node);
                idxes[j] = self.nodes.len();
                self.nodes.push(child);
            }
            self.nodes.get_mut(deepest_node).unwrap().children = Some(idxes);
        }

        let mut children: Vec<_> = children_idxes
            .into_iter()
            .flatten()
            .flat_map(|i| self.nodes.get(i))
            .map(|c| c.as_bounding_box())
            .collect();
        let node = self.nodes.get_mut(deepest_node).unwrap();

        let mut changed = node.update_bounds(&children);
        let mut parent_idx = node.parent;

        while changed && parent_idx.is_some() {
            let node = self.nodes.get(parent_idx.unwrap()).unwrap();
            children_idxes = node.children;
            let next_parent_idx = node.parent;

            children.clear();
            children.extend(
                children_idxes
                    .into_iter()
                    .flatten()
                    .flat_map(|i| self.nodes.get(i))
                    .map(|c| c.as_bounding_box()),
            );

            let node = self.nodes.get_mut(parent_idx.unwrap()).unwrap();
            changed = node.update_bounds(&children);
            parent_idx = next_parent_idx;
        }
    }

    pub fn overlaps_iter<Q: Span2D<DimType1 = V1, DimType2 = V2> + 'a>(
        &'a self,
        query: Q,
    ) -> QueryIter<'a, V1, V2, T, Q, OverlapPredicate<V1, V2, T, Q>> {
        QueryIter::new(self, query)
    }

    pub fn iter(&self) -> PreorderIter<'_, V1, V2, T> {
        PreorderIter::new(self)
    }

    pub fn overlaps(&self, item: &impl Span2D<DimType1 = V1, DimType2 = V2>) -> Vec<&T> {
        self.nodes_overlaps(item)
            .into_iter()
            .flat_map(|(i, node)| node.members.iter())
            .filter(|m| m.overlaps(item))
            .collect()
    }

    fn nodes_overlaps(
        &self,
        item: &impl Span2D<DimType1 = V1, DimType2 = V2>,
    ) -> Vec<(usize, &QuadTreeNode<V1, V2, T>)> {
        if self.is_empty() {
            return Vec::new();
        }
        let root = self.root();
        let mut queue = VecDeque::new();
        let mut nodes = Vec::new();

        if root.overlaps(item) {
            queue.push_back((0usize, root));
        }

        while let Some((i, node)) = queue.pop_front() {
            if node.contains_interval(item) {
                queue.extend(
                    node.children
                        .into_iter()
                        .flatten()
                        .filter_map(|i| self.nodes.get(i).map(|node| (i, node))),
                );
                nodes.push((i, node));
            }
        }

        nodes
    }

    pub fn as_bounding_box(&self) -> BoundingBox<V1, V2> {
        if self.is_empty() {
            BoundingBox::new((V1::zero(), V2::zero()), (V1::zero(), V2::zero()))
        } else {
            self.root().as_bounding_box()
        }
    }
}

#[derive(Debug)]
pub struct QueryIter<
    'a,
    V1: Num + HasProximity,
    V2: Num + HasProximity,
    T: Span2D<DimType1 = V1, DimType2 = V2>,
    Q: Span2D<DimType1 = V1, DimType2 = V2>,
    P: NodeSelectorCriterion<V1, V2, T, Q>,
> {
    ivtree: &'a QuadTree<V1, V2, T>,
    queue: VecDeque<&'a QuadTreeNode<V1, V2, T>>,
    current_node: Option<&'a QuadTreeNode<V1, V2, T>>,
    query: Q,
    predicate: P,
    i: usize,
}

impl<
        V1: Num + HasProximity,
        V2: Num + HasProximity,
        T: Span2D<DimType1 = V1, DimType2 = V2>,
        Q: Span2D<DimType1 = V1, DimType2 = V2>,
        P: NodeSelectorCriterion<V1, V2, T, Q>,
    > FusedIterator for QueryIter<'_, V1, V2, T, Q, P>
{
}

impl<
        'a,
        V1: Num + HasProximity,
        V2: Num + HasProximity,
        T: Span2D<DimType1 = V1, DimType2 = V2>,
        Q: Span2D<DimType1 = V1, DimType2 = V2>,
        P: NodeSelectorCriterion<V1, V2, T, Q>,
    > Iterator for QueryIter<'a, V1, V2, T, Q, P>
{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        self.find_next_value()
    }
}

impl<
        'a,
        V1: Num + HasProximity,
        V2: Num + HasProximity,
        T: Span2D<DimType1 = V1, DimType2 = V2>,
        Q: Span2D<DimType1 = V1, DimType2 = V2>,
        P: NodeSelectorCriterion<V1, V2, T, Q>,
    > QueryIter<'a, V1, V2, T, Q, P>
{
    pub fn new(ivtree: &'a QuadTree<V1, V2, T>, query: Q) -> Self {
        let queue = VecDeque::from(vec![ivtree.root()]);
        let mut this = Self {
            ivtree,
            queue,
            current_node: None,
            query,
            i: 0,
            predicate: P::new(),
        };
        this.find_next_node();
        this
    }

    fn find_value_from_node(&mut self, node: &'a QuadTreeNode<V1, V2, T>) -> Option<&'a T> {
        if self.i < node.members.len() {
            if let Some((i, v)) = node
                .members
                .iter()
                .enumerate()
                .skip(self.i)
                .find(|(_, v)| self.predicate.item_predicate(v, &self.query))
            {
                self.i = i + 1;
                Some(v)
            } else {
                self.i = node.members.len();
                None
            }
        } else {
            None
        }
    }

    fn find_next_value(&mut self) -> Option<&'a T> {
        let mut item = None;
        if let Some(node) = self.current_node {
            if let Some(value) = self.find_value_from_node(node) {
                item = Some(value)
            } else {
                while self.find_next_node() {
                    item = self.find_value_from_node(self.current_node.unwrap());
                    if item.is_some() {
                        break;
                    }
                }
            }
        }
        item
    }

    fn find_next_node(&mut self) -> bool {
        if let Some(node) = self.current_node {
            for c in node.children.into_iter().flatten() {
                let c = &self.ivtree.nodes[c];
                if self.predicate.node_predicate(c, &self.query) {
                    self.queue.push_back(c);
                }
            }
        }
        self.current_node = None;
        self.i = 0;
        self.current_node = self.queue.pop_front();
        self.current_node.is_some()
    }
}

pub trait NodeSelectorCriterion<
    V1: Num + HasProximity,
    V2: Num + HasProximity,
    T: Span2D<DimType1 = V1, DimType2 = V2>,
    Q: Span2D<DimType1 = V1, DimType2 = V2>,
>
{
    fn new() -> Self;

    fn node_predicate(&self, node: &QuadTreeNode<V1, V2, T>, query: &Q) -> bool;

    fn item_predicate(&self, interval: &T, query: &Q) -> bool;
}

#[derive(Debug, Clone, Copy)]
pub struct OverlapPredicate<
    V1: Num + HasProximity,
    V2: Num + HasProximity,
    T: Span2D<DimType1 = V1, DimType2 = V2>,
    Q: Span2D<DimType1 = V1, DimType2 = V2>,
> {
    _v1: PhantomData<V1>,
    _v2: PhantomData<V2>,
    _t: PhantomData<T>,
    _q: PhantomData<Q>,
}

impl<
        V1: Num + HasProximity,
        V2: Num + HasProximity,
        T: Span2D<DimType1 = V1, DimType2 = V2>,
        Q: Span2D<DimType1 = V1, DimType2 = V2>,
    > NodeSelectorCriterion<V1, V2, T, Q> for OverlapPredicate<V1, V2, T, Q>
{
    #[inline(always)]
    fn node_predicate(&self, node: &QuadTreeNode<V1, V2, T>, query: &Q) -> bool {
        node.overlaps(query)
    }

    #[inline(always)]
    fn item_predicate(&self, interval: &T, query: &Q) -> bool {
        interval.overlaps(query)
    }

    fn new() -> Self {
        Self {
            _v1: PhantomData,
            _v2: PhantomData,
            _t: PhantomData,
            _q: PhantomData,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ContainsPredicate<
    V1: Num + PartialOrd + Copy,
    V2: Num + PartialOrd + Copy,
    T: Span2D<DimType1 = V1, DimType2 = V2>,
    Q: Span2D<DimType1 = V1, DimType2 = V2>,
> {
    _v1: PhantomData<V1>,
    _v2: PhantomData<V2>,
    _t: PhantomData<T>,
    _q: PhantomData<Q>,
}

impl<
        V1: Num + HasProximity,
        V2: Num + HasProximity,
        T: Span2D<DimType1 = V1, DimType2 = V2>,
        Q: Span2D<DimType1 = V1, DimType2 = V2>,
    > NodeSelectorCriterion<V1, V2, T, Q> for ContainsPredicate<V1, V2, T, Q>
{
    #[inline(always)]
    fn node_predicate(&self, node: &QuadTreeNode<V1, V2, T>, query: &Q) -> bool {
        let (v1, v2) = query.start();
        node.contains((&v1, &v2))
    }

    #[inline(always)]
    fn item_predicate(&self, interval: &T, query: &Q) -> bool {
        let (v1, v2) = query.start();
        interval.contains((&v1, &v2))
    }

    fn new() -> Self {
        Self {
            _v1: PhantomData,
            _v2: PhantomData,
            _t: PhantomData,
            _q: PhantomData,
        }
    }
}

#[derive(Debug)]
pub struct PreorderIter<
    'a,
    V1: Num + HasProximity,
    V2: Num + HasProximity,
    T: Span2D<DimType1 = V1, DimType2 = V2>,
> {
    ivtree: &'a QuadTree<V1, V2, T>,
    stack: VecDeque<usize>,
}

impl<
        'a,
        V1: Num + HasProximity,
        V2: Num + HasProximity,
        T: Span2D<DimType1 = V1, DimType2 = V2>,
    > Iterator for PreorderIter<'a, V1, V2, T>
{
    type Item = &'a QuadTreeNode<V1, V2, T>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_node()
    }
}

impl<
        'a,
        V1: Num + HasProximity,
        V2: Num + HasProximity,
        T: Span2D<DimType1 = V1, DimType2 = V2>,
    > PreorderIter<'a, V1, V2, T>
{
    pub fn new(ivtree: &'a QuadTree<V1, V2, T>) -> Self {
        let stack = if ivtree.is_empty() {
            VecDeque::new()
        } else {
            VecDeque::from(vec![0])
        };
        Self { ivtree, stack }
    }

    fn next_node(&mut self) -> Option<&'a QuadTreeNode<V1, V2, T>> {
        match self.stack.pop_back() {
            Some(index) => {
                let node = &self.ivtree.nodes[index];
                self.stack.extend(node.children.into_iter().flatten());
                Some(node)
            }
            None => None,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_bbox() {
        let x = BoundingBox::new((5.0, 3.0), (10.0, 6.0));
        assert_eq!(x.start(), (5.0, 3.0));
        assert_eq!(x.end(), (10.0, 6.0));

        assert!(x.overlaps(&x));

        let y = BoundingBox::new((7.0, 3.0), (14.0, 6.0));

        assert_eq!(y.start(), (7.0, 3.0));
        assert_eq!(y.end(), (14.0, 6.0));

        let z = x.combine(&y);
        let q = y.combine(&x);
        assert_eq!(z, q);
        assert!(z.overlaps(&q));

        assert_eq!(z.start(), (5.0, 3.0));
        assert_eq!(z.end(), (14.0, 6.0));

        assert!(z.overlaps(&x));
        let (a, b) = &x.end();
        assert!(z.contains((a, b)));
        let (a, b) = &x.centroid();
        assert!(z.contains((a, b)));
        assert!(z.contains_interval(&x));
    }

    #[test]
    fn test_quad() {
        let x = BoundingBox::new((5.0, 3.0), (10.0, 6.0));
        let y = BoundingBox::new((7.0, 3.0), (14.0, 6.0));
        let z = x.combine(&y);

        let node: QuadTreeNode<f64, f64, BoundingBox<f64, f64>> = QuadTreeNode::from(z);

        assert_eq!(node.as_bounding_box(), z);
    }

    #[test]
    fn test_tree() {
        let mut tree = {
            let x = BoundingBox::new((5.0, 3.0), (10.0, 6.0));
            let y = BoundingBox::new((7.0, 3.0), (14.0, 6.0));
            let z = BoundingBox::new((7.0, 4.0), (12.0, 5.0));
            let q = BoundingBox::new((3.0, 2.0), (8.0, 5.0));
            QuadTree::from_iter([x, y, z, q])
        };
        let bb = tree.as_bounding_box();

        let e = BoundingBox::new((8.0, 3.0), (10.0, 7.0));
        tree.insert(e);
        let bb2 = tree.as_bounding_box();
        assert_ne!(bb, bb2);

        let (a, b) = &(6.0, 4.0);
        let items = tree.contains_point((a, b));
        assert_eq!(items.len(), 2);

        assert_eq!(tree.iter().flatten().count(), 5);

        {
            let x = BoundingBox::new((5.0, 3.0), (10.0, 6.0));
            let y = BoundingBox::new((7.0, 3.0), (14.0, 6.0));
            let z = BoundingBox::new((7.0, 4.0), (12.0, 5.0));
            let q = BoundingBox::new((3.0, 2.0), (8.0, 5.0));
            tree.insert(x);
            tree.insert(y);
            tree.insert(z);
            tree.insert(q);
        }
        assert_eq!(tree.iter().flatten().count(), 9);

        let items = tree.contains_point((a, b));
        assert_eq!(items.len(), 4);
    }

    #[test]
    fn test_partition() {
        let master_bbox = BoundingBox::new((0.0, 0.0), (15.0, 15.0));

        fn split_bbox(bbox: &BoundingBox<f64, f64>) -> [BoundingBox<f64, f64>; 4] {
            let center = bbox.centroid();
            let mut upper_left = BoundingBox::new(
                (bbox.start.0, bbox.start.1),
                (center.0 + 1.0, center.1 + 1.0),
            );
            let mut lower_right = BoundingBox::new(
                ((center.0 - 1.0).max(0.0), (center.1 - 1.0).max(0.0)),
                (bbox.end.0, bbox.end.1),
            );
            let mut lower_left = BoundingBox::new(
                (bbox.start.0, (center.1 - 1.0).max(0.0)),
                (center.0 + 1.0, bbox.end.1),
            );
            let mut upper_right = BoundingBox::new(
                ((center.0 - 1.0).max(0.0), bbox.start.1),
                (bbox.end.0, center.1 + 1.0),
            );
            let quads = [upper_left, upper_right, lower_left, lower_right];
            for q in quads.iter() {
                assert!(q.start.0 >= 0.0);
                assert!(q.start.1 >= 0.0);
            }
            quads
        }

        let boxes: Vec<_> = split_bbox(&master_bbox)
            .map(|b| split_bbox(&b))
            .as_flattened()
            .iter()
            .flat_map(split_bbox)
            .flat_map(|b| split_bbox(&b))
            .flat_map(|b| split_bbox(&b))
            .flat_map(|b| split_bbox(&b))
            .collect();
        let mut tree = QuadTree::empty();
        for bbox in boxes {
            tree.insert(bbox);
        }

        assert_eq!(tree.as_bounding_box(), master_bbox);

        let ext = BoundingBox::new((16.0, 13.0), (18.0, 15.0));
        let before = tree.overlaps(&ext).len();
        tree.insert(ext);
        let after = tree.overlaps(&ext);
        assert_eq!(before + 1, after.len());
        let after = tree.overlaps_iter(&ext);
        assert_eq!(before + 1, after.count());
    }
}
