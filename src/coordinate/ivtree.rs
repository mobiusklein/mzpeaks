use std::{
    borrow::Borrow,
    collections::VecDeque,
    fmt::Debug,
    iter::{FusedIterator, Sum},
    marker::PhantomData,
    mem,
};

use num_traits::real::Real;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use super::{HasProximity, SimpleInterval, Span1D};

#[allow(unused)]
fn intervals_containing_point<V, Q: Borrow<V>, T: Span1D<DimType = V>, P: Borrow<T>>(
    intervals: &[P],
    value: Q,
) -> Vec<&T> {
    let mut result = Vec::new();
    for i in intervals.iter() {
        if i.borrow().contains(value.borrow()) {
            result.push(i.borrow());
        }
    }
    result
}

/// A node in [`IntervalTree`] over `T`
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct IntervalTreeNode<V: Real + Sum + HasProximity, T: Span1D<DimType = V>> {
    pub start: V,
    pub end: V,
    pub center: V,
    pub level: u32,
    pub members: Vec<T>,
    pub parent: Option<usize>,
    pub left_child: Option<usize>,
    pub right_child: Option<usize>,
}

impl<V: Real + Sum + HasProximity, T: Span1D<DimType = V>> IntoIterator for IntervalTreeNode<V, T> {
    type Item = T;

    type IntoIter = std::vec::IntoIter<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.members.into_iter()
    }
}

impl<'a, V: Real + Sum + HasProximity, T: Span1D<DimType = V>> IntoIterator
    for &'a IntervalTreeNode<V, T>
{
    type Item = &'a T;

    type IntoIter = std::slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.members.iter()
    }
}

impl<V: Real + Sum + HasProximity, T: Span1D<DimType = V>> Span1D for IntervalTreeNode<V, T> {
    type DimType = V;

    fn start(&self) -> Self::DimType {
        self.start
    }

    fn end(&self) -> Self::DimType {
        self.end
    }
}

impl<V: Real + Sum + HasProximity, T: Span1D<DimType = V>> IntervalTreeNode<V, T> {
    pub fn new(
        center: V,
        members: Vec<T>,
        level: u32,
        parent: Option<usize>,
        left_child: Option<usize>,
        right_child: Option<usize>,
    ) -> IntervalTreeNode<V, T> {
        let mut inst = Self {
            center,
            members,
            parent,
            left_child,
            right_child,
            level,
            start: V::max_value(),
            end: -V::max_value(),
        };

        if inst.members.is_empty() {
            inst.start = inst.center;
            inst.end = inst.center;
        } else {
            for interval in inst.members.iter() {
                let i_start = interval.start();
                if i_start < inst.start {
                    inst.start = i_start;
                }
                let i_end = interval.end();
                if i_end > inst.end {
                    inst.end = i_end;
                }
            }
        }
        inst
    }

    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.members.iter()
    }
}

enum BuildeTreeSide {
    Left,
    Right,
}

/// An interval tree over `T`
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct IntervalTree<V: Real + Sum + HasProximity, T: Span1D<DimType = V>> {
    pub(crate) nodes: Vec<IntervalTreeNode<V, T>>,
}

impl<V: Real + Sum + HasProximity, T: Span1D<DimType = V>> Extend<T> for IntervalTree<V, T> {
    fn extend<I: IntoIterator<Item = T>>(&mut self, iter: I) {
        for item in iter {
            self.insert(item);
        }
    }
}

impl<V: Real + Sum + HasProximity, T: Span1D<DimType = V>> FromIterator<T> for IntervalTree<V, T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self::new(Vec::from_iter(iter))
    }
}

impl<'members, V: Real + Sum + HasProximity, T: Span1D<DimType = V>> IntervalTree<V, T> {
    pub fn nodes(&self) -> &[IntervalTreeNode<V, T>] {
        &self.nodes
    }

    pub fn empty() -> Self {
        Self::new(vec![])
    }

    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty() || !self.nodes.iter().any(|n| !n.members.is_empty())
    }

    pub fn root(&self) -> &IntervalTreeNode<V, T> {
        &self.nodes[0]
    }

    pub fn insert(&mut self, interval: T) -> usize {
        let mut index = 0;
        let insert_in: usize;

        loop {
            let node = &self.nodes[index];
            let this_index = index;

            /*
            This node spans the beginning of the interval, meaning that both its left
            and right children might also span the interval, so we may visit both of
            them
            */
            if node.contains(&interval.start()) {
                let (left_index, left_spans) = match node.left_child {
                    Some(left_index) => {
                        let left_node = &self.nodes[left_index];
                        if left_node.contains(&interval.end()) {
                            (left_index, true)
                        } else {
                            (left_index, false)
                        }
                    }
                    None => (index, false),
                };
                let (right_index, right_spans) = match node.right_child {
                    Some(right_index) => {
                        let right_node = &self.nodes[right_index];
                        if right_node.contains(&interval.end()) {
                            (right_index, true)
                        } else {
                            (right_index, false)
                        }
                    }
                    None => (index, false),
                };
                let dest = match (left_spans, right_spans) {
                    (true, false) => left_index,
                    (false, true) => right_index,
                    (false, false) | (true, true) => index,
                };
                if dest == index {
                    insert_in = index;
                    break;
                } else {
                    index = dest;
                }
            }
            /*
            This node contains the end of the interval, which means its right child may be contain it,
            but the left child would already have been covered by the above condition
            */
            if node.contains(&interval.end()) {
                let (right_index, right_spans) = match node.right_child {
                    Some(right_index) => {
                        let right_node = &self.nodes[right_index];
                        if right_node.contains(&interval.start()) {
                            (right_index, true)
                        } else {
                            (right_index, false)
                        }
                    }
                    None => (index, false),
                };
                let dest = if right_spans { right_index } else { index };
                if dest == index {
                    insert_in = index;
                    break;
                } else {
                    index = dest;
                }
            }

            /*
            If neither child contains the interval, then it must belong in the parent node. This means
            the index won't have been updated by either of the earlier blocks, so we can exit the loop
            and insert the interval here.
            */
            if index == this_index {
                insert_in = index;
                break;
            }
        }

        let mut changed = if self.nodes[insert_in].start() > interval.start() {
            self.nodes[insert_in].start = interval.start();
            true
        } else {
            false
        };

        changed |= if self.nodes[insert_in].end() < interval.end() {
            self.nodes[insert_in].end = interval.end();
            true
        } else {
            false
        };
        self.nodes[insert_in].members.push(interval);
        if changed {
            let start = self.nodes[insert_in].start();
            let end = self.nodes[insert_in].end();
            let mut up = self.nodes[insert_in].parent;
            loop {
                match up {
                    Some(parent_index) => {
                        if self.nodes[parent_index].start > start {
                            self.nodes[parent_index].start = start;
                        }
                        if self.nodes[parent_index].end < end {
                            self.nodes[parent_index].end = end;
                        }
                        up = self.nodes[parent_index].parent;
                    }
                    None => {
                        break;
                    }
                }
            }
        }
        insert_in
    }

    pub fn contains_iter(
        &'members self,
        value: V,
    ) -> QueryIter<'members, V, T, SimpleInterval<V>, ContainsPredicate<V, T, SimpleInterval<V>>>
    {
        QueryIter::new(self, SimpleInterval::new(value, value))
    }

    /// A mutable version of [`IntervalTree::contains_iter`].
    ///
    /// See the safety warning on [`QueryIterMut`]
    pub fn contains_iter_mut(
        &'members mut self,
        value: V,
    ) -> QueryIterMut<'members, V, T, SimpleInterval<V>, ContainsPredicate<V, T, SimpleInterval<V>>>
    {
        QueryIterMut::new(self, SimpleInterval::new(value, value))
    }

    pub fn contains_point(&'members self, value: V) -> Vec<&'members T> {
        let mut results: Vec<&'members T> = Vec::new();

        if self.nodes.is_empty() {
            return results;
        }
        if !self.nodes[0].contains(&value) {
            return results;
        }

        let mut queue = VecDeque::new();
        queue.push_back(&self.nodes[0]);

        while !queue.is_empty() {
            let node = queue.pop_front().unwrap();
            for i in node.members.iter() {
                if i.contains(&value) {
                    results.push(i);
                }
            }
            if let Some(left_index) = node.left_child {
                let left = &self.nodes[left_index];
                if left.contains(&value) {
                    queue.push_back(left);
                }
            }
            if let Some(right_index) = node.right_child {
                let right = &self.nodes[right_index];
                if right.contains(&value) {
                    queue.push_back(right);
                }
            }
        }
        results
    }

    pub fn overlaps_range(&'members self, start: V, end: V) -> Vec<&'members T> {
        let q = SimpleInterval::new(start, end);
        self.overlaps(&q)
    }

    pub fn overlaps_iter<Q: Span1D<DimType = V>>(
        &'members self,
        query: Q,
    ) -> QueryIter<'members, V, T, Q, OverlapPredicate<V, T, Q>> {
        QueryIter::new(self, query)
    }

    /// A mutable version of [`IntervalTree::overlaps_iter`].
    ///
    /// See the safety warning on [`QueryIterMut`]
    pub fn overlaps_iter_mut<Q: Span1D<DimType = V>>(
        &'members mut self,
        query: Q,
    ) -> QueryIterMut<'members, V, T, Q, OverlapPredicate<V, T, Q>> {
        QueryIterMut::new(self, query)
    }

    pub fn overlaps<Q: Span1D<DimType = V>>(&'members self, span: Q) -> Vec<&'members T> {
        let mut results: Vec<&'members T> = Vec::new();

        if self.nodes.is_empty() {
            return results;
        }

        if !self.nodes[0].overlaps(&span) {
            return results;
        }

        let mut queue = VecDeque::new();
        queue.push_back(&self.nodes[0]);

        while !queue.is_empty() {
            let node = queue.pop_front().unwrap();
            for i in node.members.iter() {
                if i.overlaps(&span) {
                    results.push(i);
                }
            }
            if let Some(left_index) = node.left_child {
                let left = &self.nodes[left_index];
                if left.overlaps(&span) {
                    queue.push_back(left);
                }
            }
            if let Some(right_index) = node.right_child {
                let right = &self.nodes[right_index];
                if right.overlaps(&span) {
                    queue.push_back(right);
                }
            }
        }
        results
    }

    pub fn flatten(&'members self) -> Vec<&'members T> {
        let mut results = Vec::new();
        if self.is_empty() {
            return results;
        }
        let mut stack = VecDeque::new();
        stack.push_back(self.root());
        while !stack.is_empty() {
            let node = stack.pop_back().unwrap();
            results.extend(node.members.iter());
            if let Some(right_index) = node.right_child {
                stack.push_back(&self.nodes[right_index]);
            }
            if let Some(left_index) = node.left_child {
                stack.push_back(&self.nodes[left_index]);
            }
        }
        results
    }

    pub fn drain(&mut self) -> Vec<T> {
        let mut results = Vec::new();
        if self.is_empty() {
            return results;
        }
        let mut stack = VecDeque::new();
        stack.push_back(0);
        while !stack.is_empty() {
            let index = stack.pop_back().unwrap();
            results.append(&mut self.nodes[index].members);
            if let Some(right_index) = self.nodes[index].right_child {
                stack.push_back(right_index);
            }
            if let Some(left_index) = self.nodes[index].left_child {
                stack.push_back(left_index);
            }
        }
        results
    }

    pub fn balance(&mut self) {
        let members = self.drain();
        let new = Self::new(members);
        self.nodes = new.nodes;
    }

    pub fn iter(&self) -> PreorderIter<V, T> {
        PreorderIter::new(self)
    }

    pub fn new(intervals: Vec<T>) -> IntervalTree<V, T> {
        let root: IntervalTreeNode<V, T> =
            IntervalTreeNode::new(V::zero(), vec![], 0, None, None, None);
        if intervals.is_empty() {
            return IntervalTree { nodes: vec![root] };
        }
        let mut nodes: Vec<IntervalTreeNode<V, T>> = Vec::new();
        nodes.push(root);

        let mut stack: VecDeque<(usize, Vec<T>, BuildeTreeSide)> = VecDeque::new();
        let entry = (0, intervals, BuildeTreeSide::Left);
        stack.push_back(entry);

        while !stack.is_empty() {
            if let Some((parent, members, side)) = stack.pop_back() {
                let n = members.len();
                let center = if n > 0 {
                    let acc: V = members
                        .iter()
                        .map(|i| (i.start() + i.end()) / V::from(2.0).unwrap())
                        .sum();
                    acc / (V::from(n + 1).unwrap())
                } else {
                    V::zero()
                };
                let mut left: Vec<T> = Vec::new();
                let mut right: Vec<T> = Vec::new();
                let mut contained: Vec<T> = Vec::new();

                if n < 5 {
                    contained = members;
                } else {
                    let diff = V::from(1e-6).unwrap();
                    for rec in members {
                        // The interval is essentially a point around the center
                        if (rec.start() - center).abs() < diff && (rec.end() - center).abs() < diff
                        {
                            contained.push(rec)
                        // The interval is to the left of the center, closing before it.
                        } else if center > rec.end() {
                            left.push(rec)
                        // The interval is to the right of center, starting after it.
                        } else if center < rec.start() {
                            right.push(rec)
                        // The interval spans the center, contained in this node
                        } else {
                            contained.push(rec)
                        }
                    }
                }

                // Sometimes a region is very narrow but does't stack up exactly on the center so only
                // one side is populated. Force this to become the node-contained data.
                if contained.is_empty() {
                    if left.is_empty() && !right.is_empty() {
                        mem::swap(&mut contained, &mut right);
                    } else if !left.is_empty() && right.is_empty() {
                        mem::swap(&mut contained, &mut left);
                    }
                }

                let level = nodes[parent].level + 1;
                let node_index = nodes.len();
                let node =
                    IntervalTreeNode::new(center, contained, level, Some(parent), None, None);
                match side {
                    BuildeTreeSide::Left => {
                        let start = node.start;
                        nodes[parent].left_child = Some(node_index);
                        let mut up = parent;
                        loop {
                            let p = &mut nodes[up];
                            p.start = V::min(p.start, start);
                            if let Some(next) = p.parent {
                                if up == next {
                                    break;
                                }
                                up = next
                            } else {
                                break;
                            }
                        }
                    }
                    BuildeTreeSide::Right => {
                        let end = node.end;
                        nodes[parent].right_child = Some(node_index);
                        let mut up = parent;
                        loop {
                            let p = &mut nodes[up];
                            p.end = V::max(p.end, end);
                            if let Some(next) = p.parent {
                                if up == next {
                                    break;
                                }
                                up = next;
                            } else {
                                break;
                            }
                        }
                    }
                }
                nodes.push(node);

                if !left.is_empty() {
                    stack.push_back((node_index, left, BuildeTreeSide::Left))
                }
                if !right.is_empty() {
                    stack.push_back((node_index, right, BuildeTreeSide::Right))
                }
            }
        }

        IntervalTree { nodes }
    }
}

impl<V: Real + Sum + Default + HasProximity, T: Span1D<DimType = V>> Default
    for IntervalTree<V, T>
{
    fn default() -> Self {
        let node = IntervalTreeNode::new(V::zero(), vec![], 0, None, None, None);
        Self { nodes: vec![node] }
    }
}

impl<V: Real + Sum + HasProximity, T: Span1D<DimType = V>> Span1D for IntervalTree<V, T> {
    type DimType = V;

    fn start(&self) -> Self::DimType {
        self.nodes
            .first()
            .map(|x| x.start())
            .unwrap_or_else(|| V::zero())
    }

    fn end(&self) -> Self::DimType {
        self.nodes
            .first()
            .map(|x| x.end())
            .unwrap_or_else(|| V::zero())
    }
}

#[derive(Debug, Clone)]
pub struct PreorderIter<'a, V: Real + Sum + HasProximity, T: Span1D<DimType = V>> {
    tree: &'a IntervalTree<V, T>,
    stack: VecDeque<usize>,
}

impl<'a, V: Real + Sum + HasProximity, T: Span1D<DimType = V>> Iterator for PreorderIter<'a, V, T> {
    type Item = &'a IntervalTreeNode<V, T>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_node()
    }
}

impl<'a, V: Real + Sum + HasProximity, T: Span1D<DimType = V>> PreorderIter<'a, V, T> {
    pub fn new(tree: &'a IntervalTree<V, T>) -> Self {
        let stack = if tree.is_empty() {
            VecDeque::new()
        } else {
            VecDeque::from(vec![0])
        };
        Self { tree, stack }
    }

    pub fn next_index(&self) -> Option<usize> {
        self.stack.back().copied()
    }

    fn next_node(&mut self) -> Option<&'a IntervalTreeNode<V, T>> {
        match self.stack.pop_back() {
            Some(index) => {
                let node = &self.tree.nodes[index];
                if let Some(right_index) = node.right_child {
                    self.stack.push_back(right_index);
                };
                if let Some(left_index) = node.left_child {
                    self.stack.push_back(left_index);
                };
                Some(node)
            }
            None => None,
        }
    }
}

#[derive(Debug)]
pub struct QueryIter<
    'a,
    V: Real + Sum + HasProximity,
    T: Span1D<DimType = V>,
    Q: Span1D<DimType = V> + 'a,
    P: NodeSelectorCriterion<V, T, Q>,
> {
    ivtree: &'a IntervalTree<V, T>,
    queue: VecDeque<&'a IntervalTreeNode<V, T>>,
    current_node: Option<&'a IntervalTreeNode<V, T>>,
    query: Q,
    predicate: P,
    i: usize,
}

impl<
        'a,
        V: Real + Sum + HasProximity,
        T: Span1D<DimType = V>,
        Q: Span1D<DimType = V> + 'a,
        P: NodeSelectorCriterion<V, T, Q>,
    > FusedIterator for QueryIter<'a, V, T, Q, P>
{
}

impl<
        'a,
        V: Real + Sum + HasProximity,
        T: Span1D<DimType = V>,
        Q: Span1D<DimType = V> + 'a,
        P: NodeSelectorCriterion<V, T, Q>,
    > Iterator for QueryIter<'a, V, T, Q, P>
{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        self.find_next_value()
    }
}

impl<
        'a,
        V: Real + Sum + HasProximity,
        T: Span1D<DimType = V>,
        Q: Span1D<DimType = V> + 'a,
        P: NodeSelectorCriterion<V, T, Q>,
    > QueryIter<'a, V, T, Q, P>
{
    pub fn new(ivtree: &'a IntervalTree<V, T>, query: Q) -> Self {
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

    fn find_value_from_node(&mut self, node: &'a IntervalTreeNode<V, T>) -> Option<&'a T> {
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
            if let Some(c) = node.left_child {
                let c = &self.ivtree.nodes[c];
                if self.predicate.node_predicate(c, &self.query) {
                    self.queue.push_back(c);
                }
            }

            if let Some(c) = node.right_child {
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

/// Provides a mutable iterator over a subset of elements in an [`IntervalTree`].
///
/// # Safety
/// This iterator should not be used in a way that alters the interval bounds of the
/// element or else it risks making it unreachable without explicitly calling [`IntervalTree::balance`]
#[derive(Debug)]
pub struct QueryIterMut<
    'a,
    V: Real + Sum + HasProximity,
    T: Span1D<DimType = V>,
    Q: Span1D<DimType = V> + 'a,
    P: NodeSelectorCriterion<V, T, Q>,
> {
    ivtree: &'a mut IntervalTree<V, T>,
    queue: VecDeque<usize>,
    current_node: Option<usize>,
    query: Q,
    predicate: P,
    i: usize,
}

impl<
        'a,
        V: Real + Sum + HasProximity,
        T: Span1D<DimType = V>,
        Q: Span1D<DimType = V> + 'a,
        P: NodeSelectorCriterion<V, T, Q>,
    > QueryIterMut<'a, V, T, Q, P>
{
    pub fn new(ivtree: &'a mut IntervalTree<V, T>, query: Q) -> Self {
        let queue = VecDeque::from(vec![0]);
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

    fn find_value_from_node(&mut self, node: usize) -> Option<usize> {
        let node = &mut self.ivtree.nodes[node];
        let n = node.members.len();
        if self.i < n {
            if let Some((i, _v)) = node
                .members
                .iter()
                .enumerate()
                .skip(self.i)
                .find(|(_, v)| self.predicate.item_predicate(v, &self.query))
            {
                self.i = i + 1;
                Some(i)
            } else {
                self.i = n;
                None
            }
        } else {
            None
        }
    }

    fn find_next_value(&mut self) -> Option<(usize, usize)> {
        if let Some(node) = self.current_node {
            if let Some(value) = self.find_value_from_node(node) {
                return Some((node, value));
            }
            loop {
                if !self.find_next_node() {
                    break;
                }
                let node_i = self.current_node.unwrap();
                let item_i = self.find_value_from_node(self.current_node.unwrap());
                if item_i.is_some() {
                    return Some((node_i, item_i.unwrap()));
                }
            }
        }
        None
    }

    fn find_next_node(&mut self) -> bool {
        if let Some(node) = self.current_node {
            let node = &self.ivtree.nodes[node];
            if let Some(ci) = node.left_child {
                let c = &self.ivtree.nodes[ci];
                if self.predicate.node_predicate(c, &self.query) {
                    self.queue.push_back(ci);
                }
            }

            if let Some(ci) = node.right_child {
                let c = &self.ivtree.nodes[ci];
                if self.predicate.node_predicate(c, &self.query) {
                    self.queue.push_back(ci);
                }
            }
        }
        self.current_node = None;
        self.i = 0;
        self.current_node = self.queue.pop_front();
        self.current_node.is_some()
    }
}

/// # Safety
/// This iterator should not be used in a way that alters the interval bounds of the
/// element or else it risks making it unreachable without explicitly calling [`IntervalTree::balance`]
impl<
        'a,
        V: Real + Sum + HasProximity,
        T: Span1D<DimType = V>,
        Q: Span1D<DimType = V> + 'a,
        P: NodeSelectorCriterion<V, T, Q>,
    > Iterator for QueryIterMut<'a, V, T, Q, P>
{
    type Item = &'a mut T;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((i, j)) = self.find_next_value() {
            if let Some(x) = self
                .ivtree
                .nodes
                .get_mut(i)
                .and_then(|f| f.members.get_mut(j))
            {
                let y = unsafe { &mut *core::ptr::addr_of_mut!(*x) };
                Some(y)
            } else {
                None
            }
        } else {
            None
        }
    }
}

pub trait NodeSelectorCriterion<
    V: Real + Sum + HasProximity,
    T: Span1D<DimType = V>,
    Q: Span1D<DimType = V>,
>
{
    fn new() -> Self;

    fn node_predicate(&self, node: &IntervalTreeNode<V, T>, query: &Q) -> bool;

    fn item_predicate(&self, interval: &T, query: &Q) -> bool;
}

#[derive(Debug, Clone, Copy)]
pub struct OverlapPredicate<
    V: Real + Sum + HasProximity,
    T: Span1D<DimType = V>,
    Q: Span1D<DimType = V>,
> {
    _v: PhantomData<V>,
    _t: PhantomData<T>,
    _q: PhantomData<Q>,
}

impl<V: Real + Sum + HasProximity, T: Span1D<DimType = V>, Q: Span1D<DimType = V>>
    NodeSelectorCriterion<V, T, Q> for OverlapPredicate<V, T, Q>
{
    #[inline(always)]
    fn node_predicate(&self, node: &IntervalTreeNode<V, T>, query: &Q) -> bool {
        node.overlaps(query)
    }

    #[inline(always)]
    fn item_predicate(&self, interval: &T, query: &Q) -> bool {
        interval.overlaps(query)
    }

    fn new() -> Self {
        Self {
            _v: PhantomData,
            _t: PhantomData,
            _q: PhantomData,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ContainsPredicate<V: Real + Sum, T: Span1D<DimType = V>, Q: Span1D<DimType = V>> {
    _v: PhantomData<V>,
    _t: PhantomData<T>,
    _q: PhantomData<Q>,
}

impl<V: Real + Sum + HasProximity, T: Span1D<DimType = V>, Q: Span1D<DimType = V>>
    NodeSelectorCriterion<V, T, Q> for ContainsPredicate<V, T, Q>
{
    #[inline(always)]
    fn node_predicate(&self, node: &IntervalTreeNode<V, T>, query: &Q) -> bool {
        node.contains(&query.start())
    }

    #[inline(always)]
    fn item_predicate(&self, interval: &T, query: &Q) -> bool {
        interval.contains(&query.start())
    }

    fn new() -> Self {
        Self {
            _v: PhantomData,
            _t: PhantomData,
            _q: PhantomData,
        }
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use super::*;

    #[test]
    fn test_contains() {
        let iv = SimpleInterval {
            start: 2.0,
            end: 7.0,
        };
        assert!(!iv.contains(&0.5));
        assert!(iv.contains(&5.0));
    }

    #[test]
    fn test_intervals_contain() {
        let ivs = [
            SimpleInterval::new(0.0, 3.0),
            SimpleInterval::new(2.0, 5.0),
            SimpleInterval::new(5.0, 10.0),
        ];
        let res = intervals_containing_point(&ivs[..], 2.5f64);
        assert_eq!(res.len(), 2);
    }

    #[test]
    fn test_interval_tree_incr() {
        let ivs = vec![
            SimpleInterval::new(0.0, 3.0),
            SimpleInterval::new(2.0, 5.0),
            SimpleInterval::new(5.0, 10.0),
            SimpleInterval::new(0.5, 3.0),
            SimpleInterval::new(3.0, 5.0),
            SimpleInterval::new(5.0, 12.0),
            SimpleInterval::new(5.0, 6.0),
            SimpleInterval::new(7.0, 10.0),
            SimpleInterval::new(7.0, 12.0),
        ];

        let mut tree = IntervalTree::empty();

        // Pressure test the insert API
        for iv in ivs.iter().copied() {
            tree.insert(iv);
        }

        assert_eq!(tree.iter().flatten().count(), ivs.len());

        let spanning = tree.contains_point(1.0);
        assert_eq!(spanning.len(), 2);

        let spanning = tree.contains_point(7.0);
        assert_eq!(spanning.len(), 4);

        // Stress the trees more
        let mut tree = IntervalTree::from_iter(ivs.clone());
        for (i, mut iv) in ivs.iter().copied().enumerate() {
            if i % 2 == 0 {
                iv.start -= 2.0;
                iv.end -= 2.0;
            } else {
                iv.start += 2.0;
                iv.end += 2.0;
            }
            tree.insert(iv);
        }

        for (i, mut iv) in ivs.iter().copied().enumerate() {
            if i % 2 == 0 {
                iv.start -= 6.0;
                iv.end -= 6.0;
            } else {
                iv.start += 6.0;
                iv.end += 6.0;
            }
            tree.insert(iv);
        }

        let mut seen = HashSet::new();
        let mut iter = tree.iter();
        while let Some(i) = iter.next_index() {
            let node = iter.next().unwrap();
            eprintln!(
                "{i}|{}|{}\t{} -- {} -- {}",
                node.parent.map(|s| s.to_string()).unwrap_or_default(),
                node.members.len(),
                node.start,
                node.center,
                node.end
            );
            if let Some(parent_id) = node.parent {
                assert!(seen.contains(&parent_id))
            }
            seen.insert(i);
        }

        eprintln!("{}", "-".repeat(20));

        // Rebalancing the tree should fix the depth
        tree.balance();
        let mut seen = HashSet::new();
        let mut iter = tree.iter();
        while let Some(i) = iter.next_index() {
            let node = iter.next().unwrap();
            eprintln!(
                "{i}|{}|{}\t{} -- {} -- {}",
                node.parent.map(|s| s.to_string()).unwrap_or_default(),
                node.members.len(),
                node.start,
                node.center,
                node.end
            );
            if let Some(parent_id) = node.parent {
                assert!(seen.contains(&parent_id))
            }
            seen.insert(i);
        }

        assert_eq!(tree.start(), -6.0);
        assert_eq!(tree.end(), 18.0);
        let n = tree
            .overlaps_iter_mut(SimpleInterval::new(tree.start(), tree.end()))
            .count();
        assert_eq!(tree.iter().flatten().count(), n);
    }

    #[test]
    fn test_interval_tree() {
        let ivs = vec![
            SimpleInterval::new(0.0, 3.0),
            SimpleInterval::new(2.0, 5.0),
            SimpleInterval::new(5.0, 10.0),
            SimpleInterval::new(0.5, 3.0),
            SimpleInterval::new(3.0, 5.0),
            SimpleInterval::new(5.0, 12.0),
            SimpleInterval::new(5.0, 6.0),
            SimpleInterval::new(7.0, 10.0),
            SimpleInterval::new(7.0, 12.0),
        ];
        let tree = IntervalTree::from_iter(ivs.clone());
        for (i, node) in tree.nodes().iter().enumerate() {
            eprintln!("Node {i}: {node:?}");
        }

        assert_eq!(
            tree.flatten().len(),
            ivs.len(),
            "Flattening the tree didn't produce as many items as original collection"
        );

        // Expected to be spanned by two intervals
        let spanning = tree.contains_point(1.0);
        assert_eq!(spanning.len(), 2);
        assert_eq!(tree.contains_iter(1.0).count(), 2);

        // Expected to be spanned by four intervals
        let spanning = tree.contains_point(7.0);
        assert_eq!(spanning.len(), 4);

        // Repeat the earlier tests with an interval-by-reference tree
        let ivs2: Vec<&SimpleInterval<f64>> = ivs.iter().collect();
        let tree = IntervalTree::new(ivs2);
        let spanning = tree.contains_point(1.0);
        assert_eq!(spanning.len(), 2);

        let spanning = tree.contains_point(7.0);
        assert_eq!(spanning.len(), 4);

        let query_iv = SimpleInterval::new(2.0, 5.0);
        let items = tree.overlaps(&query_iv);
        let items_iterd: Vec<_> = tree.overlaps_iter(&query_iv).collect();
        assert_eq!(items, items_iterd);
    }
}
