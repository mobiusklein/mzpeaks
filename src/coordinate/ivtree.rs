use std::{
    borrow::Borrow,
    collections::VecDeque,
    iter::Sum,
    mem
};

use num_traits::real::Real;

use super::{Span1D, SimpleInterval};

#[allow(unused)]
pub fn intervals_containg_point<V, Q: Borrow<V>, T: Span1D<DimType = V>, P: Borrow<T>>(
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

#[allow(unused)]
pub fn intervals_overlapping<
    V,
    T: Span1D<DimType = V>,
    D: Span1D<DimType = V>,
    Q: Borrow<D>,
    P: Borrow<T>,
>(
    intervals: &[P],
    value: Q,
) -> Vec<&T> {
    let mut result = Vec::new();
    for i in intervals.iter() {
        if i.borrow().overlaps(value.borrow()) {
            result.push(i.borrow());
        }
    }
    result
}

/// A node in [`IntervalTree`] over `T`
#[derive(Debug, Clone, Default)]
pub struct IntervalTreeNode<V: Real + Copy + Sum, T: Span1D<DimType = V>> {
    pub start: V,
    pub end: V,
    pub center: V,
    pub level: u32,
    pub members: Vec<T>,
    pub parent: Option<usize>,
    pub left_child: Option<usize>,
    pub right_child: Option<usize>,
}

impl<V: Real + Copy + Sum, T: Span1D<DimType = V>> Span1D for IntervalTreeNode<V, T> {
    type DimType = V;

    fn start(&self) -> Self::DimType {
        self.start
    }

    fn end(&self) -> Self::DimType {
        self.end
    }
}

impl<V: Real + Sum, T: Span1D<DimType = V>> IntervalTreeNode<V, T> {
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
}

enum BuildeTreeSide {
    Left,
    Right,
}

/// An interval tree over `T`
#[derive(Debug, Clone)]
pub struct IntervalTree<V: Real + Copy + Sum, T: Span1D<DimType = V>> {
    pub nodes: Vec<IntervalTreeNode<V, T>>,
}

impl<'members, V: Real + Copy + Sum, T: Span1D<DimType = V>> IntervalTree<V, T> {
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    pub fn root(&self) -> &IntervalTreeNode<V, T> {
        &self.nodes[0]
    }

    pub fn insert(&mut self, interval: T) -> usize {
        if self.is_empty() {
            self.nodes.push(IntervalTreeNode::new(
                (interval.start() + interval.end()) / V::from(2.0).unwrap(),
                vec![interval],
                0,
                None,
                None,
                None,
            ));
            return 0;
        }
        let mut index = 0;
        let insert_in: usize;
        loop {
            let node = &self.nodes[index];
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

    pub fn overlaps(&'members self, start: V, end: V) -> Vec<&'members T> {
        let mut results: Vec<&'members T> = Vec::new();

        if self.nodes.is_empty() {
            return results;
        }

        let q = SimpleInterval::new(start, end);

        if self.nodes[0].overlaps(&q) {
            return results;
        }

        let mut queue = VecDeque::new();
        queue.push_back(&self.nodes[0]);

        while !queue.is_empty() {
            let node = queue.pop_front().unwrap();
            for i in node.members.iter() {
                if i.overlaps(&q) {
                    results.push(i);
                }
            }
            if let Some(left_index) = node.left_child {
                let left = &self.nodes[left_index];
                if left.overlaps(&q) {
                    queue.push_back(left);
                }
            }
            if let Some(right_index) = node.right_child {
                let right = &self.nodes[right_index];
                if right.overlaps(&q) {
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

impl<V: Real + Copy + Sum + Default, T: Span1D<DimType = V>> Default for IntervalTree<V, T> {
    fn default() -> Self {
        let node = IntervalTreeNode::new(V::zero(), vec![], 0, None, None, None);
        Self { nodes: vec![node] }
    }
}

impl<V: Real + Copy + Sum, T: Span1D<DimType = V>> Span1D for IntervalTree<V, T> {
    type DimType = V;

    fn start(&self) -> Self::DimType {
        self.nodes[0].start()
    }

    fn end(&self) -> Self::DimType {
        self.nodes[0].end()
    }
}

#[derive(Debug, Clone)]
pub struct PreorderIter<'a, V: Real + Copy + Sum, T: Span1D<DimType = V>> {
    tree: &'a IntervalTree<V, T>,
    stack: VecDeque<usize>,
}

impl<'a, V: Real + Copy + Sum, T: Span1D<DimType = V>> Iterator for PreorderIter<'a, V, T> {
    type Item = &'a IntervalTreeNode<V, T>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_node()
    }
}

impl<'a, V: Real + Copy + Sum, T: Span1D<DimType = V>> PreorderIter<'a, V, T> {
    pub fn new(tree: &'a IntervalTree<V, T>) -> Self {
        Self {
            tree,
            stack: VecDeque::from(vec![0]),
        }
    }

    pub fn next_node(&mut self) -> Option<&'a IntervalTreeNode<V, T>> {
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

#[cfg(test)]
mod test {
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
        let res = intervals_containg_point(&ivs[..], 2.5f64);
        assert_eq!(res.len(), 2);
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
        let tree = IntervalTree::new(ivs.clone());
        for (i, node) in tree.iter().enumerate() {
            eprintln!("Node {i}: {node:?}");
        }
        let spanning = tree.contains_point(1.0);
        assert_eq!(spanning.len(), 2);

        let spanning = tree.contains_point(7.0);
        assert_eq!(spanning.len(), 4);

        let ivs2: Vec<&SimpleInterval<f64>> = ivs.iter().collect();
        let tree = IntervalTree::new(ivs2);
        let spanning = tree.contains_point(1.0);
        assert_eq!(spanning.len(), 2);

        let spanning = tree.contains_point(7.0);
        assert_eq!(spanning.len(), 4);
    }
}
