// Copyright (c) 2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//               : Susan Hert <hert@mpi-sb.mpg.de>

#ifndef CGAL_RANDOM_POLYGON_2_SWEEP_H
#define CGAL_RANDOM_POLYGON_2_SWEEP_H

#include <CGAL/enum.h>
#include <CGAL/Polygon_2/polygon_assertions.h>
#include <set>
#include <vector>
#include <algorithm>
#include <CGAL/Polygon_2/Polygon_2_simplicity.h>

/*
  A polygon is called simple of no edges intersect each other, except
  consecutive edges, which intersect in their common vertex.
  The test for simplicity is implemented by means of a sweep line algorithm.
  The vertical line is swept from left to right. The edges of the polygon that
  are crossed by the sweep line are stored in a tree from bottom to top.

  We discern three types of events:
  - insertion events. When both edges of a polygon vertex extend to the right
    we need to insert both edges in the tree. We need to search with the vertex
    to find out between which edges the new edges are to be inserted.
  - deletion events. When both edges extend to the left of the vertex we need
    to remove both edges from the tree. We have to check that the vertex lies
    between the edges above and below the removed edges.
  - replacement event. In the other case we need to replace the edge that
    extends to the left by the edge that extends to the right. We need to check
    that the vertex lies between the edges above and below the current edge.

  We represent the tree by a std::set. This is not a perfect fit for the
  operations described above. In particular, the fact that we search with a
  VERTEX, while the set contains EDGES, is not directly supported. The
  insertion of edges is also complicated by the fact that we need to insert
  two edges starting at the same vertex. The order in which they are inserted
  in the tree does matter, because the edges go in separate directions.
  Because of this the set needs a special associated comparison function.
  Every edge has a flag 'is_in_tree', which is true for the edges in the tree
  but false for the edges when they are inserted. The comparison function
  treats the latter type of edge as a vertex. Another flag, is_left_to_right,
  tells which of the two vertices to take. The problem of the coinciding
  vertices is solved by the convention that the highest edges is always
  inserted first. The comparison function uses this fact.

  Vertex indices of the polygon play a double role. The number v can be used to
  identify vertex v or the edge from vertex v to vertex v+1.

*/

namespace CGAL {

namespace i_generator_polygon {
 // namespace CGAL::i_generator_polygon is used for internal functions

typedef i_polygon::Index_t        Index_t;
typedef i_polygon::Vertex_index   Vertex_index;
typedef i_polygon::Vertex_order   Vertex_order;

template <class ForwardIterator, class PolygonTraits>
class Vertex_data ;

template <class ForwardIterator, class PolygonTraits>
class Less_segments {
    Vertex_data<ForwardIterator, PolygonTraits> *m_vertex_data;
    bool less_than_in_tree(Vertex_index i, Vertex_index j);
  public:
    Less_segments(Vertex_data<ForwardIterator, PolygonTraits> *vertex_data)
    : m_vertex_data(vertex_data) {}
   bool operator()(Vertex_index i, Vertex_index j);
};

template <class ForwardIterator, class PolygonTraits>
class Vertex_data : 
       public i_polygon::Vertex_data_base<ForwardIterator, PolygonTraits>
{
public:
    typedef Less_segments<ForwardIterator, PolygonTraits>    Less_segs;
    typedef std::set<Vertex_index, Less_segs>                Tree;
    typedef i_polygon::Edge_data<Less_segs>                  Edge_data;
    typedef i_polygon::Vertex_data_base<ForwardIterator, PolygonTraits>
                                                             Base_class;

    using Base_class::next;
    using Base_class::prev;
    using Base_class::point;
    using Base_class::index_at_rank;
    using Base_class::ordered_left_to_right;
    using Base_class::orientation_2;

    std::vector<Edge_data> edges;
    Vertex_index conflict1, conflict2; // Intersecting edges.

    Vertex_data(ForwardIterator begin, ForwardIterator end,
                const PolygonTraits& pgnt);
    void init(Tree *tree);
    void left_and_right_index(Vertex_index &left, Vertex_index &right,
            Vertex_index edge);
    Vertex_index left_index(Vertex_index edge)
        { return edges[edge.as_int()].is_left_to_right ? edge : next(edge); }

    void sweep(Tree *tree);
    bool insertion_event(Tree *tree,
                Vertex_index i, Vertex_index j, Vertex_index k);
    bool replacement_event(Tree *tree,
                Vertex_index cur, Vertex_index to_insert);
    bool deletion_event(Tree *tree, Vertex_index i, Vertex_index j);
    bool on_right_side(Vertex_index vt, Vertex_index edge, bool above);
    void find_conflict(Tree *tree, Vertex_index cur_vt,
        typename Tree::iterator seg1, typename Tree::iterator seg2);
    void find_conflict_between(Tree *, Vertex_index cur_vt,
        typename Tree::iterator seg1, typename Tree::iterator seg2);
};

} // end of namespace i_generator_polygon

// ----- implementation of i_generator_polygon functions. -----

namespace i_generator_polygon {
template <class ForwardIterator, class PolygonTraits>
bool Less_segments<ForwardIterator, PolygonTraits>::
operator()(Vertex_index i, Vertex_index j)
{
    if (m_vertex_data->edges[j.as_int()].is_in_tree) {
        return less_than_in_tree(i,j);
    } else {
        return !less_than_in_tree(j,i);
    }
}

template <class ForwardIterator, class PolygonTraits>
bool Less_segments<ForwardIterator, PolygonTraits>::
less_than_in_tree(Vertex_index new_edge, Vertex_index tree_edge)
{
#if defined(CGAL_POLY_GENERATOR_DEBUG)
    std::cout << "less_than_in_tree" << std::endl;
#endif
    CGAL_polygon_precondition(
       m_vertex_data->edges[tree_edge.as_int()].is_in_tree);
    CGAL_polygon_precondition(
       !m_vertex_data->edges[new_edge.as_int()].is_in_tree);
    Vertex_index left, mid, right;
    m_vertex_data->left_and_right_index(left, right, tree_edge);
    mid = m_vertex_data->left_index(new_edge);
    if (mid.as_int() == left.as_int()) {
        return true;
    }
    switch (m_vertex_data->orientation_2( m_vertex_data->point(left),
            m_vertex_data->point(mid), m_vertex_data->point(right))) {
      case LEFT_TURN: return true;
      case RIGHT_TURN: return false;
      case COLLINEAR: break;
    }
    CGAL_polygon_assertion(m_vertex_data->less_xy_2(
                               m_vertex_data->point(left),
                               m_vertex_data->point(mid)));
    CGAL_polygon_assertion( m_vertex_data->less_xy_2(
                                m_vertex_data->point(mid),
                                m_vertex_data->point(right)));
    m_vertex_data->is_simple_result = false;
    Vertex_index mid_succ = m_vertex_data->next(mid);
    if (mid_succ.as_int() <= min(left.as_int(), right.as_int()))
    {
       //
       // have one of these two situations:
       //     x+k  x   x+k+1     OR  x+k+1  x   x+k
       //     left mid right         left   mid right
       // and need to swap x and x+1 (so reverse the range (x-1, x+k])
       //
       m_vertex_data->conflict1 = m_vertex_data->prev(mid);
    }
    else
    {
       // have one of these two situations:
       //     x    x+k x+1      OR   x+1   x+k x
       //     left mid right         left  mid right
       // and need to swap x+1 and x+2 (so reverse range (x, x+k])
       m_vertex_data->conflict1 = mid;
    }
    Vertex_index left_succ = m_vertex_data->next(left);
    if (left_succ.as_int() == right.as_int())
    {
#if defined(CGAL_POLY_GENERATOR_DEBUG)
       std::cout << "conflict2 is left" << std::endl;
#endif
       m_vertex_data->conflict2 = left;
    }
    else
    {
#if defined(CGAL_POLY_GENERATOR_DEBUG)
       std::cout << "conflict2 is left" << std::endl;
#endif
       m_vertex_data->conflict2 = right;
    }
    return true;
}

template <class ForwardIterator, class PolygonTraits>
Vertex_data<ForwardIterator, PolygonTraits>::
Vertex_data(ForwardIterator begin, ForwardIterator end,
            const PolygonTraits& pgn_traits)
  : Base_class(begin, end, pgn_traits) {}

template <class ForwardIterator, class PolygonTraits>
void Vertex_data<ForwardIterator, PolygonTraits>::init(Tree *tree)
{
    // The initialization cannot be done in the constructor,
    // otherwise we copy singular valued iterators.
    edges.insert(edges.end(), this->m_size, Edge_data(tree->end()));
}

template <class ForwardIterator, class PolygonTraits>
void Vertex_data<ForwardIterator, PolygonTraits>::
left_and_right_index(Vertex_index &left, Vertex_index &right,
            Vertex_index edge)
{
    if (edges[edge.as_int()].is_left_to_right) {
        left = edge; right = next(edge);
    } else {
        right = edge; left = next(edge);
    }
}

template <class ForwardIterator, class PolygonTraits>
bool Vertex_data<ForwardIterator, PolygonTraits>::
insertion_event(Tree *tree, Vertex_index prev_vt,
                Vertex_index mid_vt, Vertex_index next_vt)
{
#if defined(CGAL_POLY_GENERATOR_DEBUG)
    std::cout << "insertion_event" << std::endl;
#endif
    // check which endpoint is above the other
    bool left_turn;
    switch(orientation_2(point(prev_vt), point(mid_vt), point(next_vt))) {
      case LEFT_TURN: left_turn = true; break;
      case RIGHT_TURN: left_turn = false; break;
      default: //found conflict prev_vt-seg - mid_vt-seg
#if defined(CGAL_POLY_GENERATOR_DEBUG)
            std::cout << "conflict2 is next_vt" << std::endl;
#endif
            conflict1 = prev_vt;
	    conflict2 = next_vt;
            return false;
      
    }
    Edge_data
        &td_prev = edges[prev_vt.as_int()],
        &td_mid = edges[mid_vt.as_int()];
    td_prev.is_in_tree = false;
    td_prev.is_left_to_right = false;
    td_mid.is_in_tree = false;
    td_mid.is_left_to_right = true;
    // insert the highest chain first
    std::pair<typename Tree::iterator, bool> result;
    if (left_turn) {
        result = tree->insert(prev_vt);
	// CGAL_polygon_assertion(result.second)
	td_prev.tree_it = result.first;
        td_prev.is_in_tree = true;
        if (!this->is_simple_result) return false;
        result = tree->insert(mid_vt);
	// CGAL_polygon_assertion(result.second)
	td_mid.tree_it = result.first;
        td_mid.is_in_tree = true;
        if (!this->is_simple_result) return false;
    } else {
        result = tree->insert(mid_vt);
	// CGAL_polygon_assertion(result.second)
	td_mid.tree_it = result.first;
        td_mid.is_in_tree = true;
        if (!this->is_simple_result) return false;
        result = tree->insert(prev_vt);
	// CGAL_polygon_assertion(result.second)
	td_prev.tree_it = result.first;
        td_prev.is_in_tree = true;
        if (!this->is_simple_result) return false;
    }
    return true;
}

template <class ForwardIterator, class PolygonTraits>
bool Vertex_data<ForwardIterator, PolygonTraits>::
on_right_side(Vertex_index vt, Vertex_index edge_id, bool above)
{
    Orientation turn =
        orientation_2(point(edge_id), point(vt), point(next(edge_id)));
    bool left_turn = edges[edge_id.as_int()].is_left_to_right ? above : !above;
    if (left_turn) {
        if (turn != RIGHT_TURN) {
            return false;
        }
    } else {
        if (turn != LEFT_TURN) {
            return false;
        }
    }
    return true;
}

template <class ForwardIterator, class PolygonTraits>
bool Vertex_data<ForwardIterator, PolygonTraits>::
replacement_event(Tree *tree, Vertex_index cur_edge, Vertex_index next_edge)
{
#if defined(CGAL_POLY_GENERATOR_DEBUG)
    std::cout << "replacement_event" << std::endl;
#endif
    // check if continuation point is on the right side of neighbor segments
    typedef typename Tree::iterator It;
    Edge_data &td = edges[cur_edge.as_int()];
    CGAL_polygon_assertion(td.is_in_tree);
    It cur_seg = td.tree_it;
    Vertex_index cur_vt = (td.is_left_to_right) ? next_edge : cur_edge;
    if (cur_seg != tree->begin()) {
        It seg_below = cur_seg;
	--seg_below;
	if (!on_right_side(cur_vt, *seg_below, true)) {
	    // found conflict cur_seg - seg_below
#if defined(CGAL_POLY_GENERATOR_DEBUG)
            std::cout << "conflict2 is seg_below" << std::endl;
#endif
            conflict1 = *cur_seg;
	    conflict2 = *seg_below;
	    return false;
        }
    }
    It seg_above = cur_seg;
    ++ seg_above;
    if (seg_above != tree->end()) {
        if (!on_right_side(cur_vt, *seg_above, false)) {
	    // found conflict cur_seg - seg_above
#if defined(CGAL_POLY_GENERATOR_DEBUG)
            std::cout << "conflict2 is seg_above" << std::endl;
#endif
            conflict1 = *cur_seg;
	    conflict2 = *seg_above;
	    return false;
        }
    }
    // replace the segment
    Edge_data &new_td =
            edges[next_edge.as_int()];
    new_td.is_left_to_right = td.is_left_to_right;
    new_td.is_in_tree = false;
    tree->erase(cur_seg);
    td.is_in_tree = false;
    new_td.tree_it = tree->insert(seg_above, next_edge);
    new_td.is_in_tree = true;
    return this->is_simple_result;
}

template <class ForwardIterator, class PolygonTraits>
void
Vertex_data<ForwardIterator, PolygonTraits>::
find_conflict(Tree *tree, Vertex_index cur_vt,
              typename Tree::iterator seg1, typename Tree::iterator seg2)
{
    typedef typename Tree::iterator It;
    It cur = seg1;
    ++cur;
    while (cur !=  tree->end() && cur != seg2)
        ++cur;
    if (cur == seg2)
        find_conflict_between(tree, cur_vt, seg1, seg2);
    else
        find_conflict_between(tree, cur_vt, seg2, seg1);
}

template <class ForwardIterator, class PolygonTraits>
void
Vertex_data<ForwardIterator, PolygonTraits>::
find_conflict_between(Tree *, Vertex_index cur_vt,
                      typename Tree::iterator seg1, 
                      typename Tree::iterator seg2)
{
#if defined(CGAL_POLY_GENERATOR_DEBUG)
    std::cout << "find_conflict_between" << std::endl;
#endif
    typedef typename Tree::iterator It;
    It between_seg = seg1;
    ++between_seg;
    if (!on_right_side(cur_vt, *between_seg, false)) {
	// found conflict between_seg - seg1
#if defined(CGAL_POLY_GENERATOR_DEBUG)
        std::cout << "conflict1 is seg1" << std::endl;
#endif
        conflict1 = *seg1;
	conflict2 = *between_seg;
    } else {
	// found conflict between_seg - seg2
#if defined(CGAL_POLY_GENERATOR_DEBUG)
        std::cout << "conflict1 is seg2" << std::endl;
#endif
        conflict1 = *seg2;
	conflict2 = *between_seg;
    }
}

template <class ForwardIterator, class PolygonTraits>
bool Vertex_data<ForwardIterator, PolygonTraits>::
deletion_event(Tree *tree, Vertex_index prev_vt, Vertex_index mid_vt)
{
#if defined(CGAL_POLY_GENERATOR_DEBUG)
    std::cout << "deletion_event" << std::endl;
#endif
    // check if continuation point is on the right side of neighbor segments
    typedef typename Tree::iterator It;
    Edge_data
        &td_prev = edges[prev_vt.as_int()],
        &td_mid = edges[mid_vt.as_int()];
    It prev_seg = td_prev.tree_it, mid_seg = td_mid.tree_it;
    Vertex_index cur_vt = (td_prev.is_left_to_right) ? mid_vt : prev_vt;
    It seg_above = prev_seg;
    ++seg_above;

    if (seg_above == mid_seg) {
        ++seg_above;
    } else {
        // mid_seg was not above prev_seg, so prev_seg should be above mid_seg
        // We check this to see if the edges are really neighbors in the tree.
        It prev_seg_copy = mid_seg;
        ++prev_seg_copy;
        if (prev_seg_copy != prev_seg) {
	    find_conflict(tree, cur_vt, prev_seg, mid_seg);
            return false;
        }
    }
    // remove the segments
    tree->erase(prev_seg);
    td_prev.is_in_tree = false;
    tree->erase(mid_seg);
    td_mid.is_in_tree = false;
    // Check if the vertex that is removed lies between the two tree edges.
    if (seg_above != tree->end()) {
        if (!on_right_side(cur_vt, *seg_above, false)) {
	    // found conflicts prev_seg - seg_above and mid_seg - seg_above
#if defined(CGAL_POLY_GENERATOR_DEBUG)
            std::cout << "conflict2 is seg_above" << std::endl;
#endif
            conflict1 = prev_vt;
	    conflict2 = *seg_above;
	    return false;
        }
    }
    if (seg_above != tree->begin()) {
        --seg_above;  // which turns it into seg_below
        if (!on_right_side(cur_vt, *seg_above, true)) {
	    // found conflicts prev_seg - seg_below and mid_seg - seg_below
#if defined(CGAL_POLY_GENERATOR_DEBUG)
            std::cout << "conflict2 is --seg_above" << std::endl;
#endif
            conflict1 = prev_vt;
	    conflict2 = *seg_above;
	    return false;
        }
    }
    return true;
}

template <class ForwardIterator, class PolygonTraits>
void Vertex_data<ForwardIterator, PolygonTraits>::
sweep(Tree *tree)
{
    if (this->m_size < 3)
    	return;
    bool success = true;
    for (Index_t i=0; i< this->m_size; ++i) {
        Vertex_index cur = index_at_rank(Vertex_order(i));
	    Vertex_index prev_vt = prev(cur), next_vt = next(cur);
	    if (ordered_left_to_right(cur, next_vt)) {
	        if (ordered_left_to_right(cur, prev_vt))
	            success = insertion_event(tree, prev_vt, cur, next_vt);
	        else
	            success = replacement_event(tree, prev_vt, cur);
	    } else {
	        if (ordered_left_to_right(cur, prev_vt))
	            success = replacement_event(tree, cur, prev_vt);
	        else
	            success = deletion_event(tree, prev_vt, cur);
	    }
	    if (!success)
	        break;
    }
    if (!success)
    	this->is_simple_result = false;
}
}
// ----- End of implementation of i_generator_polygon functions. -----

template <class Iterator, class PolygonTraits>
std::pair<std::ptrdiff_t,std::ptrdiff_t>
check_simple_polygon(Iterator points_begin, Iterator points_end,
                     const PolygonTraits& polygon_traits)
{
    typedef Iterator ForwardIterator;
    typedef std::set<i_generator_polygon::Vertex_index,
      i_generator_polygon::Less_segments<ForwardIterator,PolygonTraits> > Tree;
    i_generator_polygon::Vertex_data<ForwardIterator, PolygonTraits>
        vertex_data(points_begin, points_end, polygon_traits);
    Tree tree(&vertex_data);
    vertex_data.init(&tree);
    vertex_data.sweep(&tree);
    std::pair<std::ptrdiff_t, std::ptrdiff_t> result;
    if (vertex_data.is_simple_result) {
        result.first = result.second = -1;
	return result;
    }
    // swap with vertex_data.conflict1, vertex_data.conflict2;
    if (vertex_data.conflict1.as_int() < vertex_data.conflict2.as_int()) {
	result.first = vertex_data.conflict1.as_int();
	result.second = vertex_data.conflict2.as_int();
    } else {
        result.first = vertex_data.conflict2.as_int();
	result.second = vertex_data.conflict1.as_int();
    }
    return result;
}

template <class Iterator, class PolygonTraits>
void make_simple_polygon(Iterator points_begin, Iterator points_end,
                         const PolygonTraits& polygon_traits)
{
  std::pair<std::ptrdiff_t,std::ptrdiff_t> swap_interval;

#if defined(CGAL_POLY_GENERATOR_DEBUG)
    Iterator it;
    std::cout << "In make_simple_polygon the points are: " << std::endl;
    int size = 0;
    for (it = points_begin; it != points_end; it++, size++)
      std::cout << *it << " ";
    std::cout << std::endl;
#endif


    do {
        swap_interval = check_simple_polygon(points_begin,
	                                     points_end, polygon_traits);
#if defined(CGAL_POLY_GENERATOR_DEBUG)
        std::cout << swap_interval.first << " "
                  << swap_interval.second << std::endl; 
        CGAL_polygon_assertion(swap_interval.first >= -1 && 
                               swap_interval.second >= -1 &&
                               swap_interval.first < size &&
                               swap_interval.second < size);
#endif
        // will break out when a negative nonsense value is selected
        // or with -1 is given indicating that the polygon was simple.
        // For positive nonsense values, one needs to know the how
        // many points there are...
	if (swap_interval.first <= -1 || swap_interval.second <= -1)
	    break;
        // swap with vertex_data.conflict1, vertex_data.conflict2;
	Iterator b = points_begin;
	std::advance(b, swap_interval.first+1);
	Iterator e = b;
	std::advance(e, swap_interval.second-swap_interval.first);
        std::reverse(b, e);
    } while (true);
}

} // end of namespace CGAL

#endif // CGAL_RANDOM_POLYGON_2_SWEEP_H
