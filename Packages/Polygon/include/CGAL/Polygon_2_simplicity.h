// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  :
//
// file          : include/CGAL/Polygon_2_simplicity.h
// source        :
// author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_POLYGON_2_SIMPLICITY_H
#define CGAL_POLYGON_2_SIMPLICITY_H

#include <CGAL/enum.h>
#include <CGAL/polygon_assertions.h>
#include <set>
#include <vector>
#include <algorithm>

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

namespace i_polygon {
 // namespace CGAL::i_polygon is used for internal functions

typedef std::vector<int>::size_type Index_t;

struct Vertex_index {
    Vertex_index() {}
    explicit Vertex_index(Index_t i): m_i(i) {}
    Index_t as_int() const {return m_i;}
    Vertex_index operator++() {++m_i; return *this; }
private:
    Index_t m_i;
};

struct Vertex_order {
    explicit Vertex_order(Index_t i): m_i(i) {}
    Index_t as_int() {return m_i;}
private:
    Index_t m_i;
};

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

// The data in Edge_data is attached to an edge when it is (about to be)
// inserted in the tree.
// Although conceptually this data belongs in the tree, it is stored with
// the vertices in the Vertex_data structure.

template <class ForwardIterator, class PolygonTraits>
struct Edge_data {
    typedef std::set<Vertex_index,
            Less_segments<ForwardIterator,PolygonTraits> > Tree;
    Edge_data() : is_in_tree(false) {}
    typename Tree::iterator tree_it; // The iterator of the edge in the tree.
                                     // Needed for cross reference. If edge j
				     // is in the tree: *edges[j].tree_it == j
    bool is_in_tree :1;              // Must be set -after- inserting the edge
                                     // in the tree. Plays a role in the
				     // comparison function of the tree.
    bool is_left_to_right :1;        // Direction of edge from vertex v to v+1
};

template <class ForwardIterator, class PolygonTraits>
class Vertex_data {
public:
    typedef std::set<Vertex_index,
                    Less_segments<ForwardIterator,PolygonTraits> > Tree;

    typedef typename PolygonTraits::Point_2 Point_2;
    Vertex_data(ForwardIterator begin, ForwardIterator end,PolygonTraits pgnt);
//    ForwardIterator points_start;
    std::vector<ForwardIterator> iterators;
    std::vector<Vertex_order> m_order_of;
    std::vector<Vertex_index> m_idx_at_rank;
    std::vector<Edge_data<ForwardIterator, PolygonTraits> > edges;
    std::vector<Vertex_index>::size_type m_size;
    typename PolygonTraits::Orientation_2 orientation_2;
    typename PolygonTraits::Less_xy_2 less_xy_2;
    bool is_simple_result;

    bool ordered_left_to_right(Vertex_index v1, Vertex_index v2)
        { return  m_order_of[v1.as_int()].as_int() <
	m_order_of[v2.as_int()].as_int();}
    Vertex_index index_at_rank(Vertex_order vo) const
        { return m_idx_at_rank[vo.as_int()];}
    Vertex_index next(Vertex_index k) const
        { ++k; return k.as_int() == m_size ? Vertex_index(0) : k;}
    Vertex_index prev(Vertex_index k) const
        { return k.as_int() == 0
	       ?  Vertex_index(m_size-1)
	       : Vertex_index(k.as_int()-1);
	}
    void left_and_right_index(Vertex_index &left, Vertex_index &right,
            Vertex_index edge);
    Vertex_index left_index(Vertex_index edge)
        { return edges[edge.as_int()].is_left_to_right ? edge : next(edge); }
    Point_2 point(Vertex_index i)
//    { return points_start[i.as_int()];}
        { return *iterators[i.as_int()];}
    void sweep(Tree *tree);
    bool insertion_event(Tree *tree,
                Vertex_index i, Vertex_index j, Vertex_index k);
    bool replacement_event(Tree *tree,
                Vertex_index cur, Vertex_index to_insert);
    bool deletion_event(Tree *tree, Vertex_index i, Vertex_index j);
    bool on_right_side(Vertex_index vt, Vertex_index edge, bool above);
};

template <class ForwardIterator, class PolygonTraits>
class Less_vertex_data {
    Vertex_data<ForwardIterator, PolygonTraits> *m_vertex_data;
public:
    Less_vertex_data(Vertex_data<ForwardIterator, PolygonTraits> *vd)
    : m_vertex_data(vd) {}
    bool operator()(Vertex_index i, Vertex_index j);
};

} // end of namespace i_polygon

// ----- implementation of i_polygon functions. -----

namespace i_polygon {
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
      case LEFTTURN: return true;
      case RIGHTTURN: return false;
      case COLLINEAR: break;
    }
    m_vertex_data->is_simple_result = false;
    return true;
}

template <class ForwardIterator, class PolygonTraits>
bool Less_vertex_data<ForwardIterator, PolygonTraits>::
operator()(Vertex_index i, Vertex_index j)
{
    return m_vertex_data->less_xy_2(
            m_vertex_data->point(i), m_vertex_data->point(j));
}

template <class ForwardIterator, class PolygonTraits>
Vertex_data<ForwardIterator, PolygonTraits>::
Vertex_data(ForwardIterator begin, ForwardIterator end,
            PolygonTraits pgn_traits)
: // points_start(begin),
  orientation_2(pgn_traits.orientation_2_object()),
  less_xy_2(pgn_traits.less_xy_2_object())
{
    m_size = std::distance(begin, end);
    is_simple_result = true;
    m_idx_at_rank.reserve(m_size);
    iterators.reserve(m_size);
    m_order_of.insert(m_order_of.end(), m_size, Vertex_order(0));
    edges.insert(edges.end(), m_size,
            Edge_data<ForwardIterator, PolygonTraits>());
    for (Index_t i = 0; i< m_size; ++i, ++begin) {
        m_idx_at_rank.push_back(Vertex_index(i));
	iterators.push_back(begin);
    }
    std::sort(m_idx_at_rank.begin(), m_idx_at_rank.end(),
        Less_vertex_data<ForwardIterator, PolygonTraits>(this));
    for (Index_t j = 0; j < m_size; ++j) {
	Vertex_order vo(j);
        m_order_of[index_at_rank(vo).as_int()] = vo;
    }
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
    // check which endpoint is above the other
    bool left_turn;
    switch(orientation_2(point(prev_vt), point(mid_vt), point(next_vt))) {
      case LEFTTURN: left_turn = true; break;
      case RIGHTTURN: left_turn = false; break;
      case COLLINEAR: return false;
      
    }
    Edge_data<ForwardIterator, PolygonTraits>
        &td_prev = edges[prev_vt.as_int()],
        &td_mid = edges[mid_vt.as_int()];
    td_prev.is_in_tree = false;
    td_prev.is_left_to_right = false;
    td_mid.is_in_tree = false;
    td_mid.is_left_to_right = true;
    // insert the highest chain first
    std::pair<CGAL_TYPENAME_MSVC_NULL Tree::iterator, bool> result;
    if (left_turn) {
        result = tree->insert(prev_vt);
	// assert(result.second)
	td_prev.tree_it = result.first;
        td_prev.is_in_tree = true;
        result = tree->insert(mid_vt);
	// assert(result.second)
	td_mid.tree_it = result.first;
        td_mid.is_in_tree = true;
    } else {
        result = tree->insert(mid_vt);
	// assert(result.second)
	td_mid.tree_it = result.first;
        td_mid.is_in_tree = true;
        result = tree->insert(prev_vt);
	// assert(result.second)
	td_prev.tree_it = result.first;
        td_prev.is_in_tree = true;
    }
    return true;
}

template <class ForwardIterator, class PolygonTraits>
bool Vertex_data<ForwardIterator, PolygonTraits>::
on_right_side(Vertex_index vt, Vertex_index edge_id, bool above)
{
    Orientation turn =
        orientation_2(point(edge_id), point(vt), point(next(edge_id)));
    bool leftturn = edges[edge_id.as_int()].is_left_to_right ? above : !above;
    if (leftturn) {
        if (turn != RIGHTTURN) {
            return false;
        }
    } else {
        if (turn != LEFTTURN) {
            return false;
        }
    }
    return true;
}

template <class ForwardIterator, class PolygonTraits>
bool Vertex_data<ForwardIterator, PolygonTraits>::
replacement_event(Tree *tree, Vertex_index cur_edge, Vertex_index next_edge)
{
    // check if continuation point is on the right side of neighbor segments
    typedef typename Tree::iterator It;
    Edge_data<ForwardIterator, PolygonTraits> &td = edges[cur_edge.as_int()];
    CGAL_polygon_assertion(td.is_in_tree);
    It cur_seg = td.tree_it;
    Vertex_index cur_vt = (td.is_left_to_right) ? next_edge : cur_edge;
    if (cur_seg != tree->begin()) {
        It seg_below = cur_seg;
	--seg_below;
	if (!on_right_side(cur_vt, *seg_below, true)) {
	    return false;
        }
    }
    It seg_above = cur_seg;
    ++ seg_above;
    if (seg_above != tree->end()) {
        if (!on_right_side(cur_vt, *seg_above, false)) {
	    return false;
        }
    }
    // replace the segment
    Edge_data<ForwardIterator, PolygonTraits> &new_td =
            edges[next_edge.as_int()];
    new_td.is_left_to_right = td.is_left_to_right;
    new_td.is_in_tree = false;
    tree->erase(cur_seg);
    td.is_in_tree = false;
    new_td.tree_it = tree->insert(seg_above, next_edge);
    new_td.is_in_tree = true;
    return true;
}

template <class ForwardIterator, class PolygonTraits>
bool Vertex_data<ForwardIterator, PolygonTraits>::
deletion_event(Tree *tree, Vertex_index prev_vt, Vertex_index mid_vt)
{
    // check if continuation point is on the right side of neighbor segments
    typedef typename Tree::iterator It;
    Edge_data<ForwardIterator, PolygonTraits>
        &td_prev = edges[prev_vt.as_int()],
        &td_mid = edges[mid_vt.as_int()];
    It prev_seg = td_prev.tree_it, mid_seg = td_mid.tree_it;
    Vertex_index cur_vt = (td_prev.is_left_to_right) ? mid_vt : prev_vt;
    It seg_above = prev_seg;
    ++seg_above;
    if (seg_above == mid_seg) ++seg_above;
    tree->erase(prev_seg);
    td_prev.is_in_tree = false;
    tree->erase(mid_seg);
    td_mid.is_in_tree = false;
    if (seg_above != tree->end()) {
        if (!on_right_side(cur_vt, *seg_above, false))
	    return false;
    }
    // remove the segments
    if (seg_above != tree->begin()) {
        --seg_above; // which turns it in seg_below
        if (!on_right_side(cur_vt, *seg_above, true))
	    return false;
    }
    return true;
}

template <class ForwardIterator, class PolygonTraits>
void Vertex_data<ForwardIterator, PolygonTraits>::
sweep(Tree *tree)
{
    if (m_size < 3)
    	return;
    bool succes = true;
    for (Index_t i=0; i< m_size; ++i) {
        Vertex_index cur = index_at_rank(Vertex_order(i));
	    Vertex_index prev_vt = prev(cur), next_vt = next(cur);
	    if (ordered_left_to_right(cur, next_vt)) {
	        if (ordered_left_to_right(cur, prev_vt))
	            succes = insertion_event(tree, prev_vt, cur, next_vt);
	        else
	            succes = replacement_event(tree, prev_vt, cur);
	    } else {
	        if (ordered_left_to_right(cur, prev_vt))
	            succes = replacement_event(tree, cur, prev_vt);
	        else
	            succes = deletion_event(tree, prev_vt, cur);
	    }
	    if (!succes)
	        break;
    }
    if (!succes)
    	is_simple_result = false;
}
}
// ----- End of implementation of i_polygon functions. -----


template <class Iterator, class PolygonTraits>
bool is_simple_polygon(Iterator points_begin, Iterator points_end,
        PolygonTraits polygon_traits)
{
    typedef Iterator ForwardIterator;
    typedef std::set<i_polygon::Vertex_index,
            i_polygon::Less_segments<ForwardIterator,PolygonTraits> > Tree;
    i_polygon::Vertex_data<ForwardIterator, PolygonTraits>
        vertex_data(points_begin, points_end, polygon_traits);
    Tree tree(&vertex_data);
    vertex_data.sweep(&tree);
    return vertex_data.is_simple_result;
}

} // end of namespace CGAL

#endif
