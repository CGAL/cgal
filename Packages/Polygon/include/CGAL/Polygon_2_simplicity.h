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

// #define GJ_DEBUG_43

namespace CGAL {

namespace i_polygon {

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

template <class RandomAccessIt, class PolygonTraits>
class Vertex_data ;

template <class RandomAccessIt, class PolygonTraits>
class Less_segments {
    Vertex_data<RandomAccessIt, PolygonTraits> *m_vertex_data;
    bool less_than_in_tree(Vertex_index i, Vertex_index j);
  public:
    Less_segments(Vertex_data<RandomAccessIt, PolygonTraits> *vertex_data)
    : m_vertex_data(vertex_data) {}
    bool operator()(Vertex_index i, Vertex_index j);
};



template <class RandomAccessIt, class PolygonTraits>
struct Edge_data {
    typedef std::set<Vertex_index, Less_segments<RandomAccessIt,PolygonTraits> >
            Tree;
    Edge_data() : is_in_tree(false) {}
    typename Tree::iterator tree_it;
    bool is_in_tree :1;
    bool is_in_order :1;
};

template <class RandomAccessIt, class PolygonTraits>
class Vertex_data {
public:
    typedef std::set<Vertex_index, Less_segments<RandomAccessIt,PolygonTraits> >
	            Tree;

    typedef typename PolygonTraits::Point_2 Point_2;
    Vertex_data(RandomAccessIt begin, RandomAccessIt end, PolygonTraits pgnt);
    RandomAccessIt points_start;
    std::vector<Vertex_order> m_order_of;
    std::vector<Vertex_index> m_idx_at_rank;
    std::vector<Edge_data<RandomAccessIt, PolygonTraits> > edges;
    std::vector<Vertex_index>::size_type m_size;
    typename PolygonTraits::Orientation_2 orientation_2;
    typename PolygonTraits::Less_xy_2 less_xy_2;
    bool is_simple_result;

    Vertex_order xy_order_of(Vertex_index vi) const
    	{ return m_order_of[vi.as_int()];}
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
        { return edges[edge.as_int()].is_in_order ? edge : next(edge); }
    Point_2 point(Vertex_index i) { return points_start[i.as_int()];}
    void sweep(Tree *tree);
    bool chain_start(Tree *tree, Vertex_index i, Vertex_index j, Vertex_index k);
    bool chain_continuation(Tree *tree, Vertex_index cur, Vertex_index to_insert);
    bool chain_end(Tree *tree, Vertex_index i, Vertex_index j);
    bool on_right_side(Vertex_index vt, Vertex_index edge, bool above);
#ifdef GJ_DEBUG_43
   void print_tree(Tree *tree);
#endif
};

#ifdef GJ_DEBUG_43
template <class RandomAccessIt, class PolygonTraits>
void Vertex_data<RandomAccessIt, PolygonTraits>::
print_tree(Tree *tree)
{
    typedef typename Tree::iterator Tree_it;
    for (Tree_it cur = tree->begin(); cur != tree->end(); ++cur) {
	Vertex_index nb = next(*cur);
        std::cout << (*cur).as_int() << ' ' << nb.as_int() <<'\n';
    }
    std::cout << "-----\n";
}
#endif

template <class RandomAccessIt, class PolygonTraits>
class Less_vertex_data {
    Vertex_data<RandomAccessIt, PolygonTraits> *m_vertex_data;
public:
    Less_vertex_data(Vertex_data<RandomAccessIt, PolygonTraits> *vd)
    : m_vertex_data(vd) {}
    bool operator()(Vertex_index i, Vertex_index j);
};

} // end of namespace i_polygon

// ----- implementation of i_polygon functions. -----

template <class RandomAccessIt, class PolygonTraits>
bool i_polygon::Less_segments<RandomAccessIt, PolygonTraits>::
operator()(Vertex_index i, Vertex_index j)
{
    if (m_vertex_data->edges[j.as_int()].is_in_tree) {
        return less_than_in_tree(i,j);
    } else {
        return !less_than_in_tree(j,i);
    }
}

template <class RandomAccessIt, class PolygonTraits>
bool i_polygon::Less_segments<RandomAccessIt, PolygonTraits>::
less_than_in_tree(Vertex_index new_edge, Vertex_index tree_edge)
{
    CGAL_polygon_precondition(m_vertex_data->edges[tree_edge.as_int()].is_in_tree);
    CGAL_polygon_precondition(!m_vertex_data->edges[new_edge.as_int()].is_in_tree);
    Vertex_index left, mid, right;
    m_vertex_data->left_and_right_index(left, right, tree_edge);
    mid = m_vertex_data->left_index(new_edge);
#ifdef GJ_DEBUG_43
if (new_edge.as_int() == 14 || new_edge.as_int() == 15)
        std::cout << "Checking  "<< new_edge.as_int() << " and " 
	    << tree_edge.as_int() <<'\n';
    std::cout<< "left: "<<left.as_int()<<"  right: "<<right.as_int()<<'\n';
#endif
    if (mid.as_int() == left.as_int()) {
#ifdef GJ_DEBUG_43
        std::cout << "Insertion detected. "<< new_edge.as_int() << "<" 
	    << tree_edge.as_int() <<'\n';
#endif
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

template <class RandomAccessIt, class PolygonTraits>
bool i_polygon::Less_vertex_data<RandomAccessIt, PolygonTraits>::
operator()(Vertex_index i, Vertex_index j)
{
    return m_vertex_data->less_xy_2(
            m_vertex_data->point(i), m_vertex_data->point(j));
}

template <class RandomAccessIt, class PolygonTraits>
i_polygon::Vertex_data<RandomAccessIt, PolygonTraits>::
Vertex_data(RandomAccessIt begin, RandomAccessIt end, PolygonTraits pgn_traits)
: points_start(begin),
  orientation_2(pgn_traits.orientation_2_object()),
  less_xy_2(pgn_traits.less_xy_2_object())
{
    m_size = end - begin;
    is_simple_result = true;
    m_idx_at_rank.reserve(m_size);
    m_order_of.insert(m_order_of.end(), m_size, Vertex_order(0));
    edges.insert(edges.end(), m_size, Edge_data<RandomAccessIt, PolygonTraits>());
    for (Index_t i = 0; i< m_size; ++i)
        m_idx_at_rank.push_back(Vertex_index(i));
    std::sort(m_idx_at_rank.begin(), m_idx_at_rank.end(),
        Less_vertex_data<RandomAccessIt, PolygonTraits>(this));
    for (Index_t j = 0; j < m_size; ++j) {
	Vertex_order vo(j);
        m_order_of[index_at_rank(vo).as_int()] = vo;
    }
}

template <class RandomAccessIt, class PolygonTraits>
void i_polygon::Vertex_data<RandomAccessIt, PolygonTraits>::
left_and_right_index(Vertex_index &left, Vertex_index &right,
            Vertex_index edge)
{
    if (edges[edge.as_int()].is_in_order) {
        left = edge; right = next(edge);
    } else {
        right = edge; left = next(edge);
    }
}

template <class RandomAccessIt, class PolygonTraits>
bool i_polygon::Vertex_data<RandomAccessIt, PolygonTraits>::
chain_start(Tree *tree, Vertex_index prev_vt,
            Vertex_index mid_vt, Vertex_index next_vt)
{
    // check which endpoint is above the other
    bool left_turn;
    switch(orientation_2(point(prev_vt), point(mid_vt), point(next_vt))) {
      case LEFTTURN: left_turn = true; break;
      case RIGHTTURN: left_turn = false; break;
      case COLLINEAR: return false;
      
    }
    Edge_data<RandomAccessIt, PolygonTraits>
        &td_prev = edges[prev_vt.as_int()],
        &td_mid = edges[mid_vt.as_int()];
    td_prev.is_in_tree = false;
    td_prev.is_in_order = false;
    td_mid.is_in_tree = false;
    td_mid.is_in_order = true;
    // insert the highest chain first
#ifdef GJ_DEBUG_43
if (left_turn)
  std::cout << "Ins " << prev_vt.as_int() << " and " << mid_vt.as_int() <<'\n';
else
  std::cout << "Ins " << mid_vt.as_int() << " and " << prev_vt.as_int() <<'\n';
#endif
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

template <class RandomAccessIt, class PolygonTraits>
bool i_polygon::Vertex_data<RandomAccessIt, PolygonTraits>::
on_right_side(Vertex_index vt, Vertex_index edge_id, bool above)
{
    Orientation turn =
        orientation_2(point(edge_id), point(vt), point(next(edge_id)));
    bool leftturn = edges[edge_id.as_int()].is_in_order ? above : !above;
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

template <class RandomAccessIt, class PolygonTraits>
bool i_polygon::Vertex_data<RandomAccessIt, PolygonTraits>::
chain_continuation(Tree *tree, Vertex_index cur_edge, Vertex_index next_edge)
{
    // check if continuation point is on the right side of neighbor segments
    typedef typename Tree::iterator It;
    Edge_data<RandomAccessIt, PolygonTraits> &td = edges[cur_edge.as_int()];
    CGAL_polygon_assertion(td.is_in_tree);
    It cur_seg = td.tree_it;
    Vertex_index cur_vt = (td.is_in_order) ? next_edge : cur_edge;
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
    Edge_data<RandomAccessIt, PolygonTraits> &new_td =
            edges[next_edge.as_int()];
    new_td.is_in_order = td.is_in_order;
    new_td.is_in_tree = false;
    tree->erase(cur_seg);
    td.is_in_tree = false;
    new_td.tree_it = tree->insert(seg_above, next_edge);
    new_td.is_in_tree = true;
    return true;
}

template <class RandomAccessIt, class PolygonTraits>
bool i_polygon::Vertex_data<RandomAccessIt, PolygonTraits>::
chain_end(Tree *tree, Vertex_index prev_vt, Vertex_index mid_vt)
{
    // check if continuation point is on the right side of neighbor segments
    typedef typename Tree::iterator It;
    Edge_data<RandomAccessIt, PolygonTraits>
        &td_prev = edges[prev_vt.as_int()],
        &td_mid = edges[mid_vt.as_int()];
    It prev_seg = td_prev.tree_it, mid_seg = td_mid.tree_it;
    Vertex_index cur_vt = (td_prev.is_in_order) ? mid_vt : prev_vt;
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

template <class RandomAccessIt, class PolygonTraits>
void i_polygon::Vertex_data<RandomAccessIt, PolygonTraits>::
sweep(Tree *tree)
{
    if (m_size < 3)
    	return;
    bool succes = true;
    for (Index_t i=0; i< m_size; ++i) {
        Vertex_index cur = index_at_rank(Vertex_order(i));
	Vertex_index prev_vt = prev(cur), next_vt = next(cur);
	if (xy_order_of(next_vt).as_int() > i) {
	    if (xy_order_of(prev_vt).as_int() > i)
	        succes = chain_start(tree, prev_vt, cur, next_vt);
	    else
	        succes = chain_continuation(tree, prev_vt, cur);
	} else {
	    if (xy_order_of(prev_vt).as_int() > i)
	        succes = chain_continuation(tree, cur, prev_vt);
	    else
	        succes = chain_end(tree, prev_vt, cur);
	}
#ifdef GJ_DEBUG_43
  std::cout << "after treating " << cur.as_int() << ":\n";
  print_tree(tree);
#endif
	if (!succes)
	    break;
    }
    if (!succes)
    	is_simple_result = false;
}

// ----- End of implementation of i_polygon functions. -----


template <class Iterator, class PolygonTraits>
bool is_simple_polygon(Iterator points_begin, Iterator points_end,
        PolygonTraits polygon_traits)
{
    typedef Iterator RandomAccessIt;
    typedef std::set<i_polygon::Vertex_index,
            i_polygon::Less_segments<RandomAccessIt,PolygonTraits> > Tree;
    i_polygon::Vertex_data<RandomAccessIt, PolygonTraits>
        vertex_data(points_begin, points_end, polygon_traits);
    Tree tree(&vertex_data);
    vertex_data.sweep(&tree);
    return vertex_data.is_simple_result;
}

} // end of namespace CGAL

#endif
