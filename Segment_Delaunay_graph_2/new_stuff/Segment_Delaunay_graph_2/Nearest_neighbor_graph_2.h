// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_NEAREST_NEIGHBOR_GRAPH_2_H
#define CGAL_NEAREST_NEIGHBOR_GRAPH_2_H 1

#include <CGAL/basic.h>
#include <set>
#include <map>

#include "Directed_graph.h"

namespace CGAL {

template<class NIT, class NH>
class Node_iterator_wrapper
{
private:
  typedef NIT   Node_iterator;
  typedef NH    Node_handle;

  typedef Node_iterator_wrapper<Node_iterator,Node_handle>  Self;

  typedef std::iterator_traits<Node_iterator>  Traits;

public:
  typedef typename std::size_t                 size_type;
  typedef typename Traits::difference_type     difference_type;
  typedef typename Traits::iterator_category   iterator_category;
  typedef typename Traits::pointer             pointer;
  typedef typename Traits::reference           reference;
  typedef typename Traits::value_type          value_type;

public:
  Node_iterator_wrapper() : cur_() {}
  Node_iterator_wrapper(Node_iterator nit) : cur_(nit) {}

  operator Node_handle() const {
    return *cur_;
  }

  reference operator*() const {
    return *cur_;
  }

  pointer operator->() const {
    return cur_.operator->();
  }

  Self& operator++() {
    ++cur_;
    return *this;
  }

  Self& operator--() {
    --cur_;
    return *this;
  }

  Self operator++(int) {
    return cur_++;
  }

  Self operator--(int) {
    return cur_--;
  }

  bool operator==(const Self& other) const {
    return cur_ == other.cur_;
  }

  bool operator!=(const Self& other) const {
    return cur_ != other.cur_;
  }

private:
  Node_iterator cur_;
};



template<class DG, class TR>
class Nearest_neighbor_graph_2
  : private Directed_graph<typename DG::Vertex_handle>
{
public:
  typedef  DG   Delaunay_graph;
  typedef  TR   Traits;

  typedef typename Delaunay_graph::Vertex_handle   Vertex_handle;

private:
  typedef Directed_graph<Vertex_handle>            Base;

  typedef typename Delaunay_graph::Finite_vertices_iterator
  Finite_vertices_iterator;

  typedef typename Delaunay_graph::Vertex_circulator
  Vertex_circulator;

public:
  typedef typename Base::Node_handle     Node_handle;
  typedef typename Base::Node_iterator   Node_iterator_base;

  typedef Node_iterator_wrapper<Node_iterator_base,Node_handle>
  Node_iterator;

  typedef typename Base::size_type       size_type;

private:
  void initialize_nng(const Delaunay_graph& dg)
  {
    CGAL_precondition( dg.number_of_vertices() >= 2 );

    std::map<Vertex_handle,Node_handle> v2n_map;

    for (Finite_vertices_iterator vit = dg.finite_vertices_begin();
	 vit != dg.finite_vertices_end(); ++vit) {
      Vertex_handle v(vit);
      Node_handle n = this->add_node(v);
      v2n_map[v] = n;
    }

    typename Traits::Distance_comparator_2 comparator =
      tr_.distance_comparator_2_object();

    typename Traits::Access_site_2 accessor =
      tr_.access_site_2_object();

    for (Finite_vertices_iterator vit = dg.finite_vertices_begin();
	 vit != dg.finite_vertices_end(); ++vit) {
      Vertex_handle v(vit);
      typename Traits::Site_2 t = accessor(v);

      //      CGAL_assertion( t.is_point() && t.is_input() );

      std::set<Vertex_handle> nearest_neighbors;

      Vertex_circulator vc_start = dg.incident_vertices(v);
      Vertex_circulator vc = vc_start;

      while ( dg.is_infinite(vc) ) { ++vc; }
      Vertex_handle nn = vc;
      nearest_neighbors.clear();
      nearest_neighbors.insert(nn);
      ++vc;
      do {
	if ( !dg.is_infinite(vc) ) {
	  CGAL::Comparison_result cr =
	    comparator(t, accessor(nn), accessor(vc));
	  if ( cr == CGAL::EQUAL ) {
	    nearest_neighbors.insert(vc);
	  } else if ( cr == CGAL::LARGER ) {
	    nn = vc;
	    nearest_neighbors.clear();
	    nearest_neighbors.insert(nn);
	  }
	}
	++vc;
      } while (vc != vc_start);

      for (typename std::set<Vertex_handle>::iterator
	     it = nearest_neighbors.begin();
	   it != nearest_neighbors.end(); ++it) {
	this->add_directed_edge(v2n_map[v], v2n_map[*it]);
      }
      // end of while-loop for incident vertices
    } // end of for-loop for all vertices

  }

public:
  Nearest_neighbor_graph_2(const Delaunay_graph& dg)
    : tr_(dg.geom_traits())
  {
    initialize_nng(dg);
  }

  Nearest_neighbor_graph_2(const Delaunay_graph& dg, const Traits& tr)
    : tr_(tr)
  {
    initialize_nng(dg);
  }

  inline const Traits& traits() const { return tr_; }

  inline Node_iterator nodes_begin() const {
    return Base::nodes_begin();
  }

  inline Node_iterator nodes_end() const {
    return Base::nodes_end();
  }

  inline Node_iterator out_neighbors_begin(Node_handle n) const {
    return n->out_neighbors_begin();
  }

  inline Node_iterator out_neighbors_end(Node_handle n) const {
    return n->out_neighbors_end();
  }

  inline Node_iterator in_neighbors_begin(Node_handle n) const {
    return n->in_neighbors_begin();
  }

  inline Node_iterator in_neighbors_end(Node_handle n) const {
    return n->in_neighbors_end();
  }

  inline size_type in_degree(Node_handle n) const {
    return n->in_degree();
  }

  inline size_type out_degree(Node_handle n) const {
    return n->out_degree();
  }

  inline Vertex_handle operator[](Node_handle n) const {
    return n->data();
  }

private:
  Traits tr_;
};

} //namespace CGAL

#endif // CGAL_NEAREST_NEIGHBOR_GRAPH_2_H
