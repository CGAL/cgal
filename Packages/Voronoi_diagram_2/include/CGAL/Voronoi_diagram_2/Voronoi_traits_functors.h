// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_VORONOI_TRAITS_FUNCTORS_H
#define CGAL_VORONOI_DIAGRAM_2_VORONOI_TRAITS_FUNCTORS_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

struct Null_functor
{
  Null_functor() {}
  template<typename T> Null_functor(T t) {}
};

//=========================================================================

template<class Functor>
struct Functor_exists
{
  typedef Tag_true  Value;
};

template<>
struct Functor_exists<Null_functor>
{
  typedef Tag_false Value;
};

template<class VT, class SI> class Default_caching_site_inserter;

template<class VT>
struct Functor_exists< Default_caching_site_inserter<VT,Null_functor> >
{
  typedef Tag_false Value;
};

//=========================================================================
//=========================================================================

struct Degeneracy_tester_base
{
  inline void clear() {}
  inline void swap(Degeneracy_tester_base&) {}
  inline bool is_valid() const { return true; }
};

//=========================================================================
//=========================================================================

template<class DG>
struct Identity_edge_degeneracy_tester
  : public Degeneracy_tester_base
{
  typedef DG                  Delaunay_graph;
  typedef bool                result_type;
  typedef Arity_tag<2>        Arity;

  typedef typename Delaunay_graph::Edge                   Edge;
  typedef typename Delaunay_graph::Face_handle            Face_handle;
  typedef typename Delaunay_graph::Edge_circulator        Edge_circulator;
  typedef typename Delaunay_graph::All_edges_iterator     All_edges_iterator;
  typedef typename Delaunay_graph::Finite_edges_iterator  Finite_edges_iterator;

  bool operator()(const Delaunay_graph& dual,
		  const Face_handle& f, int i) const {
    return false;
  }

  bool operator()(const Delaunay_graph& dual, const Edge& e) const {
    return false;
  }

  bool operator()(const Delaunay_graph& dual,
		  const All_edges_iterator& eit) const {
    return false;
  }

  bool operator()(const Delaunay_graph& dual,
		  const Finite_edges_iterator& eit) const {
    return false;
  }

  bool operator()(const Delaunay_graph& dual,
		  const Edge_circulator& ec) const {
    return false;
  }
};

//=========================================================================
//=========================================================================

template<class DG>
struct Identity_face_degeneracy_tester
  : public Degeneracy_tester_base
{
  typedef DG                  Delaunay_graph;
  typedef bool                result_type;
  typedef Arity_tag<2>        Arity;

  typedef typename Delaunay_graph::Vertex_handle  Vertex_handle;

  bool operator()(const Delaunay_graph&, const Vertex_handle&) const {
    return false;
  }
};

//=========================================================================
//=========================================================================


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_VORONOI_TRAITS_FUNCTORS_H
