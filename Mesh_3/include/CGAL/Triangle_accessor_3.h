// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_TRIANGLE_ACCESSOR_3_H
#define CGAL_TRIANGLE_ACCESSOR_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh/Surface_mesh_fwd.h>
#include <CGAL/boost/graph/Graph_with_descriptor_with_graph.h>

namespace CGAL {

template <typename Polyhedron, typename K>
class Triangle_accessor_3 {
  // This class template should never be instantiated, but only its
  // specializations.
  typedef typename Polyhedron::Error_bad_match Error_bad_match;
};


template < class K,class Items,
           template < class T, class I, class A>
           class T_HDS,
           class Alloc>
class Triangle_accessor_3<Polyhedron_3<K,Items,T_HDS,Alloc>, K >
{
  typedef Polyhedron_3<K,Items,T_HDS,Alloc> Polyhedron;
public:
  typedef typename K::Triangle_3                    Triangle_3;
  typedef typename Polyhedron::Facet_const_iterator Triangle_iterator;
  typedef typename Polyhedron::Facet_const_handle   Triangle_handle;

  Triangle_accessor_3() { }

  Triangle_iterator triangles_begin(const Polyhedron& p) const
  {
    return p.facets_begin();
  }

  Triangle_iterator triangles_end(const Polyhedron& p) const
  {
    return p.facets_end();
  }

  Triangle_3 triangle(const Triangle_handle& handle) const
  {
    typedef typename K::Point_3 Point;
    const Point& a = handle->halfedge()->vertex()->point();
    const Point& b = handle->halfedge()->next()->vertex()->point();
    const Point& c = handle->halfedge()->next()->next()->vertex()->point();
    return Triangle_3(a,b,c);
  }
};


  template <class P, class K>
class Triangle_accessor_3<Graph_with_descriptor_with_graph<Surface_mesh<P> >, K >
{
  typedef  Graph_with_descriptor_with_graph<Surface_mesh<P> > Polyhedron;
public:
  typedef typename K::Triangle_3                    Triangle_3;
  typedef typename boost::graph_traits<Polyhedron>::face_iterator Triangle_iterator;
  typedef typename boost::graph_traits<Polyhedron>::face_iterator   Triangle_handle;

  Triangle_accessor_3() { }

  Triangle_iterator triangles_begin(const Polyhedron& p) const
  {
    return faces(p).first;
  }

  Triangle_iterator triangles_end(const Polyhedron& p) const
  {
    return faces(p).second;
  }

  Triangle_3 triangle(const Triangle_handle& thandle) const
  {
    typedef typename boost::graph_traits<Polyhedron>::face_descriptor::Graph Graph;
    typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;

    typename boost::graph_traits<Polyhedron>::face_descriptor handle = *thandle;
    const Graph& g = *(handle.graph);
    halfedge_descriptor hd = halfedge(handle.descriptor, g);
    typedef typename boost::property_map<Graph,vertex_point_t>::type Vpm;
    Vpm vpm = get(vertex_point, g);
    typedef typename K::Point_3 Point;
    const Point& a = get(vpm, target(hd,g));
    const Point& b = get(vpm, target(next(hd,g),g));
    const Point& c = get(vpm, target(next(next(hd,g),g),g));
    return Triangle_3(a,b,c);
  }
};


} // end namespace CGAL


#endif // POLYHEDRON_TRIANGLE_ACCESSOR_H
