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

#ifndef CGAL_VORONOI_DIAGRAM_2_LOCATE_RESULT_H
#define CGAL_VORONOI_DIAGRAM_2_LOCATE_RESULT_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//====================================================================
//====================================================================

template<class V, class F, class E>
class Locate_result_base
{
 protected:
  enum Type { UNDEFINED = -1, VERTEX, FACE, EDGE };

  typedef V  Vertex_handle;
  typedef F  Face_handle;
  typedef E  Edge_rep;

  typedef Locate_result_base<Vertex_handle,Face_handle,Edge_rep>  Self;

  Locate_result_base() : t(UNDEFINED) {}

  void set(const Vertex_handle& v) { t = VERTEX;  v_ = v; }
  void set(const Face_handle& f)   { t = FACE;    f_ = f; }
  void set(const Edge_rep& e)      { t = EDGE;    e_ = e; }

 public:
  bool is_vertex() const { return t == VERTEX; }
  bool is_face() const { return t == FACE; }
  bool is_edge() const { return t == EDGE; }

  operator const Vertex_handle&() const {
    CGAL_precondition( is_vertex() );
    return v_;
  }

  operator const Face_handle&() const {
    CGAL_precondition( is_face() );
    return f_;
  }

  operator const Edge_rep&() const {
    CGAL_precondition( is_edge() );
    return e_;
  }

  bool operator==(const Self& o) const {
    if ( t != o.t ) { return false; }
    switch (t) {
    case VERTEX:
      return v_ == o.v_;
    case EDGE:
      return f_ == o.f_;
    case FACE:
      return e_ == o.e_;
    default:
      CGAL_assertion( t == UNDEFINED && o.t == UNDEFINED );
      return true;
    }
  }

  bool operator!=(const Self& o) const {
    return !((*this) == o);
  }

 protected:
  Type t;
  Vertex_handle v_;
  Face_handle   f_;
  Edge_rep      e_;
};

//====================================================================
//====================================================================


template<class C, bool is_vd> class Locate_result;
template<class C, bool is_vd> class Locate_result_accessor;

template<class VDA>
class Locate_result<VDA,true>
  : public Locate_result_base<typename VDA::Vertex_handle,
			      typename VDA::Face_handle,
			      typename VDA::Halfedge_handle>
{
  friend class Locate_result_accessor<VDA,true>;

 public:
  typedef typename VDA::Vertex_handle    Vertex_handle;
  typedef typename VDA::Face_handle      Face_handle;
  typedef typename VDA::Halfedge_handle  Halfedge_handle;

 protected:
  typedef Locate_result_base<Vertex_handle,Face_handle,Halfedge_handle> Base;
  typedef Locate_result<VDA,true> Self;
};


template<class DG>
class Locate_result<DG,false>
  : public Locate_result_base<typename DG::Vertex_handle,
			      typename DG::Face_handle,
			      typename DG::Edge>
{
  friend class Locate_result_accessor<DG,false>;

 public:
  typedef DG                          Delaunay_graph;
  typedef typename DG::Vertex_handle  Vertex_handle;
  typedef typename DG::Face_handle    Face_handle;
  typedef typename DG::Edge           Edge;

 protected:
  typedef Locate_result_base<Vertex_handle,Face_handle,Edge> Base;
  typedef Locate_result<DG,false> Self;
};


//====================================================================
//====================================================================

template<class C, bool b>
struct Locate_result_accessor
{
  typedef Locate_result<C,b> Locate_result;

  template<class T>
  static Locate_result make_locate_result(const T& t) {
    Locate_result lr;
    lr.set(t);
    return lr;
  }
};

//====================================================================
//====================================================================

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_LOCATE_RESULT_H
