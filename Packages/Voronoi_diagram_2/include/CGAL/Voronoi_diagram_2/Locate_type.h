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

#ifndef CGAL_VORONOI_DIAGRAM_2_LOCATE_TYPE_H
#define CGAL_VORONOI_DIAGRAM_2_LOCATE_TYPE_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//====================================================================
//====================================================================

template<class V, class F, class E>
class Locate_type_base
{
 protected:
  enum Type { UNDEFINED = -1, VERTEX, FACE, EDGE };

  typedef V  Vertex_handle;
  typedef F  Face_handle;
  typedef E  Edge_rep;

  Locate_type_base() : t(UNDEFINED) {}
  Locate_type_base(const Vertex_handle& v) : t(VERTEX), v_(v) {}
  Locate_type_base(const Face_handle& f) : t(FACE), f_(f) {}
  Locate_type_base(const Edge_rep& e) : t(EDGE), e_(e) {}

 public:
  bool is_vertex() const { return t == VERTEX; }
  bool is_face() const { return t == FACE; }
  bool is_edge() const { return t == EDGE; }

  const Vertex_handle& vertex() const {
    CGAL_precondition( is_vertex() );
    return v_;
  }

  const Face_handle& face() const {
    CGAL_precondition( is_face() );
    return f_;
  }

  const Edge_rep& edge() const {
    CGAL_precondition( is_edge() );
    return e_;
  }

 protected:
  Type t;
  Vertex_handle v_;
  Face_handle   f_;
  Edge_rep      e_;
};

//====================================================================
//====================================================================


template<class C, bool is_vd> class Locate_type;
template<class C, bool is_vd> class Locate_type_accessor;

template<class VDA>
class Locate_type<VDA,true>
  : public Locate_type_base<typename VDA::Vertex_handle,
			    typename VDA::Face_handle,
			    typename VDA::Halfedge_handle>
{
  friend class Locate_type_accessor<VDA,true>;

 public:
  typedef typename VDA::Vertex_handle    Vertex_handle;
  typedef typename VDA::Face_handle      Face_handle;
  typedef typename VDA::Halfedge_handle  Halfedge_handle;

 protected:
  typedef Locate_type_base<Vertex_handle,Face_handle,Halfedge_handle> Base;
  typedef Locate_type<VDA,true> Self;

  Locate_type(const Vertex_handle& v) : Base(v) {}
  Locate_type(const Face_handle& f) : Base(f) {}
  Locate_type(const Halfedge_handle& e) : Base(e) {}

 public:
  Locate_type() : Base() {}
};


template<class DG>
class Locate_type<DG,false>
  : public Locate_type_base<typename DG::Vertex_handle,
			    typename DG::Face_handle,
			    typename DG::Edge>
{
  friend class Locate_type_accessor<DG,false>;

 public:
  typedef DG                          Dual_graph;
  typedef typename DG::Vertex_handle  Vertex_handle;
  typedef typename DG::Face_handle    Face_handle;
  typedef typename DG::Edge           Edge;

 protected:
  typedef Locate_type_base<Vertex_handle,Face_handle,Edge> Base;
  typedef Locate_type<DG,false> Self;

  Locate_type(const Vertex_handle& v) : Base(v) {}
  Locate_type(const Face_handle& f) : Base(f) {}
  Locate_type(const Edge& e) : Base(e) {}

 public:
  Locate_type() : Base() {}
};


//====================================================================
//====================================================================

template<class C, bool b>
struct Locate_type_accessor
{
  typedef Locate_type<C,b> Locate_type;

  template<class T>
  Locate_type make_locate_type(const T& t) const {
    return Locate_type(t);
  }
};

//====================================================================
//====================================================================

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_LOCATE_TYPE_H
