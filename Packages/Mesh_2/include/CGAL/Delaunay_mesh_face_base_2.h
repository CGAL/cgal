// Copyright (c) 2001-2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DELAUNAY_FACE_BASE_2_H
#define CGAL_DELAUNAY_FACE_BASE_2_H

#include <CGAL/Constrained_triangulation_face_base_2.h>

namespace CGAL {

template <class Gt,
          class Fb = Constrained_triangulation_face_base_2<Gt> >
class Delaunay_mesh_face_base_2 : public Fb
{
public:
  typedef Gt Geom_traits;
  typedef Constrained_triangulation_face_base_2<Gt> CTFb;
  typedef typename Fb::Vertex_handle Vertex_handle;
  typedef typename Fb::Face_handle Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename CTFb::template Rebind_TDS<TDS2>::Other Fb2;
    typedef Delaunay_mesh_face_base_2<Gt,Fb2> Other;
  };

protected:
  bool marked;

public:
  Delaunay_mesh_face_base_2(): Fb(), marked(false) {};

  Delaunay_mesh_face_base_2(Vertex_handle v0, 
			    Vertex_handle v1, 
			    Vertex_handle v2)
    : Fb(v0,v1,v2), marked(false) {};

  Delaunay_mesh_face_base_2(Vertex_handle v0, 
			    Vertex_handle v1, 
			    Vertex_handle v2,
			    Face_handle n0, 
			    Face_handle n1, 
			    Face_handle n2)
    : Fb(v0,v1,v2,n0,n1,n2), marked(false) {};

  inline
  bool is_marked() const { return marked; };

  inline
  void set_marked(const bool b) { marked=b; };
};

}; // namespace CGAL

#endif
