// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium

// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Triangulation_face_base_2.h
// package       : Triangulation_2 
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_TRIANGULATION_FACE_BASE_2_H
#define CGAL_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_ds_face_base_2.h>

CGAL_BEGIN_NAMESPACE 

template < typename Gt, typename Fb = Triangulation_ds_face_base_2<> >
class Triangulation_face_base_2 
  : public Fb
{
public:
  typedef typename Fb::Vertex_handle                   Vertex_handle;
  typedef typename Fb::Face_handle                     Face_handle;

  typedef GT                                           Geom_traits;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Triangulation_face_base_3<GT, Fb2>             Other;
  };

public:
  Triangulation_ds_face_base_2()
       : Fb() {}

  Triangulation_ds_face_base_2(Vertex_handle v0, 
			       Vertex_handle v1, 
			       Vertex_handle v2)
    : Fb(v0,v1,v2) {}

  Triangulation_ds_face_base_2(Vertex_handle v0, 
			       Vertex_handle v1, 
			       Vertex_handle v2,
			       Face_handle n0, 
			       Face_handle n1, 
			       Face_handle n2)
    : Fb(v0,v1,v2,n0,n1,n2) {}
};


CGAL_END_NAMESPACE 

#endif //CGAL_TRIANGULATION_FACE_BASE_2_H
