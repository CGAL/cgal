// ============================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Triangulation_face_base_with_info_2.h
// revision      : $Revision$
// author(s)     : Mariette Yvinec,Sylvain Pion
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

// face of a triangulation of any dimension <=3

#ifndef CGAL_TRIANGULATION_FACE_BASE_WITH_INFO_2_H
#define CGAL_TRIANGULATION_FACE_BASE_WITH_INFO_2_H

#include <CGAL/Triangulation_face_base_2.h>

CGAL_BEGIN_NAMESPACE

template < typename Info_, typename GT,
           typename Fb = Triangulation_face_base_2<GT> >
class Triangulation_face_base_with_info_2
  : public Fb
{
  Info_ _info;
public:
  typedef typename Fb::Vertex_handle                   Vertex_handle;
  typedef typename Fb::Face_handle                     Face_handle;
  typedef Info_                                        Info;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other       Fb2;
    typedef Triangulation_face_base_with_info_2<Info, GT, Fb2>  Other;
  };

  Triangulation_face_base_with_info_2()
    : Fb() {}

  Triangulation_face_base_with_info_2(Vertex_handle v0, 
				      Vertex_handle v1,
                                      Vertex_handle v2)
    : Fb(v0, v1, v2) {}

  Triangulation_face_base_with_info_2(Vertex_handle v0, 
				      Vertex_handle v1,
                                      Vertex_handle v2, 
                                      Face_handle   n0, 
				      Face_handle   n1,
                                      Face_handle   n2 )
    : Fb(v0, v1, v2, v3, n0, n1, n2, n3) {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_FACE_BASE_WITH_INFO_2_H
