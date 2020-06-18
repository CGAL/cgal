// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mariette Yvinec,Sylvain Pion

// face of a triangulation of any dimension <=3

#ifndef CGAL_TRIANGULATION_FACE_BASE_WITH_INFO_2_H
#define CGAL_TRIANGULATION_FACE_BASE_WITH_INFO_2_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/Triangulation_face_base_2.h>

namespace CGAL {

template < typename Info_, typename GT,
           typename Fb_ = Triangulation_face_base_2<GT> >
class Triangulation_face_base_with_info_2
  : public Fb_
{
  Info_ _info;
public:
  typedef typename Fb_::Vertex_handle                   Vertex_handle;
  typedef typename Fb_::Face_handle                     Face_handle;
  typedef Info_                                        Info;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb_::template Rebind_TDS<TDS2>::Other       Fb2;
    typedef Triangulation_face_base_with_info_2<Info, GT, Fb2>  Other;
  };

  Triangulation_face_base_with_info_2()
    : Fb_() {}

  Triangulation_face_base_with_info_2(Vertex_handle v0,
                                      Vertex_handle v1,
                                      Vertex_handle v2)
    : Fb_(v0, v1, v2) {}

  Triangulation_face_base_with_info_2(Vertex_handle v0,
                                      Vertex_handle v1,
                                      Vertex_handle v2,
                                      Face_handle   n0,
                                      Face_handle   n1,
                                      Face_handle   n2 )
    : Fb_(v0, v1, v2, n0, n1, n2) {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_FACE_BASE_WITH_INFO_2_H
