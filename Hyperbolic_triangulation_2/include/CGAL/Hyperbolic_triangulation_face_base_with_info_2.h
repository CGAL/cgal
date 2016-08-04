// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud

// (based on Triangulation_face_base_with_info_2.h)

#ifndef CGAL_HYPERBOLIC_TRIANGULATION_FACE_BASE_WITH_INFO_2_H
#define CGAL_HYPERBOLIC_TRIANGULATION_FACE_BASE_WITH_INFO_2_H

#include <CGAL/Hyperbolic_triangulation_face_base_2.h>

namespace CGAL {

template < typename Info_, typename GT,
           typename Fb_ = Hyperbolic_triangulation_face_base_2<GT> >
class Hyperbolic_triangulation_face_base_with_info_2
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
    typedef Hyperbolic_triangulation_face_base_with_info_2<Info, GT, Fb2>  Other;
  };

  Hyperbolic_triangulation_face_base_with_info_2()
    : Fb_()
   {}

  Hyperbolic_triangulation_face_base_with_info_2(Vertex_handle v0, 
				      Vertex_handle v1,
                                      Vertex_handle v2)
    : Fb_(v0, v1, v2)
   {}

  Hyperbolic_triangulation_face_base_with_info_2(Vertex_handle v0, 
				      Vertex_handle v1,
                                      Vertex_handle v2, 
                                      Face_handle   n0, 
				      Face_handle   n1,
                                      Face_handle   n2 )
    : Fb_(v0, v1, v2, n0, n1, n2)
   {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

} //namespace CGAL

#endif // CGAL_HYPERBOLIC_TRIANGULATION_FACE_BASE_WITH_INFO_2_H
