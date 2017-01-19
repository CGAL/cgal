// Copyright (c) 2011   Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Andreas Fabri

// face of a triangulation of any dimension <=3

#ifndef CGAL_CONVEX_HULL_FACE_BASE_2_H
#define CGAL_CONVEX_HULL_FACE_BASE_2_H

#include <CGAL/license/Convex_hull_3.h>


#include <list>
#include <CGAL/Triangulation_face_base_2.h>

namespace CGAL {

template < typename Info_, typename GT,
           typename Fb = Triangulation_face_base_2<GT> >
class Convex_hull_face_base_2
  : public Fb
{
  Info_ _info;
public:
  typedef typename Fb::Vertex_handle                   Vertex_handle;
  typedef typename Fb::Face_handle                     Face_handle;
  typedef Info_                                        Info;

  typename std::list<Face_handle>::iterator it;
  std::list<typename GT::Point_3> points;
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other       Fb2;
    typedef Convex_hull_face_base_2<Info, GT, Fb2>  Other;
  };

  Convex_hull_face_base_2()
    : Fb(), _info(0) {}

  Convex_hull_face_base_2(Vertex_handle v0, 
				      Vertex_handle v1,
                                      Vertex_handle v2)
    : Fb(v0, v1, v2), _info(0) {}

  Convex_hull_face_base_2(Vertex_handle v0, 
				      Vertex_handle v1,
                                      Vertex_handle v2, 
                                      Face_handle   n0, 
				      Face_handle   n1,
                                      Face_handle   n2 )
    : Fb(v0, v1, v2, n0, n1, n2), _info(0) {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

} //namespace CGAL

#endif // CGAL_CONVEX_HULL_FACE_BASE_2_H
