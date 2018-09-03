// Copyright (c) 2016  INRIA Nancy (France).
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
// $URL: 
// $Id: 
// 
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@inria.fr>
//                 Mikhail Bogdanov

#ifndef CGAL_HYPERBOLIC_TRIANGULATION_FACE_BASE_2_H
#define CGAL_HYPERBOLIC_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/Object.h>

namespace CGAL { 

template < typename Gt, typename Fb = Triangulation_ds_face_base_2<> >
class Hyperbolic_triangulation_face_base_2 
  : public Fb
{
public:
  typedef Gt                                           Geom_traits;
  typedef typename Fb::Vertex_handle                   Vertex_handle;
  typedef typename Fb::Face_handle                     Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Hyperbolic_triangulation_face_base_2<Gt, Fb2>             Other;
  };

public:
  Hyperbolic_triangulation_face_base_2()
    : Fb()
  {}

  Hyperbolic_triangulation_face_base_2(Vertex_handle v0, 
			               Vertex_handle v1, 
			               Vertex_handle v2)
    : Fb(v0,v1,v2)
  {}

  Hyperbolic_triangulation_face_base_2(Vertex_handle v0, 
			               Vertex_handle v1, 
			               Vertex_handle v2,
			               Face_handle n0, 
			               Face_handle n1, 
			               Face_handle n2)
    : Fb(v0,v1,v2,n0,n1,n2)
  {}

  static int ccw(int i) {return Triangulation_cw_ccw_2::ccw(i);}
  static int  cw(int i) {return Triangulation_cw_ccw_2::cw(i);}

  CGAL::Object& tds_data() {
    return this->_tds_data; 
  }

  const CGAL::Object& tds_data() const {
    return this->_tds_data;
  }  


private:

  CGAL::Object 	_tds_data;
};

} //namespace CGAL 

#endif //CGAL_HYPERBOLIC_TRIANGULATION_FACE_BASE_2_H
