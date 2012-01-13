// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>

#ifndef CGAL_ALPHA_SHAPE_FACE_BASE_2_H
#define CGAL_ALPHA_SHAPE_FACE_BASE_2_H

#include <CGAL/utility.h>

namespace CGAL {

template < class Gt, class Fb = Triangulation_face_base_2<Gt> >
class Alpha_shape_face_base_2 : public Fb
{
  typedef typename Fb::Triangulation_data_structure  TDS;
public:
  typedef TDS                                Triangulation_data_structure;
  typedef typename TDS::Vertex_handle        Vertex_handle;
  typedef typename TDS::Face_handle          Face_handle;

  typedef typename Gt::FT                    Coord_type;
  typedef Triple<Coord_type, Coord_type, Coord_type> Interval_3;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other    Fb2;
    typedef Alpha_shape_face_base_2<Gt,Fb2>         Other;
  };


private:
  Interval_3 vec_edge[3];
  Coord_type A;
 
public:
  Alpha_shape_face_base_2()  : Fb()     {}
  
  Alpha_shape_face_base_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2)
    : Fb(v0, v1, v2)     {}
  
  Alpha_shape_face_base_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
			  Face_handle n0, Face_handle n1, Face_handle n2)
    : Fb(v0, v1, v2, n0, n1, n2)
    {}

  const Coord_type & get_alpha() const
    {
      return A;
    }

  void set_alpha(const Coord_type & AA)
    {
      A = AA;
    }

  const Interval_3 & get_ranges(int i) const
    {
      return vec_edge[i];
    }

  void set_ranges(int i, const Interval_3& Inter)
    {
      vec_edge[i]=Inter;
    }


};

} //namespace CGAL

#endif //ALPHA_SHAPE_FACE_BASE_2_H
