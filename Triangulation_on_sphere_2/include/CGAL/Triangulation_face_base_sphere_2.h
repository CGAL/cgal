// Copyright (c) 1997, 2012, 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the Licenxse, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mariette Yvinec, Claudia Werner

#ifndef CGAL_TRIANGULATION_FACE_BASE_SPHERE_2_H
#define CGAL_TRIANGULATION_FACE_BASE_SPHERE_2_H

#include <CGAL/Triangulation_ds_face_base_2.h>

namespace CGAL {

template < typename Gt, typename Fb = Triangulation_ds_face_base_2<> >
class Triangulation_face_base_sphere_2
  : public Fb
{
public:
  typedef Gt                                           Geom_traits;
  typedef typename Fb::Vertex_handle                   Vertex_handle;
  typedef typename Fb::Face_handle                     Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Triangulation_face_base_sphere_2<Gt, Fb2>      Other;
  };

  unsigned char _in_conflict_flag;

public:
  void set_in_conflict_flag(unsigned char f) { _in_conflict_flag = f; }
  unsigned char get_in_conflict_flag() const { return _in_conflict_flag; }

public:
  Triangulation_face_base_sphere_2()
    : Fb(), _ghost(false)
  {
    set_in_conflict_flag(0);
  }

  Triangulation_face_base_sphere_2(Vertex_handle v0,
                                   Vertex_handle v1,
                                   Vertex_handle v2)
    : Fb(v0, v1, v2), _ghost(false)
  {
    set_in_conflict_flag(0);
  }

  Triangulation_face_base_sphere_2(Vertex_handle v0,
                                   Vertex_handle v1,
                                   Vertex_handle v2,
                                   Face_handle n0,
                                   Face_handle n1,
                                   Face_handle n2)
    : Fb(v0, v1, v2, n0, n1, n2), _ghost(false)
  {
    set_in_conflict_flag(0);
  }

  const bool& is_ghost() const { return _ghost; }
  bool& ghost() { return _ghost; }

protected:
  bool _ghost;
};

} // namespace CGAL

#endif //CGAL_Triangulation_face_base_sphere_2_H
