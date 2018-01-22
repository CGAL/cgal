// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//                 Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//                 Nico Kruithof  <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_HIERARCHY_VERTEX_BASE_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_HIERARCHY_VERTEX_BASE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>


#include <CGAL/basic.h>
#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>

#include <iostream>

namespace CGAL
{

template < class Gt, class Vb = CGAL::Periodic_2_triangulation_vertex_base_2<Gt> >
class Periodic_2_triangulation_hierarchy_vertex_base_2
  : public Vb
{
  typedef Vb Base;
  typedef typename Vb::Triangulation_data_structure    Tds;
  typedef Gt                                    Geom_traits;
  typedef typename Gt::Point_2                  Point;
  typedef Tds                                   Triangulation_data_structure;
  typedef typename Tds::Face_handle             Face_handle;
  typedef typename Tds::Vertex_handle           Vertex_handle;

  template < typename Tds2 >
  struct Rebind_Tds
  {
    typedef typename Vb::template Rebind_Tds<Tds2>::Other  Vb2;
    typedef Triangulation_vertex_base_2<Gt, Vb2>           Other;
  };
  typedef typename Vb::Offset                                  Offset;

  Periodic_2_triangulation_hierarchy_vertex_base_2()
    : Base(), _up(0), _down(0)
  {}
  Periodic_2_triangulation_hierarchy_vertex_base_2(const Point & p, Face_handle f)
    : Base(p, f), _up(0), _down(0)
  {}
  Periodic_2_triangulation_hierarchy_vertex_base_2(const Point & p)
    : Base(p), _up(0), _down(0)
  {}

  Vertex_handle up()
  {
    return _up;
  }
  Vertex_handle down()
  {
    return _down;
  }
  void set_up(Vertex_handle u)
  {
    _up = u;
  }
  void set_down(Vertex_handle d)
  {
    _down = d;
  }


private:
  Vertex_handle  _up;    // same vertex one level above
  Vertex_handle  _down;  // same vertex one level below
};

template < class Tds >
inline
std::istream&
operator>>(std::istream &is, Periodic_2_triangulation_hierarchy_vertex_base_2<Tds> &)
// no combinatorial information.
{
  return is;
}

template < class Tds >
inline
std::ostream&
operator<<(std::ostream &os,
           const Periodic_2_triangulation_hierarchy_vertex_base_2<Tds> &)
// no combinatorial information.
{
  return os;
}

} //namespace CGAL

#endif // CGAL_PERIODIC_2_TRIANGULATION_HIERARCHY_VERTEX_BASE_2_H
