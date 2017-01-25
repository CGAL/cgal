// Copyright (c) 1998, 2001, 2003  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Olivier Devillers <Olivier.Devillers@sophia.inria.fr>
//                 Sylvain Pion

#ifndef CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_3_H
#define CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_3_H

#include <CGAL/license/Triangulation_3.h>


#include <CGAL/basic.h>

namespace CGAL {

template < class Vbb >
class Triangulation_hierarchy_vertex_base_3
  : public Vbb
{
  typedef Vbb                                           Base;
  typedef typename Base::Triangulation_data_structure   Tds;
public:
  typedef typename Tds::Vertex_handle                   Vertex_handle;
  typedef typename Tds::Cell_handle                     Cell_handle;
  typedef typename Base::Point                          Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<TDS2>::Other      Vb2;
    typedef Triangulation_hierarchy_vertex_base_3<Vb2>          Other;
  };

  Triangulation_hierarchy_vertex_base_3()
    : Base(), _up(), _down() {}

  Triangulation_hierarchy_vertex_base_3(const Point & p, Cell_handle f)
    : Base(p,f), _up(), _down() {}

  Triangulation_hierarchy_vertex_base_3(const Point & p)
    : Base(p), _up(), _down() {}

  Vertex_handle up()   const { return _up; }
  Vertex_handle down() const { return _down; }
  void set_up(Vertex_handle u)   { _up=u; }
  void set_down(Vertex_handle d) { _down=d; }

private:
  Vertex_handle _up;    // same vertex one level above
  Vertex_handle _down;  // same vertex one level below
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_3_H
