// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Olivier Devillers <Olivivier.Devillers@sophia.inria.fr>
//                 Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_2_H
#define CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_2_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/basic.h>

namespace CGAL {

template < class Vbb>
class Triangulation_hierarchy_vertex_base_2
 : public Vbb
{
  typedef Vbb                                           Base;
  typedef typename Base::Triangulation_data_structure   Tds;

public:
  typedef typename Base::Point            Point;
  typedef Tds                             Triangulation_data_structure;
  typedef typename Tds::Vertex_handle     Vertex_handle;
  typedef typename Tds::Face_handle       Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<TDS2>::Other      Vb2;
    typedef Triangulation_hierarchy_vertex_base_2<Vb2>          Other;
  };

  Triangulation_hierarchy_vertex_base_2()
    : Base()
    {}
  Triangulation_hierarchy_vertex_base_2(const Point & p, Face_handle f)
    : Base(p,f)
    {}
  Triangulation_hierarchy_vertex_base_2(const Point & p)
    : Base(p)
    {}

  Vertex_handle up() {return _up;}
  Vertex_handle down() {return _down;}
  void set_up(Vertex_handle u) {_up=u;}
  void set_down(Vertex_handle d) { _down=d;}


 private:
  Vertex_handle  _up   = nullptr;    // same vertex one level above
  Vertex_handle  _down = nullptr;  // same vertex one level below
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_HIERARCHY_VERTEX_BASE_2_H
