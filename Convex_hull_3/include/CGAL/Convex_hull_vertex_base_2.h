// Copyright (c) 2011   Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé

// vertex of a triangulation of any dimension <= 3

#ifndef CGAL_CONVEX_HULL_VERTEX_BASE_2_H
#define CGAL_CONVEX_HULL_VERTEX_BASE_2_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/Triangulation_ds_vertex_base_2.h>
#include <CGAL/IO/io.h>

#include <iostream>

namespace CGAL {

template < typename GT,
           typename Vb = Triangulation_ds_vertex_base_2< > >
class Convex_hull_vertex_base_2
  : public Vb
{
public:
  typedef typename GT::Point_2                                       Point;

  typedef typename Vb::Face_handle                                  Face_handle;
  typedef typename Vb::Vertex_handle                                 Vertex_handle;

private:
  int _info = 0;
  Point _p;

public:
  template < typename TDS2 >
  struct Rebind_TDS
  {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other            Vb2;
    typedef Convex_hull_vertex_base_2<GT, Vb2>                       Other;
  };

  Convex_hull_vertex_base_2()
    : Vb() {}

  Convex_hull_vertex_base_2(const Point& p)
    : Vb(), _p(p) {}

  Convex_hull_vertex_base_2(const Point& p, Face_handle f)
    : Vb(f), _p(p) {}

  Convex_hull_vertex_base_2(Face_handle f)
    : Vb(f) {}

  void set_point(const Point& p) { _p = p; }
  const Point&  point() const { return _p; }
  Point& point() { return _p; }

  const int& info() const { return _info; }
  int&       info()       { return _info; }
};

template <typename GT, typename Vb>
std::istream&
operator>>(std::istream &is, Convex_hull_vertex_base_2<GT, Vb>& v)
{
  return is >> static_cast<Vb&>(v) >> v.point();
}

template <typename GT, typename Vb>
std::ostream&
operator<<(std::ostream &os, const Convex_hull_vertex_base_2<GT, Vb>& v)
{
  return os << static_cast<const Vb&>(v) << IO::serialize(v.point());
}

} //namespace CGAL

#endif // CGAL_CONVEX_HULL_VERTEX_BASE_2_H
