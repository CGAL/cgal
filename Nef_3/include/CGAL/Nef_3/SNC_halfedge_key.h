// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SNC_HALFEDGE_KEY_H
#define CGAL_SNC_HALFEDGE_KEY_H

#include <CGAL/license/Nef_3.h>

#include <CGAL/Kernel/global_functions.h>

namespace CGAL {

struct int_lt {
  bool operator()(const int& i1, const int& i2) const { return i1<i2; }
};

template <typename Edge_handle>
struct Halfedge_key_lt4 {

  bool operator()(const Edge_handle& e1, const Edge_handle& e2) const {
    if(CGAL::sign(e1->point().x()) != 0) {
      if(e1->source() != e2->source())
        return CGAL::compare_x(e1->source()->point(), e2->source()->point()) < 0;
      else
        return e1->point().x() < 0;
    }
    if(CGAL::sign(e1->point().y()) != 0) {
      if(e1->source() != e2->source())
        return CGAL::compare_y(e1->source()->point(), e2->source()->point()) < 0;
      else
        return e1->point().y() < 0;
    }
    if(e1->source() != e2->source())
      return CGAL::compare_z(e1->source()->point(), e2->source()->point()) < 0;
    return e1->point().z() < 0;
  }
};

template <typename Edge_handle>
struct Halfedge_key_lt3 {

  bool operator()(const Edge_handle& e1, const Edge_handle& e2) const {
    if(e1->source() != e2->source())
      return CGAL::lexicographically_xyz_smaller(e1->source()->point(), e2->source()->point());
    if(CGAL::sign(e1->point().x()) != 0)
      return e1->point().x() < 0;
    if(CGAL::sign(e1->point().y()) != 0)
      return e1->point().y() < 0;
    return e1->point().z() < 0;
  }
};

template <typename Point, typename Edge>
struct Halfedge_key {
  typedef Halfedge_key<Point,Edge> Self;
  Point p; int i; Edge e;
  Halfedge_key(Point pi, int ii, Edge ei) :
    p(pi), i(ii), e(ei) {}
  Halfedge_key(const Self& k) : p(k.p), i(k.i), e(k.e) {}
  Self& operator=(const Self& k) { p=k.p; i=k.i; e=k.e; return *this; }
  bool operator==(const Self& k) const { return p==k.p && i==k.i; }
  bool operator!=(const Self& k) const { return !operator==(k); }
};

template <typename Point, typename Edge, class Decorator>
struct Halfedge_key_lt {
  typedef Halfedge_key<Point,Edge> Key;
  typedef typename Point::R R;
  typedef typename R::Vector_3 Vector;
  typedef typename R::Direction_3 Direction;
  bool operator()( const Key& k1, const Key& k2) const {
    if( k1.e->source() == k2.e->source())
      return (k1.i < k2.i);
    Direction l(k1.e->vector());
    if( k1.i < 0) l = -l;
    return (Direction( k2.p - k1.p) == l);
  }
};

template <typename Point, typename Edge>
std::ostream& operator<<(std::ostream& os,
                         const Halfedge_key<Point,Edge>& k )
{ os << k.p << " " << k.i; return os; }

}
#endif //CGAL_SNC_HALFEDGE_KEY_H
