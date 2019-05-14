// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_APOLLONIUS_SITE_2_H
#define CGAL_APOLLONIUS_SITE_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <iostream>
#include <CGAL/basic.h>

namespace CGAL {


template < class K >
class Apollonius_site_2
{
public:
  typedef K                       Kernel;
  typedef typename K::Point_2     Point_2;
  typedef Apollonius_site_2<K>    Self;
  typedef typename K::FT          Weight;
  typedef typename K::RT          RT;
  typedef Weight                  FT;


public:
  Apollonius_site_2(const Point_2& p = Point_2(),
		    const Weight& w = Weight(0))
    : _p(p), _w(w) {}

  Apollonius_site_2(const Apollonius_site_2& other)
    : _p(other._p), _w(other._w) {}

  
  const Point_2& point() const { return _p; }
  const Weight&  weight() const { return _w; }
  RT      x() const { return _p.x(); }
  RT      y() const { return _p.y(); }

  bool operator==(const Apollonius_site_2& other) const
  {
    return ( _p == other._p && _w == other._w );
  }

private:
  Point_2 _p;
  Weight  _w;

};

template <class K>
std::ostream&
operator<<(std::ostream &os, const Apollonius_site_2<K>& wp)
{
  return os << wp.point() << " " << wp.weight();
}

template <class K>
std::istream&
operator>>(std::istream& is, Apollonius_site_2<K>& wp)
{
  typename Apollonius_site_2<K>::Weight   weight (0.);
  typename Apollonius_site_2<K>::Point_2  p (0., 0.);
  is >> p >> weight;
  if(is) {
    wp = Apollonius_site_2<K>(p, weight);
  }
  return is;
}


} //namespace CGAL

#endif // CGAL_APOLLONIUS_SITE_2_H
