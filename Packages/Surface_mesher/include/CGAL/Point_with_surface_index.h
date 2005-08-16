// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POINT_WITH_SURFACE_INDEX_H
#define CGAL_POINT_WITH_SURFACE_INDEX_H

#include <CGAL/Point_traits.h>

namespace CGAL {

template <class Point>
class Point_with_surface_index : public Point
{
  typedef Point_traits<Point> Point_traits;
  typedef typename Point_traits::Bare_point Bare_point;

public:
  Point_with_surface_index() : Point(), index(0) {}

  Point_with_surface_index(const Point& p) : Point(p), index(0) {}

//   Point_with_surface_index(const Bare_point& bp) : Point(bp), index(0) {}

//   Point_with_surface_index(const Bare_point& bp, 
// 			   typename Point::Weight weight,
// 			   int i) 
//     : Point(bp,weight), index(i) {}
			   

  Point_with_surface_index(const Point_with_surface_index& pi)
    : Point(pi), index(pi.surface_index()) {}

  

//   Point_with_surface_index operator=(const Point& p)
//   {
//     static_cast<Point*>(this)->operator=(p);
//     index = 0;
//   }

//   Point_with_surface_index operator=(const Bare_point& bp)
//   {
//     static_cast<Point*>(this)->operator=(bp);
//     index = 0;
//   }

  template <typename RT>
  Point_with_surface_index(const RT& x, const RT& y, const RT& z)
    : Point(Point_traits().point(Bare_point(x, y, z))), index(0) {}

  int surface_index() const
  {
    return index;
  }

  void set_surface_index(const int i)
  {
    index = i;
  }
private:
  int index;
}; // end class Point_with_surface_index

template <class GT>
class Point_with_surface_index_geom_traits : public GT
{
  typedef typename GT::Point_3 Old_Point_3;

public:
  typedef Point_with_surface_index<Old_Point_3> Point_3;
  typedef Point_with_surface_index<Old_Point_3> Weighted_point_3;

};  // end Point_with_surface_index_geom_traits

} // end namespace CGAL

#endif
