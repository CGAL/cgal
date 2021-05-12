// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POINT_WITH_SURFACE_INDEX_H
#define CGAL_POINT_WITH_SURFACE_INDEX_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Point_traits.h>

#include <string>

namespace CGAL {

template <class Point>
class Point_with_surface_index : public Point
{
  typedef CGAL::Point_traits<Point> Point_traits;
  typedef typename Point_traits::Bare_point Bare_point;
  typedef typename Kernel_traits<Bare_point>::Kernel Kernel;
  typedef typename Kernel::FT FT;

public:
  Point_with_surface_index() : Point(), index(0) {}

  Point_with_surface_index(const Point& p) : Point(p), index(0) {}

  Point_with_surface_index(const FT& x, const FT& y, const FT& z, const FT& w = FT(1))
    : Point(Point_traits().point(Bare_point(x, y, z, w))), index(0) {}

  int surface_index() const
  {
    return index;
  }

  void set_surface_index(const int i)
  {
    index = i;
  }

#ifdef CGAL_MESH_3_IO_H
  static
  std::string io_signature()
  {
    return Get_io_signature<Point>()() + "+i";
  }
#endif
private:
  int index;
}; // end class Point_with_surface_index

template <class Point>
std::ostream&
operator<<(std::ostream &os, const Point_with_surface_index<Point>& p)
{
  os << static_cast<const Point&>(p);
  if(is_ascii(os))
    os << ' ' << p.surface_index();
  else
    write(os, p.surface_index());
  return os;
}

template <class Point>
std::istream&
operator>>(std::istream &is, Point_with_surface_index<Point>& p)
{
  is >>  static_cast<Point&>(p);
  int index;
  if(is_ascii(is))
    is >> index;
  else
    read(is, index);
  p.set_surface_index(index);
  return is;
}

} // end namespace CGAL

#endif
