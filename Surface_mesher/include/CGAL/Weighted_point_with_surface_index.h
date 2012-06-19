// Copyright (c) 2005-2007  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_WEIGHTED_POINT_WITH_SURFACE_INDEX_H
#define CGAL_WEIGHTED_POINT_WITH_SURFACE_INDEX_H

#include <CGAL/Point_traits.h>
#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel_traits.h>

#include <string>

namespace CGAL {

template <class Weighted_point>
class Weighted_point_with_surface_index : public Weighted_point
{
  typedef CGAL::Point_traits<Weighted_point> Point_traits;
  typedef typename Point_traits::Bare_point Bare_point;
  typedef typename Kernel_traits<Bare_point>::Kernel::FT FT;

  CGAL_static_assertion((Is_weighted<Weighted_point>::value));
  CGAL_static_assertion((::boost::is_same<typename Point_traits::Is_weighted,
                                        Tag_true>::value));
  CGAL_static_assertion((::boost::is_same<Bare_point,
                                    typename Weighted_point::Point>::value));

public:
  Weighted_point_with_surface_index() : Weighted_point(), index(0) {}

  Weighted_point_with_surface_index(const Weighted_point& p)
    : Weighted_point(p), index(0) {}

  Weighted_point_with_surface_index(const Bare_point& bp)
    : Weighted_point(bp), index(0) {}

  Weighted_point_with_surface_index(const FT& x,
				    const FT& y,
				    const FT& z,
				    const FT& w)
    : Weighted_point(Bare_point(x, y, z, w)), index(0) {}

  Weighted_point_with_surface_index(const Bare_point& bp, 
                                    typename Weighted_point::Weight weight,
                                    int i) 
    : Weighted_point(bp,weight), index(i) {}

  Weighted_point_with_surface_index(const Weighted_point_with_surface_index& 
                                                                            pi)
    : Weighted_point(pi), index(pi.surface_index()) {}

  template <typename RT>
  Weighted_point_with_surface_index(const RT& x, const RT& y, const RT& z)
    : Weighted_point(Bare_point(x, y, z)), index(0) {}

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
    return Get_io_signature<Weighted_point>()() + "+i";
  }
#endif

private:
  int index;
}; // end class Weighted_point_with_surface_index

template <class P>
struct Is_weighted< Weighted_point_with_surface_index<P> > 
  : public Is_weighted<P> {};

template <class Point>
std::ostream&
operator<<(std::ostream &os, const Weighted_point_with_surface_index<Point>& p)
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
operator>>(std::istream &is, Weighted_point_with_surface_index<Point>& p)
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
