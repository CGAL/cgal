// Copyright (c) 1999-2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec
//                 Sylvain Pion

#ifndef CGAL_WEIGHTED_POINT_H
#define CGAL_WEIGHTED_POINT_H

#include <CGAL/license/Triangulation_2.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Weighted_point.h>"
#define CGAL_DEPRECATED_MESSAGE_DETAILS \
   "Weighted points are now part of the concept Kernel. One should therefore "\
   "use `CGAL::Weighted_point_2<K>` and `CGAL::Weighted_point_3<K>`."
#include <CGAL/internal/deprecation_warning.h>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Dimension.h>
#include <CGAL/IO/io.h>
#include <iostream>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/logical.hpp>

namespace CGAL {

template < class Pt, class We >
class Weighted_point : public Pt
{
  typedef typename Kernel_traits<Pt>::Kernel::FT FT;
public:
  typedef We Weight;
  typedef Pt Point;

  Weighted_point ()
      : Point(), _weight(0) {}

  //explicit
  Weighted_point (const Point &p)
      : Point(p), _weight(0)
  {
    // CGAL_error_msg( "Warning : truncated weight !!!");
  }

  Weighted_point (const Point &p, const Weight &w)
      : Point(p), _weight(w) {}


  // Constructors from coordinates are also provided for convenience, except
  // that they are only from Cartesian coordinates, and with no weight, so as
  // to avoid any potential ambiguity between the homogeneous weight and the
  // power weight (it should be easy enough to pass a Point explicitly in those
  // cases).
  // The enable_if complexity comes from the fact that we separate dimension 2 and 3.

  template < typename Tx, typename Ty >
  Weighted_point (const Tx &x, const Ty &y,
	          typename boost::enable_if< boost::mpl::and_<boost::is_convertible<Tx, FT>,
					                      boost::is_convertible<Ty, FT>,
							      boost::mpl::bool_<CGAL::Ambient_dimension<Point>::value == 2> > >::type* = 0)
      : Point(x, y), _weight(0) {}

  template < typename Tx, typename Ty, typename Tz >
  Weighted_point (const Tx &x, const Ty &y, const Tz &z,
	          typename boost::enable_if< boost::mpl::and_<boost::is_convertible<Tx, FT>,
					                      boost::is_convertible<Ty, FT>,
					                      boost::is_convertible<Tz, FT>,
							      boost::mpl::bool_<CGAL::Ambient_dimension<Point>::value == 3> > >::type* = 0)
      : Point(x, y, z), _weight(0) {}

  const Point & point() const
  {
      return *this;
  }

  const Weight & weight() const
  {
      return _weight;
  }

// The following power() member functions are not used at the moment.
// They belong to the traits class anyway.
//
//  Weight power(const Point &p)
//  {	
//      return squared_distance(*this, p) - weight();
//  }
// 
//  Weight power(const Weighted_point &p)
//  {	
//      return squared_distance(*this, p) - weight() - p.weight();
//  }

private:
  Weight _weight;
};


template < class Point, class Weight >
std::ostream &
operator<<(std::ostream &os, const Weighted_point<Point,Weight> &p)
{
  switch(get_mode(os))
  {
  case IO::ASCII :
    return os << p.point() <<  " " << p.weight();
  case IO::BINARY :
    os << p.point();
    write(os, p.weight());
    return os;
  default:
    return os << "Weighted_point(" << p.point() << ", " << p.weight() << ")";
  }
}

template < class Point, class Weight >
std::istream &
operator>>(std::istream &is, Weighted_point<Point,Weight> &wp)
{
  Weight w;
  Point p;
  is >> p;
  if(!is) return is;
  if(is_ascii(is))
    is >> w;
  else
    read(is, w);
  if (is)
    wp = Weighted_point<Point,Weight>(p,w);
  return is;
}

} // namespace CGAL

#endif // CGAL_WEIGHTED_POINT_H
