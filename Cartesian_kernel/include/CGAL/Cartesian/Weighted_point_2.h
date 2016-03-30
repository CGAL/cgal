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
// 
//
// Author(s)     : Mariette Yvinec
//                 Sylvain Pion

#ifndef CGAL_CARTESIAN_WEIGHTED_POINT_2_H
#define CGAL_CARTESIAN_WEIGHTED_POINT_2_H

#include <iostream>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Dimension.h>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/logical.hpp>

namespace CGAL {

template < class R_ >
class Weighted_pointC2 : public R_::Point_2
{
  typedef typename R_::FT FT;
public:
  typedef FT Weight;
  typedef typename R_::Point_2 Point;

  Weighted_pointC2 ()
      : Point(), _weight(0) {}

  //explicit
  Weighted_pointC2 (const Point &p)
      : Point(p), _weight(0)
  {
    // CGAL_error_msg( "Warning : truncated weight !!!");
  }

  Weighted_pointC2 (const Point &p, const Weight &w)
      : Point(p), _weight(w) {}


  // Constructors from coordinates are also provided for convenience, except
  // that they are only from Cartesian coordinates, and with no weight, so as
  // to avoid any potential ambiguity between the homogeneous weight and the
  // power weight (it should be easy enough to pass a Point explicitly in those
  // cases).
  // The enable_if complexity comes from the fact that we separate dimension 2 and 3.


  template < typename Tx, typename Ty>
  Weighted_pointC2 (const Tx &x, const Ty &y,
	          typename boost::enable_if< boost::mpl::and_<boost::is_convertible<Tx, FT>,
					                      boost::is_convertible<Ty, FT>,
							      boost::mpl::bool_<CGAL::Ambient_dimension<Point>::value == 2> > >::type* = 0)
      : Point(x, y), _weight(0) {}

  const Point & point() const
  {
      return *this;
  }

  const Weight & weight() const
  {
      return _weight;
  }


private:
  Weight _weight;
};


template < class R_ >
std::ostream &
operator<<(std::ostream &os, const Weighted_pointC2<R_> &p)
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

template < class R_ >
std::istream &
operator>>(std::istream &is, Weighted_pointC2<R_> &wp)
{
  typename Weighted_pointC2<R_>::Weight w;
  typename Weighted_pointC2<R_>::Point p;
  is >> p;
  if(!is) return is;
  if(is_ascii(is))
    is >> w;
  else
    read(is, w);
  if (is)
    wp = Weighted_point_2<R_>(p,w);
  return is;
}

} // namespace CGAL

#endif // CGAL_CARTESIAN_WEIGHTED_POINT_2_H
