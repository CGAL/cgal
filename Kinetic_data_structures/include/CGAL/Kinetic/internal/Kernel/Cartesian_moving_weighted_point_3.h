// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_CARTESIAN_WEIGHTED_MOVING_POINT_3_H_
#define CGAL_KINETIC_CARTESIAN_WEIGHTED_MOVING_POINT_3_H_
#include <CGAL/Kinetic/basic.h>
#include <iostream>
#include <CGAL/Kinetic/internal/Kernel/Cartesian_moving_point_3.h>

namespace CGAL { namespace Kinetic { namespace internal {

template <class Coordinate_t>
class Cartesian_moving_weighted_point_3
{
protected:
  typedef Cartesian_moving_weighted_point_3<Coordinate_t> This;
  typedef Cartesian_moving_point_3<Coordinate_t> Point;
public:
  typedef Point Bare_point;
  typedef typename Point::Coordinate Coordinate;
  //! initialize it from polys
  Cartesian_moving_weighted_point_3(const Point &pt,
				    const Coordinate &w): point_(pt), weight_(w) {
  }

  Cartesian_moving_weighted_point_3(const Coordinate &x,
				    const Coordinate &y,
				    const Coordinate &z): point_(x,y,z), weight_(0) {
  }

  //! initialize it from a still point
  template <class Static_point>
  Cartesian_moving_weighted_point_3(const Static_point &pt): point_(pt.point()),
							     weight_(pt.weight()) {
  }

  bool operator==(const This & o) const {
    return weight() == o.weight() && point() == o.point();
  }

  //! null
  Cartesian_moving_weighted_point_3(){}

  const Coordinate &weight() const
  {
    return weight_;
  }

  const Point &point() const
  {
    return point_;
  }

  bool is_constant() const
  {
    if (weight_.degree() >0) return false;
    return point_.is_constant();
  }

  //! Reverse the motion
  template <class NV>
  This transformed_coordinates(const NV &nv) const
  {
    This ret(point_.transformed_coordinates(nv), nv(weight_));
    return ret;
  }

  template <class SK>
  struct Static_traits
  {
    typedef typename SK::Weighted_point_3 Static_type;
    static Static_type to_static(const This &o, const typename SK::FT &t, const SK&) {
      return Static_type(typename SK::Point_3(o.point().x()(t),
                                              o.point().y()(t),
                                              o.point().z()(t)),
                         o.weight()(t));
    }
  };

  template <class Converter>
  struct Coordinate_converter
  {
    Coordinate_converter(const Converter &c): c_(c), pc_(c){}
    typedef Cartesian_moving_weighted_point_3<typename Converter::argument_type> argument_type;
    typedef Cartesian_moving_weighted_point_3<typename Converter::result_type> result_type;

    result_type operator()(const argument_type &i) const
    {
      return result_type(pc_(i.point()), c_(i.weight()));
    }

    Converter c_;
    typename Bare_point::template Coordinate_converter<Converter> pc_;
  };

protected:
  Point point_;
  Coordinate weight_;
};

template <class Coordinate>
std::ostream &operator<<(std::ostream &out, const Cartesian_moving_weighted_point_3<Coordinate> &point)
{
  out << point.point() << ", " << point.weight();
  return out;
}


template <class Coordinate>
std::istream &operator>>(std::istream &in,
			 Cartesian_moving_weighted_point_3<Coordinate> &point)
{
  Coordinate w;
  typename Cartesian_moving_weighted_point_3<Coordinate>::Bare_point p;
  in >> p;
  char c;
  do {
    in >> c;
  }
   while (std::isspace(c,std::locale::classic() ));

  if (c != ',') {
    in.setstate(std::ios_base::failbit);
    return in;
  }
  in >> w;
  point= Cartesian_moving_weighted_point_3<Coordinate>(p,w);
  return in;
}


} } } //namespace CGAL::Kinetic::internal

/*namespace CGAL { namespace Kinetic {;

template <>
template <class Coord, class SK>
class To_static<typename internal::Cartesian_moving_weighted_point_3<Coord>, typename SK::Weighted_point>:
public To_static_base<typename Coord::NT,
typename internal::Cartesian_moving_weighted_point_3<Coord>,
typename SK::Weighted_point> {
typedef To_static_base<typename Coord::NT,
typename internal::Cartesian_moving_weighted_point_3<Coord>,
typename SK::Weighted_point>  P;
public:
To_static(){}
typename P::result_type operator()(const typename P::argument_type &arg) const {
return typename P::result_type(typename SK::Bare_point(arg.point().x()(P::time()),
arg.point().y()(P::time()),
arg.point().z()(P::time())),
arg.weight()(t_));
}
};
} } //namespace CGAL::Kinetic*/
#endif
