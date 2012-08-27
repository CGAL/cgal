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

#ifndef CGAL_KINETIC_CARTESIAN_MOVING_POINT_3_H_
#define CGAL_KINETIC_CARTESIAN_MOVING_POINT_3_H_
#include <CGAL/Kinetic/basic.h>
#include <iostream>

namespace CGAL { namespace Kinetic { namespace internal {

template <class Coordinate_t>
class Cartesian_moving_point_3
{
protected:
  typedef Cartesian_moving_point_3<Coordinate_t> This;

  //! What should I do for this
  typedef typename Coordinate_t::NT NT;
public:
  //typedef Static_point_t Static_point;

  //! The cartesian coordinate type
  typedef Coordinate_t Coordinate;

  //! initialize it from polys
  Cartesian_moving_point_3(const Coordinate &x, const Coordinate &y,
			   const Coordinate &z) {
    _coords[0]=x;
    _coords[1]=y;
    _coords[2]=z;
  }

  //! initialize it from a still point
  template <class Static_point>
  explicit Cartesian_moving_point_3(const Static_point &pt) {
    _coords[0]=Coordinate(typename Coordinate::NT(pt.x()));
    _coords[1]=Coordinate(typename Coordinate::NT(pt.y()));
    _coords[2]=Coordinate(typename Coordinate::NT(pt.z()));
  }

  //! null
  Cartesian_moving_point_3(){}

  //! homogeneous x
  const Coordinate &hx() const
  {
    return _coords[0];
  }

  //! homogeneous y
  const Coordinate &hy() const
  {
    return _coords[1];
  }

  //! homogeneous z
  const Coordinate &hz() const
  {
    return _coords[2];
  }

  //! homogeneous w
  const Coordinate hw() const
  {
    return Coordinate(1);
  }

  //! x
  const Coordinate &x() const
  {
    return _coords[0];
  }

  //! y
  const Coordinate &y() const
  {
    return _coords[1];
  }

  //! z
  const Coordinate &z() const
  {
    return _coords[2];
  }

  bool is_constant() const
  {
    for (unsigned int i=0; i< 3; ++i) {
      if (_coords[i].degree()>0) return false;
    }
    return true;
  }
  bool operator==(const This &o) const
  {
    return x()==o.x() && y()==o.y() && z()==o.z();
  }
#if 0
  //! Returns the value at time t.
  /*!
   */
  Static_point operator()(const NT &t) const
  {
    return Static_point(hx(t), hy(t), hz(t));
  }

  //! Non-operator version of operator()
  Static_point value_at( NT time) {
    return operator()(time);
  }
#endif


  //! Reverse the motion
  template <class NV>
  This transformed_coordinates(const NV &nv) const
  {
    return This(nv(_coords[0]), nv(_coords[1]), nv(_coords[2]));
  }

  template <class SK>
  struct Static_traits
  {
    typedef typename SK::Point_3 Static_type;
    static Static_type to_static(const This &o, const typename SK::FT &t, const SK &) {
      return Static_type(o.x()(t), o.y()(t), o.z()(t));
    }
  };
  template <class Converter>
  struct Coordinate_converter
  {
    Coordinate_converter(const Converter &c): c_(c){}
    typedef Cartesian_moving_point_3<typename Converter::argument_type> argument_type;
    typedef Cartesian_moving_point_3<typename Converter::result_type> result_type;

    result_type operator()(const argument_type &i) const
    {
      return result_type(c_(i.x()), c_(i.y()), c_(i.z()));
    }

    Converter c_;
  };
protected:
  Coordinate _coords[3];
};

template <class Coordinate>
std::ostream &operator<<(std::ostream &out, const Cartesian_moving_point_3<Coordinate> &point)
{
  out << point.x() << ", " << point.y() << ", " << point.z();
  return out;
}


template <class Coordinate>
std::istream &operator>>(std::istream &in,
			 Cartesian_moving_point_3<Coordinate> &point)
{
  Coordinate x, y, z;
  in >> x;
  char c;
  do {
    in >> c;
  } 
  while (std::isspace(c,std::locale::classic() ));

  if (c != ',') {
    in.setstate(std::ios_base::failbit);
    return in;
  }
  in >> y;
  do {
    in >> c;
  } while (std::isspace(c, std::locale::classic()));
  if (c != ',') {
    in.setstate(std::ios_base::failbit);
    return in;
  }
  in >> z;
  point= Cartesian_moving_point_3<Coordinate>(x,y,z);
  return in;
}


} } } //namespace CGAL::Kinetic::internal

/*namespace CGAL { namespace Kinetic {;

template <>
template <class Coord, class SK>
class To_static<typename internal::Cartesian_moving_point_3<Coord>, SK>:
public To_static_base<typename Coord::NT> {
public:
To_static(){}
To
typedef typename internal::Cartesian_moving_point_3<Coord> argument_type;
typedef typename SK::Point_3 result_type;
result_type operator()(const argument_type &arg) const {
return result_type(arg.x()(time()),
arg.y()(time()),
arg.z()(time()));
}
};
} } //namespace CGAL::Kinetic*/
#endif
