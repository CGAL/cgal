// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Stefan Schirra

#ifndef CGAL_DIRECTION_2_H
#define CGAL_DIRECTION_2_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R_>
class Direction_2 : public R_::Kernel_base::Direction_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Vector_2              Vector_2;
  typedef typename R_::Line_2                Line_2;
  typedef typename R_::Ray_2                 Ray_2;
  typedef typename R_::Segment_2             Segment_2;
  typedef typename R_::Aff_transformation_2  Aff_transformation_2;
  typedef typename R_::Kernel_base::Direction_2      RDirection_2;

  typedef Direction_2                        Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Direction_2>::value));

public:

  typedef Dimension_tag<2>  Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  typedef RDirection_2 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_   R;

  Direction_2() {}

  Direction_2(const RDirection_2& d)
    : RDirection_2(d) {}

  explicit Direction_2(const Vector_2& v)
    : RDirection_2(typename R::Construct_direction_2()(Return_base_tag(), v)) {}

  explicit Direction_2(const Line_2& l)
    : RDirection_2(typename R::Construct_direction_2()(Return_base_tag(), l)) {}

  explicit Direction_2(const Ray_2& r)
    : RDirection_2(typename R::Construct_direction_2()(Return_base_tag(), r)) {}

  explicit Direction_2(const Segment_2& s)
    : RDirection_2(typename R::Construct_direction_2()(Return_base_tag(), s)) {}

  Direction_2(const RT &x, const RT &y)
    :  RDirection_2(typename R::Construct_direction_2()(Return_base_tag(), x,y)) {}

  typename R::Boolean
  counterclockwise_in_between(const Direction_2 &d1,
			      const Direction_2 &d2) const
  {
    return R().counterclockwise_in_between_2_object()(*this, d1, d2);
  }

  Direction_2 perpendicular(const Orientation &o) const
  {
    return R().construct_perpendicular_direction_2_object()(*this,o);
  }

  typename cpp11::result_of<typename R::Compute_dx_2( Direction_2)>::type
  dx() const
  {
    return R().compute_dx_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_dy_2( Direction_2)>::type
  dy() const
  {
    return R().compute_dy_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_dx_2( Direction_2)>::type
  delta(int i) const
  {
    CGAL_kernel_precondition( ( i == 0 ) || ( i == 1 ) );
    return (i==0) ? dx() : dy();
  }

  typename R::Boolean
  operator<(const Direction_2 &d) const
  {
    return R().compare_angle_with_x_axis_2_object()(*this, d) == SMALLER;
  }


  typename R::Boolean
  operator>(const Direction_2 &d) const
  {
    return d < *this;
  }


  typename R::Boolean
  operator>=(const Direction_2 &d) const
  {
    return R().compare_angle_with_x_axis_2_object()(*this, d) != SMALLER;
  }


  typename R::Boolean
  operator<=(const Direction_2 &d) const
  {
    return R().compare_angle_with_x_axis_2_object()(*this, d) != LARGER;
  }

  Direction_2
  operator-() const
  {
    return R().construct_opposite_direction_2_object()(*this);
  }

  Vector_2 vector() const
  {
    return R().construct_vector_2_object()(*this);
  }

  Vector_2 to_vector() const
  {
    return this->vector();
  }

  typename R::Boolean
  operator==(const Direction_2& d) const
  {
    return R().equal_2_object()(*this, d);
  }

  typename R::Boolean
  operator!=(const Direction_2& d) const
  {
    return !(*this == d);
  }

  Direction_2 transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }

};



template <class R >
std::ostream&
insert(std::ostream& os, const Direction_2<R>& d, const Cartesian_tag&)
{
    typename R::Vector_2 v = d.to_vector();
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << v.x() << ' ' << v.y();
    case IO::BINARY :
        write(os, v.x());
        write(os, v.y());
        return os;
    default:
        return os << "DirectionC2(" << v.x() << ", " << v.y() << ')';
    }
}

template <class R >
std::ostream&
insert(std::ostream& os, const Direction_2<R>& d, const Homogeneous_tag&)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << d.dx() << ' ' << d.dy();
    case IO::BINARY :
        write(os, d.dx());
        write(os, d.dy());
        return os;
    default:
        return os << "DirectionH2(" << d.dx() << ", "
                                    << d.dy() << ')';
  }
}

template < class R >
std::ostream&
operator<<(std::ostream& os, const Direction_2<R>& d)
{
  return insert(os, d, typename R::Kernel_tag() );
}


template <class R >
std::istream&
extract(std::istream& is, Direction_2<R>& d, const Cartesian_tag&)
{
    typename R::FT x, y;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> iformat(x) >> iformat(y);
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        std::cerr << std::endl << "Stream must be in ascii or binary mode"
                  << std::endl;
        break;
    }
    if (is)
        d = Direction_2<R>(x, y);
    return is;
}

template <class R >
std::istream&
extract(std::istream& is, Direction_2<R>& d, const Homogeneous_tag&)
{
  typename R::RT x, y;
  switch(is.iword(IO::mode))
  {
    case IO::ASCII :
        is >> iformat(x) >> iformat(y);
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  d = Direction_2<R>(x, y);
  return is;
}

template < class R >
std::istream&
operator>>(std::istream& is, Direction_2<R>& d)
{
  return extract(is, d, typename R::Kernel_tag() );
}

} //namespace CGAL

#endif // CGAL_DIRECTION_2_H
