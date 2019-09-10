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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Andreas Fabri
//                 Stefan Schirra

#ifndef CGAL_LINE_3_H
#define CGAL_LINE_3_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>
#include <CGAL/IO/io.h>

namespace CGAL {

template <class R_>
class Line_3 : public R_::Kernel_base::Line_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Plane_3               Plane_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;

  typedef Line_3                             Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Line_3>::value));

public:

  typedef Dimension_tag<3>  Ambient_dimension;
  typedef Dimension_tag<1>  Feature_dimension;

  typedef typename R_::Kernel_base::Line_3   Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  Line_3() {}

  // conversion impl -> interface class
  Line_3(const Rep& l)
      : Rep(l) {}

  Line_3(const Point_3 & p, const Point_3 & q)
      : Rep(typename R::Construct_line_3()(Return_base_tag(), p, q)) {}

  explicit Line_3(const Segment_3 & s)
      : Rep(typename R::Construct_line_3()(Return_base_tag(), s)) {}

  explicit Line_3(const Ray_3 & r)
      : Rep(typename R::Construct_line_3()(Return_base_tag(), r)) {}

  Line_3(const Point_3 & p, const Direction_3 & d)
      : Rep(typename R::Construct_line_3()(Return_base_tag(), p, d)) {}

  Line_3(const Point_3 & p, const Vector_3 & v)
      : Rep(typename R::Construct_line_3()(Return_base_tag(), p, v)) {}

  Line_3 transform(const Aff_transformation_3 &t) const
  {
    return Line_3(t.transform(this->point()), t.transform(this->direction()));
  }

  Vector_3 to_vector() const
  {
    return R().construct_vector_3_object()(*this);
  }

  Direction_3 direction() const
  {
    return R().construct_direction_3_object()(*this);
  }

  bool has_on(const Point_3 &p) const 
  { 
    return R().has_on_3_object()(*this, p);
    //return has_on_boundary(p); 
  }

  Point_3 point() const
  { 
    return R().construct_point_on_3_object()(*this, 0);
  }

  Point_3 point(int i) const
  { 
    return R().construct_point_on_3_object()(*this, i);
  }

  Line_3 opposite() const
  {
    return R().construct_opposite_line_3_object()(*this);
  }

  Plane_3 perpendicular_plane(const Point_3 &p) const
  {
    return R().construct_perpendicular_plane_3_object()(*this, p);
  }

  Point_3 projection(const Point_3 &p) const
  {
    return R().construct_projected_point_3_object()(*this, p);
  }

  bool is_degenerate() const
  {
    return R().is_degenerate_3_object()(*this);
  }
  
};

template < class R >
std::ostream &
operator<<(std::ostream &os, const Line_3<R> &l)
{
    switch(get_mode(os)) {
    case IO::ASCII :
        return os << l.point(0) << ' ' << l.point(1);
    case IO::BINARY :
        return os << l.point(0) <<  l.point(1);
    default:
        return  os << "Line_3(" << l.point(0) << ", " << l.point(1) << ")";
    }
}

template < class R >
std::istream &
operator>>(std::istream &is, Line_3<R> &l)
{
    typename R::Point_3 p, q;
    is >> p >> q;
    if (is)
        l = Line_3<R>(p, q);
    return is;
}

} //namespace CGAL

#endif // CGAL_LINE_3_H
