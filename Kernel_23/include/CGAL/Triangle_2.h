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
// Author(s)     : Andreas Fabri

#ifndef CGAL_TRIANGLE_2_H
#define CGAL_TRIANGLE_2_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R_>
class Triangle_2 : public R_::Kernel_base::Triangle_2
{
  typedef typename R_::Point_2          Point_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;
  typedef typename R_::Kernel_base::Triangle_2  RTriangle_2;

  typedef Triangle_2                            Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Triangle_2>::value));

public:

  typedef Dimension_tag<2>  Ambient_dimension;
  typedef Dimension_tag<2>  Feature_dimension;

  typedef RTriangle_2 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_                          R;
  typedef typename R::FT               FT;

  Triangle_2() {}

  Triangle_2(const RTriangle_2& t)
      : RTriangle_2(t) {}

  Triangle_2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
      : RTriangle_2(typename R::Construct_triangle_2()(Return_base_tag(), p,q,r)) {}

  FT
  area() const
  {
    return R().compute_area_2_object()(vertex(0), vertex(1), vertex(2));
  }

  typename R::Orientation
  orientation() const
  {
    return R().orientation_2_object()(vertex(0), vertex(1), vertex(2));
  }

  typename R::Bounded_side
  bounded_side(const Point_2 &p) const
  {
    return R().bounded_side_2_object()(*this,p);
  }

  typename R::Oriented_side
  oriented_side(const Point_2 &p) const
  {
    return R().oriented_side_2_object()(*this,p);
  }

  typename R::Boolean
  operator==(const Triangle_2 &t) const
  {
    return R().equal_2_object()(*this,t);
  }

  typename R::Boolean
  operator!=(const Triangle_2 &t) const
  {
    return !(*this == t);
  }

  typename cpp11::result_of<typename R::Construct_vertex_2( Triangle_2, int)>::type
  vertex(int i) const
  {
    return R().construct_vertex_2_object()(*this,i);
  }

  typename cpp11::result_of<typename R::Construct_vertex_2( Triangle_2, int)>::type
  operator[](int i) const
  {
    return vertex(i);
  }

  typename R::Boolean
  has_on_bounded_side(const Point_2 &p) const
  {
    return bounded_side(p) == ON_BOUNDED_SIDE;
  }

  typename R::Boolean
  has_on_unbounded_side(const Point_2 &p) const
  {
    return bounded_side(p) == ON_UNBOUNDED_SIDE;
  }

  typename R::Boolean
  has_on_boundary(const Point_2 &p) const
  {
    return bounded_side(p) == ON_BOUNDARY;
  }

  typename R::Boolean
  has_on_negative_side(const Point_2 &p) const
  {
    return oriented_side(p) == ON_NEGATIVE_SIDE;
  }

  typename R::Boolean
  has_on_positive_side(const Point_2 &p) const
  {
    return oriented_side(p) == ON_POSITIVE_SIDE;
  }

  typename R::Boolean
  is_degenerate() const
  {
    return R().collinear_2_object()(vertex(0), vertex(1), vertex(2));
  }

  Bbox_2
  bbox() const
  {
    return R().construct_bbox_2_object()(*this);
  }

  Triangle_2
  opposite() const
  {
    return R().construct_opposite_triangle_2_object()(*this);
  }

  Triangle_2
  transform(const Aff_transformation_2 &t) const
  {
    return Triangle_2(t.transform(vertex(0)),
                      t.transform(vertex(1)),
                      t.transform(vertex(2)));
  }

};


template < class R >
std::ostream &
operator<<(std::ostream &os, const Triangle_2<R> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0] << t[1]  << t[2];
    default:
        return os<< "Triangle_2(" << t[0] << ", "
                 << t[1] << ", " << t[2] <<")";
    }
}

template < class R >
std::istream &
operator>>(std::istream &is, Triangle_2<R> &t)
{
    typename R::Point_2 p, q, r;

    is >> p >> q >> r;

    if (is)
        t = Triangle_2<R>(p, q, r);
    return is;
}

} //namespace CGAL

#endif // CGAL_TRIANGLE_2_H
