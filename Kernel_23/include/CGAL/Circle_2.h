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
//                 Sven Schoenherr

#ifndef CGAL_CIRCLE_2_H
#define CGAL_CIRCLE_2_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Dimension.h>
#include <CGAL/number_utils.h>
#include <CGAL/result_of.h>

namespace CGAL {

template <class R_>
class Circle_2 : public R_::Kernel_base::Circle_2
{
  typedef typename R_::FT                    FT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Kernel_base::Circle_2 RCircle_2;
  typedef typename R_::Aff_transformation_2  Aff_transformation_2;

  typedef Circle_2                           Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Circle_2>::value));

public:

  typedef Dimension_tag<2>  Ambient_dimension;
  typedef Dimension_tag<1>  Feature_dimension;

  typedef RCircle_2 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef  R_   R;

  Circle_2() {}

  Circle_2(const RCircle_2& t)
    : RCircle_2(t) {}

  Circle_2(const Point_2 &center, const FT &squared_radius,
	   const Orientation &orientation)
    : RCircle_2(typename R::Construct_circle_2()(Return_base_tag(), center, squared_radius, orientation)) {}

  Circle_2(const Point_2 &center, const FT &squared_radius)
    : RCircle_2(typename R::Construct_circle_2()(Return_base_tag(), center, squared_radius, COUNTERCLOCKWISE)) {}

  Circle_2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
    : RCircle_2(typename R::Construct_circle_2()(Return_base_tag(), p, q, r)) {}

  Circle_2(const Point_2 & p, const Point_2 & q,
	   const Orientation &orientation)
    : RCircle_2(typename R::Construct_circle_2()(Return_base_tag(), p, q, orientation)) {}

  Circle_2(const Point_2 & p, const Point_2 & q)
    : RCircle_2(typename R::Construct_circle_2()(Return_base_tag(), p, q, COUNTERCLOCKWISE)) {}

  Circle_2(const Point_2 & p, const Point_2 & q, const FT &bulge)
    : RCircle_2(typename R::Construct_circle_2()(Return_base_tag(), p, q, bulge)) {}

  Circle_2(const Point_2 & center, const Orientation& orientation)
    : RCircle_2(typename R::Construct_circle_2()(Return_base_tag(), center, FT(0), orientation)) {}

  Circle_2(const Point_2 & center)
    : RCircle_2(typename R::Construct_circle_2()(Return_base_tag(), center, FT(0), COUNTERCLOCKWISE)) {}

  typename cpp11::result_of<typename R::Construct_center_2(Circle_2)>::type
  center() const
  {
    return R().construct_center_2_object()(*this);
  }

  typename cpp11::result_of<typename R::Compute_squared_radius_2(Circle_2)>::type
  squared_radius() const
  {
    return R().compute_squared_radius_2_object()(*this);
  }

  Orientation orientation() const
  {
    // This make_certain(), the uncertain orientation of circles, the orientation
    // of circles, are all yucky.
    return make_certain(R().orientation_2_object()(*this));
  }


  typename R::Bounded_side
  bounded_side(const Point_2 &p) const
  {
    return R().bounded_side_2_object()(*this, p);
  }

  typename R::Oriented_side
  oriented_side(const Point_2 &p) const
  {
    return R().oriented_side_2_object()(*this, p);
  }

  typename R::Boolean
  has_on_boundary(const Point_2 &p) const
  {
    return bounded_side(p) == ON_BOUNDARY;
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
  has_on_negative_side(const Point_2 &p) const
  {
    if (orientation() == COUNTERCLOCKWISE)
      return has_on_unbounded_side(p);
    return has_on_bounded_side(p);
  }

  typename R::Boolean
  has_on_positive_side(const Point_2 &p) const
  {
    if (orientation() == COUNTERCLOCKWISE)
      return has_on_bounded_side(p);
    return has_on_unbounded_side(p);
  }

  typename R::Boolean
  is_degenerate() const
  {
    return CGAL_NTS is_zero(squared_radius());
  }

  Circle_2
  opposite() const
  {
    //return R().construct_opposite_circle_2_object()(*this);
    return Circle_2(center(),
		    squared_radius(),
		    CGAL::opposite(orientation()) );
  }

  Bbox_2
  bbox() const
  {
    return R().construct_bbox_2_object()(*this);
  }

  typename R::Boolean
  operator==(const Circle_2 &c) const
  {
    return R().equal_2_object()(*this, c);
  }

  typename R::Boolean
  operator!=(const Circle_2 &c) const
  {
    return !(*this == c);
  }

  Circle_2 transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }

  Circle_2 orthogonal_transform(const Aff_transformation_2 &t) const;


};

template <class R_>
Circle_2<R_>
Circle_2<R_>::
orthogonal_transform(const typename R_::Aff_transformation_2& t) const
{
  typedef typename  R_::RT  RT;
  typedef typename  R_::FT  FT;
  typedef typename  R_::Vector_2  Vector_2;

  Vector_2 vec(RT(1), RT(0) );   // unit vector // AF: was FT
  vec = vec.transform(t);             // transformed
  FT sq_scale = vec.squared_length();       // squared scaling factor

  return Circle_2(t.transform(center()),
                               sq_scale * squared_radius(),
                               t.is_even() ? orientation()
                                           : CGAL::opposite(orientation()));
}


template <class R >
std::ostream&
insert(std::ostream& os, const Circle_2<R>& c)
{
    switch(get_mode(os)) {
    case IO::ASCII :
        os << c.center() << ' ' << c.squared_radius() << ' '
           << static_cast<int>(c.orientation());
        break;
    case IO::BINARY :
        os << c.center();
        write(os, c.squared_radius());
        write(os, static_cast<int>(c.orientation()));
        break;
    default:
        os << "Circle_2(" << c.center() <<  ", " << c.squared_radius() ;
        switch (c.orientation()) {
        case CLOCKWISE:
            os << ", clockwise)";
            break;
        case COUNTERCLOCKWISE:
            os << ", counterclockwise)";
            break;
        default:
            os << ", collinear)";
            break;
        }
        break;
    }
    return os;
}

template < class R >
std::ostream &
operator<<(std::ostream &os, const Circle_2<R> &c)
{
  return insert(os, c);
}


template <class R >
std::istream&
extract(std::istream& is, Circle_2<R>& c)
{
    typename R::Point_2 center;
    typename R::FT squared_radius(0);
    int o=0;
    switch(get_mode(is)) {
    case IO::ASCII :
        is >> center >> iformat(squared_radius) >> o;
        break;
    case IO::BINARY :
        is >> center;
        read(is, squared_radius);
        is >> o;
        break;
    default:
        is.setstate(std::ios::failbit);
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	c = Circle_2<R>(center, squared_radius, static_cast<Orientation>(o));
    return is;
}

template < class R >
std::istream &
operator>>(std::istream &is, Circle_2<R> &c)
{
  return extract(is,c);
}

} //namespace CGAL

#endif  // CGAL_CIRCLE_2_H
