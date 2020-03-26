// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado    <tashimir@gmail.com>
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_CIRCLE_3_H
#define CGAL_CIRCLE_3_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R_>
  class Circle_3
  : public R_::Kernel_base::Circle_3
{
  typedef typename R_::RT                    RT;
  typedef typename R_::FT                    FT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Plane_3               Plane_3;
  typedef typename R_::Vector_3              Vector_3;
  typedef typename R_::Sphere_3              Sphere_3;
  typedef typename R_::Direction_3           Direction_3;

  typedef Circle_3                           Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Circle_3>::value));

public:

  typedef Dimension_tag<3>  Ambient_dimension;
  typedef Dimension_tag<1>  Feature_dimension;

  typedef typename R_::Kernel_base::Circle_3 Rep;
  typedef R_                                 R;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  Circle_3() {}

  Circle_3(const Sphere_3& s, const Plane_3& p, int a)
    : Rep(typename R::Construct_circle_3()(s,p,a)) {}

  Circle_3(const Point_3& c, const FT& sr, const Plane_3& p)
    : Rep(typename R::Construct_circle_3()(c,sr,p)) {}

  Circle_3(const Point_3& c, const FT& sr, const Direction_3& d)
    : Rep(typename R::Construct_circle_3()(c,sr,d)) {}

  Circle_3(const Point_3& c, const FT& sr, const Vector_3& v)
    : Rep(typename R::Construct_circle_3()(c,sr,v)) {}

  Circle_3(const Sphere_3& s1, const Sphere_3& s2)
    : Rep(typename R::Construct_circle_3()(s1,s2)) {}

  Circle_3(const Sphere_3& s, const Plane_3& p)
    : Rep(typename R::Construct_circle_3()(s,p)) {}

  Circle_3(const Plane_3& p, const Sphere_3& s)
    : Rep(typename R::Construct_circle_3()(p,s)) {}

  Circle_3(const Point_3& p1, const Point_3& p2, const Point_3& p3)
    : Rep(typename R::Construct_circle_3()(p1,p2,p3)) {}

  Circle_3(const Rep& r)
    : Rep(r) {}

  typename cpp11::result_of
  <typename R::Construct_sphere_3( Circle_3)>::type
  diametral_sphere() const
  {
    return typename R::Construct_sphere_3()(*this);
  }

  Point_3 center() const
  {
    return typename R::Construct_sphere_3()(*this).center();
  }

  FT squared_radius() const
  {
    return typename R::Construct_sphere_3()(*this).squared_radius();
  }

  typename cpp11::result_of
  <typename R::Construct_plane_3( Circle_3)>::type
  supporting_plane() const
  {
    return typename R::Construct_plane_3()(*this);
  }

  Bbox_3 bbox() const
  {
    return typename R::Construct_bbox_3()(*this);
  }

        FT area_divided_by_pi() const
        {
          return typename R::Compute_area_divided_by_pi_3()(*this);
  }

  double approximate_area() const
  {
          return typename R::Compute_approximate_area_3()(*this);
        }

        FT squared_length_divided_by_pi_square() const
        {
          return typename R::Compute_squared_length_divided_by_pi_square_3()(*this);
  }

  double approximate_squared_length() const
  {
          return typename R::Compute_approximate_squared_length_3()(*this);
        }

        typename R::Boolean
  has_on(const Point_3 &p) const
  {
    return typename R::Has_on_3()(*this, p);
  }

};

template < typename R >
inline
bool
operator==(const Circle_3<R> &p,
           const Circle_3<R> &q)
{
  return R().equal_3_object()(p, q);
}

template < typename R >
inline
bool
operator!=(const Circle_3<R> &p,
           const Circle_3<R> &q)
{
  return ! (p == q);
}

template < typename R >
std::ostream &
operator<<(std::ostream & os, const Circle_3<R> &c)
{
  return os << c.supporting_plane() << " "
    << c.diametral_sphere() << " ";
}

template < typename R >
std::istream &
operator>>(std::istream & is, Circle_3<R> &c)
{
  typename R::Plane_3 p;
  typename R::Sphere_3 s;

  is >> p >> s ;
  if (is)
    c = Circle_3<R>(p, s);
  return is;
}

} //namespace CGAL

#endif // CGAL_CIRCLE_3_H
