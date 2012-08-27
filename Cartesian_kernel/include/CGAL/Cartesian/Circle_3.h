// Copyright (c) 2000  
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
// Author(s)     : Monique Teillaud, Pedro Machado, Sebastien Loriot

#ifndef CGAL_CARTESIAN_CIRCLEC3_H
#define CGAL_CARTESIAN_CIRCLEC3_H

#include <CGAL/Interval_nt.h>

namespace CGAL {

template <class R_ >
class CircleC3 {
  typedef typename R_::Sphere_3                 Sphere_3;
  typedef typename R_::Plane_3                  Plane_3;
  typedef typename R_::Point_3                  Point_3;
  typedef typename R_::Vector_3                 Vector_3;
  typedef typename R_::Direction_3              Direction_3;
  typedef typename R_::FT                       FT;

  typedef std::pair<Sphere_3, Plane_3>             Rep;
  typedef typename R_::template Handle<Rep>::type  Base;
  Base base;  

public:
  typedef R_                                     R;

  CircleC3() {}

  CircleC3(const Point_3& center, const FT& squared_r, const Direction_3& d)
  {
    CGAL_kernel_assertion(squared_r >= FT(0));
    // non-degenerated Direction
    CGAL_kernel_assertion((d.dx() != FT(0)) || (d.dy() != FT(0)) || (d.dz() != FT(0)));
    base = Rep(Sphere_3(center,squared_r),
                plane_from_point_direction(center, d));
  }

  CircleC3(const Point_3& center, const FT& squared_r, const Vector_3& normal) 
  {
    CGAL_kernel_assertion(squared_r >= FT(0));
    // non-degenerated Vector
    CGAL_kernel_assertion((normal.x() != FT(0)) ||
                          (normal.y() != FT(0)) ||
                          (normal.z() != FT(0)));
    base = Rep(Sphere_3(center,squared_r),
                Plane_3(center, normal.direction()));
  }

  CircleC3(const Point_3& center, const FT& squared_r, const Plane_3& p)
  {
    // the plane contains the center and it is not degenerate
    CGAL_kernel_assertion(!R().is_degenerate_3_object()(p));
    CGAL_kernel_assertion((p.a() * center.x() +
                           p.b() * center.y() +
                           p.c() * center.z() +
                           p.d()) == CGAL::ZERO);
    CGAL_kernel_assertion(squared_r >= FT(0));
    base = Rep(Sphere_3(center,squared_r), p);
  }

  CircleC3(const Sphere_3 &s1, const Sphere_3 &s2) {
    Object obj = R().intersect_3_object()(s1, s2);
    // s1,s2 must intersect
    CGAL_kernel_precondition(!(obj.is_empty()));
    const typename R::Circle_3* circle_ptr=object_cast<typename R::Circle_3>(&obj);
    if(circle_ptr!=NULL)
      base = Rep(circle_ptr->diametral_sphere(), circle_ptr->supporting_plane());
    else {
      const typename R::Point_3* point=object_cast<typename R::Point_3>(&obj);
      CGAL_kernel_precondition(point!=NULL);
      CircleC3 circle = CircleC3(*point, FT(0), Vector_3(FT(1),FT(0),FT(0)));
      base = Rep(circle.diametral_sphere(), circle.supporting_plane());
    }
  }

  CircleC3(const Plane_3 &p, const Sphere_3 &s, int) : base(s, p) {}

  CircleC3(const Plane_3 &p, const Sphere_3 &s) {
    Object obj = R().intersect_3_object()(p, s);
    // s1,s2 must intersect
    CGAL_kernel_precondition(!(obj.is_empty()));
    const typename R::Circle_3* circle_ptr=object_cast<typename R::Circle_3>(&obj);
    if(circle_ptr!=NULL)
      base = Rep(circle_ptr->diametral_sphere(), circle_ptr->supporting_plane());
    else {
      const typename R::Point_3* point=object_cast<typename R::Point_3>(&obj);
      CGAL_kernel_precondition(point!=NULL);
      CircleC3 circle = CircleC3(*point, FT(0), Vector_3(FT(1),FT(0),FT(0)));
      base = Rep(circle.diametral_sphere(), circle.supporting_plane());
    }
  }

  CircleC3(const Point_3 &p, const Point_3 &q, const Point_3 &r) {
	  // p, q, r are not collinear
	  CGAL_kernel_precondition(!R().collinear_3_object()(p, q, r));
		Plane_3 p1 = R().construct_plane_3_object()(p, q, r);
    Plane_3 p2 = R().construct_bisector_3_object()(p, q);
    Plane_3 p3 = R().construct_bisector_3_object()(p, r);
    Object obj = R().intersect_3_object()(p1, p2, p3);
    // must be a point, otherwise they are collinear
    const Point_3& center=*object_cast<Point_3>(&obj);
		FT sqr = R().compute_squared_distance_3_object()(center, r);
		Sphere_3 s = R().construct_sphere_3_object()(center, sqr);
		base = Rep(s, p1);
  }

  const Plane_3& supporting_plane() const
  {
    return get(base).second;
  }

  const Sphere_3& supporting_sphere() const
  {
    return diametral_sphere();
  }

  Point_3 center() const
  {
    return diametral_sphere().center();
  }

  FT squared_radius() const
  {
    return diametral_sphere().squared_radius();
  }

  const Sphere_3& diametral_sphere() const
  {
    return get(base).first;
  }

  double approximate_area() const
  {
    return CGAL_PI * to_double(squared_radius());
  }

  double approximate_squared_length() const
  {
    return CGAL_PI * CGAL_PI * 4.0 * to_double(squared_radius());
  }

  FT area_divided_by_pi() const
  {
    return squared_radius();
  }

  FT squared_length_divided_by_pi_square() const
  {
    return 4 * squared_radius();
  }

  // this bbox function
  // can be optimize by doing different cases
  // for each variable = 0 (cases with is_zero)
  CGAL::Bbox_3 bbox() const
  {
    typedef CGAL::Interval_nt<false> Interval;
    CGAL::Interval_nt<false>::Protector ip;
    const Sphere_3 &s = diametral_sphere();
    const FT &sq_r = s.squared_radius();
    const Point_3 &p = s.center();
    if(sq_r == FT(0)) return p.bbox();
    const Plane_3 &plane = supporting_plane();
    const Interval a = CGAL::to_interval(plane.a());
    const Interval b = CGAL::to_interval(plane.b());
    const Interval c = CGAL::to_interval(plane.c());
    const Interval x = CGAL::to_interval(p.x());
    const Interval y = CGAL::to_interval(p.y());
    const Interval z = CGAL::to_interval(p.z());
    const Interval r2 = CGAL::to_interval(sq_r);
    const Interval r = CGAL::sqrt(r2); // maybe we can work with r2
                                       // in order to save this operation
                                       // but if the coefficients are to high
                                       // the multiplication would lead to inf
                                       // results
    const Interval a2 = CGAL::square(a);
    const Interval b2 = CGAL::square(b);
    const Interval c2 = CGAL::square(c);
    const Interval sqr_sum = a2 + b2 + c2;
    const Interval mx = r * CGAL::sqrt((sqr_sum - a2)/sqr_sum);
    const Interval my = r * CGAL::sqrt((sqr_sum - b2)/sqr_sum);
    const Interval mz = r * CGAL::sqrt((sqr_sum - c2)/sqr_sum);
    return CGAL::Bbox_3((x-mx).inf(),(y-my).inf(),(z-mz).inf(),
                        (x+mx).sup(),(y+my).sup(),(z+mz).sup());
  }

  bool operator==(const CircleC3 &) const;
  bool operator!=(const CircleC3 &) const;

  bool has_on(const Point_3 &p) const;
  bool has_on_bounded_side(const Point_3 &p) const;
  bool has_on_unbounded_side(const Point_3 &p) const;
  Bounded_side bounded_side(const Point_3 &p) const;

  bool is_degenerate() const
  {
    return diametral_sphere().is_degenerate();
  }

};

template < class R >
inline
bool
CircleC3<R>::
has_on(const typename CircleC3<R>::Point_3 &p) const
{
  return R().has_on_3_object()(diametral_sphere(),p) &&
         R().has_on_3_object()(supporting_plane(),p);
}

template < class R >
inline
bool
CircleC3<R>::
has_on_bounded_side(const typename CircleC3<R>::Point_3 &p) const
{
  CGAL_kernel_precondition(R().has_on_3_object()(supporting_plane(), p));
  return squared_distance(center(),p) < squared_radius();
}

template < class R >
inline
bool
CircleC3<R>::
has_on_unbounded_side(const typename CircleC3<R>::Point_3 &p) const
{
  CGAL_kernel_precondition(R().has_on_3_object()(supporting_plane(), p));
  return squared_distance(center(),p) > squared_radius();
}

template < class R >
CGAL_KERNEL_INLINE
Bounded_side
CircleC3<R>::
bounded_side(const typename CircleC3<R>::Point_3 &p) const
{
  CGAL_kernel_precondition(is_degenerate() || R().has_on_3_object()(supporting_plane(), p));
  return diametral_sphere().bounded_side(p);
}

template < class R >
CGAL_KERNEL_INLINE
bool
CircleC3<R>::operator==(const CircleC3<R> &t) const
{
  if (CGAL::identical(base, t.base))
    return true;
  if(!(center() == t.center() &&
       squared_radius() == t.squared_radius())) return false;

  const typename R::Plane_3 p1 = supporting_plane();
  const typename R::Plane_3 p2 = t.supporting_plane();

  if(is_zero(p1.a())) {
    if(!is_zero(p2.a())) return false;
    if(is_zero(p1.b())) {
      if(!is_zero(p2.b())) return false;
      return p1.c() * p2.d() == p1.d() * p2.c();
    }
    return (p2.c() * p1.b() == p1.c() * p2.b()) &&
           (p2.d() * p1.b() == p1.d() * p2.b());
  }
  return (p2.b() * p1.a() == p1.b() * p2.a()) &&
         (p2.c() * p1.a() == p1.c() * p2.a()) &&
         (p2.d() * p1.a() == p1.d() * p2.a());
}

template < class R >
CGAL_KERNEL_INLINE
bool
CircleC3<R>::operator!=(const CircleC3<R> &t) const
{
  return !(*this == t);
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_CIRCLEC3_H
