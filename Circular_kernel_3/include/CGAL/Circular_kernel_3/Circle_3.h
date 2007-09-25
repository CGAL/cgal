// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//             Pedro Machado    <tashimir@gmail.com>
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_SPHERICAL_KERNEL_CIRCLE_3_H
#define CGAL_SPHERICAL_KERNEL_CIRCLE_3_H

#include <CGAL/Circular_kernel_3/internal_functions_on_sphere_3.h>

namespace CGAL {
  namespace CGALi{
    template <class SK,class Sphere= typename CGAL::Sphere_3<SK> > class Circle_3 {

      typedef typename SK::Plane_3      Plane_3;
      //~ typedef typename SK::Sphere_3     Sphere_3;
      typedef Sphere                    Sphere_3;
      typedef typename SK::Point_3      Point_3;
      typedef typename SK::Vector_3     Vector_3;
      typedef typename SK::Direction_3  Direction_3;
      typedef typename SK::FT           FT;

    private:
      typedef std::pair<Sphere_3, Plane_3>  Rep;
      typedef typename SK::template Handle<Rep>::type  Base;

      Base base;

    public:
      Circle_3() {}

      Circle_3(const Point_3& center, const FT& squared_r, const Direction_3& d) 
      {
        // It is not allowed non-positive radius 
        // Should we keep this pre-condition?
        CGAL_kernel_assertion(squared_r > 0);
        base = Rep(Sphere_3(center,squared_r), 
                   plane_from_point_direction(center, d));
      }

      Circle_3(const Point_3& center, const FT& squared_r, const Vector_3& normal) 
      {
        // It is not allowed non-positive radius 
        // Should we keep this pre-condition?
        CGAL_kernel_assertion(squared_r > 0);
        base = Rep(Sphere_3(center,squared_r), 
                   plane_from_point_direction(center, normal.direction()));
      }

      Circle_3(const Point_3& center, const FT& squared_r, const Plane_3& p)
      {
        // the plane contains the center
  CGAL_kernel_assertion((p.a() * center.x() +
                               p.b() * center.y() +
                               p.c() * center.z() +
                               p.d()) == ZERO);
        // It is not allowed non-positive radius 
        // Should we keep this pre-condition?
        CGAL_kernel_assertion(squared_r > 0);
        base = Rep(Sphere_3(center,squared_r), p);
      }

      Circle_3(const Sphere_3 &s1, 
               const Sphere_3 &s2) {
         std::vector<Object> sols;
         SK().intersect_3_object()(s1, s2, std::back_inserter(sols));
         // s1,s2 must intersect
         CGAL_kernel_precondition(sols.size() != 0);
         typename SK::Circle_3 circle;
         // the intersection must be a circle (no point allowed)
         CGAL_kernel_precondition(assign(circle,sols[0]));
         assign(circle,sols[0]);
         *this = circle.rep();
      }

      Circle_3(const Plane_3 &p, 
               const Sphere_3 &s) {
         std::vector<Object> sols;
         SK().intersect_3_object()(p, s, std::back_inserter(sols));
         // s1,s2 must intersect
         CGAL_kernel_precondition(sols.size() != 0);
         typename SK::Circle_3 circle;
         // the intersection must be a circle (no point allowed)
         CGAL_kernel_precondition(assign(circle,sols[0]));
         assign(circle,sols[0]);
         *this = circle.rep();
      }

      const Plane_3& supporting_plane() const {
        return get(base).second;
      }

      Point_3 center() const {
        return diametral_sphere().center();
      }

      FT squared_radius() const {
        return diametral_sphere().squared_radius();
      }

      const Sphere_3& diametral_sphere() const {
        return get(base).first;
      }

      FT area_divided_by_pi() const {
        return squared_radius();
      }

      FT squared_length_divided_by_pi_square() const {
        return 4 * squared_radius();
      }

      static double pi;

      double approximate_area() const {
        return pi * to_double(squared_radius());
      }

      double approximate_squared_length() const {
        return pi * pi * 4.0 * to_double(squared_radius());
      }

      // this bbox function
      // can be optimize by doing different cases
      // for each variable = 0 (cases with is_zero)
      CGAL::Bbox_3 bbox() const {
        typedef CGAL::Interval_nt<false> Interval;
        CGAL::Interval_nt<false>::Protector ip;
        const Plane_3 &plane = supporting_plane();
        const Sphere_3 &s = diametral_sphere();
        const Point_3 &p = s.center();
        const FT &sq_r = s.squared_radius();
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
                                           // the multiplication would lead to inf results
        const Interval a2 = CGAL::square(a);
        const Interval b2 = CGAL::square(b);
        const Interval c2 = CGAL::square(c);
        const Interval sqr_sum = a2 + b2 + c2;
        const Interval mx = r * CGAL::sqrt((sqr_sum - a2)/sqr_sum);
        const Interval my = r * CGAL::sqrt((sqr_sum - b2)/sqr_sum);
        const Interval mz = r * CGAL::sqrt((sqr_sum - c2)/sqr_sum);
        return CGAL::Bbox_3((x-mx).inf(),(y-my).inf(), (z-mz).inf(),
                            (x+mx).sup(),(y+my).sup(), (z+mz).sup());
      }

      bool operator==(const Circle_3 &) const;
      bool operator!=(const Circle_3 &) const;

    };

    template < class SK, class Sphere >
    double Circle_3<SK,Sphere>::pi = CGAL_PI;
  //#PB : operateur for Circle_3 and Circle_on_sphere_3
    template < class SK, class Sphere >
    CGAL_KERNEL_INLINE
    bool
    Circle_3<SK,Sphere>::operator==(const Circle_3<SK,Sphere> &t) const
    {
      if (CGAL::identical(base, t.base))
        return true;
      return CGAL::SphericalFunctors::non_oriented_equal<SK>(*this, t);
    } 

    template < class SK, class Sphere >
    CGAL_KERNEL_INLINE
    bool
    Circle_3<SK,Sphere>::operator!=(const Circle_3<SK,Sphere> &t) const
    {
      return !(*this == t);
    }

  }
}

#endif
