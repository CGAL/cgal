// Copyright (c) 2014  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Jocelyn Meyron
//

#ifndef CGAL_CONVEX_HULL_TRAITS_DUAL_3_H
#define CGAL_CONVEX_HULL_TRAITS_DUAL_3_H

#include <CGAL/license/Convex_hull_3.h>


#include <CGAL/Convex_hull_3/dual/predicates.h>
#include <CGAL/Convex_hull_3/dual/Convex_hull_traits_dual_2.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Cartesian_converter.h>

// Traits class used during the computation of the dual convex
// hull for the intersection of halfspaces.

// This traits class base on the concept ConvexHullTraits_3

namespace CGAL
{
  namespace Convex_hull_3
  {
    // Base traits class for dual predicates
    template <class R_>
      class Convex_hull_traits_base_dual_3
      {
        private:
          // Origin : point inside the halfspaces intersection
          typedef typename R_::Point_3 Primal_point_3;
          Primal_point_3 origin;

        public:

          Convex_hull_traits_base_dual_3 (Primal_point_3 o =
                                          Primal_point_3(0, 0, 0)) : origin(o)
          {}

          // Types
          typedef R_                                     R;
          typedef Convex_hull_traits_base_dual_3<R>      Self;
          typedef typename R::RT                         RT;

          // Dual
          typedef typename R::Plane_3         Point_3;
          typedef Plane_dual<R>               Plane_3;
          typedef Segment_dual<R>             Segment_3;
          typedef Plane_dual<R>               Triangle_3;

          // Traits used by convex_hull_2
          typedef typename CGAL::Convex_hull_3::Traits_xy_dual<R> Traits_xy_3;
          typedef typename CGAL::Convex_hull_3::Traits_yz_dual<R> Traits_yz_3;
          typedef typename CGAL::Convex_hull_3::Traits_xz_dual<R> Traits_xz_3;

          // Construct objects
          // Segment_3
          class Construct_segment_3 {
            public:
              Segment_3 operator ()(const Point_3& p, const Point_3& q)
              {
                return Segment_3(p, q);
              }
          };

          // Triangle_3
          class Construct_triangle_3 {
            public:
              Triangle_3 operator ()(const Point_3& p,
                                     const Point_3& q,
                                     const Point_3& r)
              {
                return Triangle_3(p, q, r);
              }
          };

          // Plane_3
          class Construct_plane_3 {
              public:
                  Plane_3 operator ()(const Point_3& p,
                                      const Point_3& q,
                                      const Point_3& r)
                  {
                      return Plane_3(p,q,r);
                  }
          };

          // Predicates
          typedef Equal_3_dual_point<R>                     Equal_3;
          typedef Collinear_3_dual_point<R>                 Collinear_3;
          typedef Coplanar_3_dual_point<R>                  Coplanar_3;
          typedef Has_on_positive_side_3_dual_point<R>      Has_on_positive_side_3;
          typedef Less_distance_to_point_3_dual_point<R>    Less_distance_to_point_3;
          typedef Less_signed_distance_to_plane_3_dual_point<R> Less_signed_distance_to_plane_3;

          Construct_segment_3
              construct_segment_3_object() const
              { return Construct_segment_3(); }

          Construct_plane_3
              construct_plane_3_object() const
              { return Construct_plane_3(); }

          Construct_triangle_3
              construct_triangle_3_object() const
              { return Construct_triangle_3(); }

          Collinear_3
              collinear_3_object() const
              { return Collinear_3(origin); }

          Coplanar_3
              coplanar_3_object() const
              { return Coplanar_3(origin); }

          Less_distance_to_point_3
              less_distance_to_point_3_object() const
              { return Less_distance_to_point_3(origin); }

          Has_on_positive_side_3
              has_on_positive_side_3_object() const
              { return Has_on_positive_side_3(origin); }

          Equal_3
              equal_3_object() const
              { return Equal_3(origin); }

          Less_signed_distance_to_plane_3
              less_signed_distance_to_plane_3_object() const
              { return Less_signed_distance_to_plane_3(origin); }

      };

    // Non-filtered traits class
    template <class R_, bool Has_filtered_predicates = R_::Has_filtered_predicates >
        class Convex_hull_traits_dual_3
        : public Convex_hull_traits_base_dual_3<R_>
        {
            private:
                typedef typename R_::Point_3 Primal_point_3;

            public:
                Convex_hull_traits_dual_3 (Primal_point_3 o =
                                           Primal_point_3(0, 0, 0)) : Convex_hull_traits_base_dual_3<R_>(o)
                {}
        } ;

    // Converter for dual planes
    template <class K1, class K2>
        struct Cartesian_converter_dual : public Cartesian_converter<K1, K2>
    {
        using CGAL::Cartesian_converter<K1, K2>::operator();

        Plane_dual<K2> operator() (const Plane_dual<K1> &in) const
        {
            return Plane_dual<K2>(operator()(in.p1),
                                  operator()(in.p2),
                                  operator()(in.p3));
        }
    };

    // Filtered traits
    template <typename R_>
        class Convex_hull_filtered_traits_dual_3
        : public Convex_hull_traits_base_dual_3<R_>
        {
            private:
                // Origin
                typedef typename R_::Point_3 Primal_point_3;
                Primal_point_3 origin;

            public:
                Convex_hull_filtered_traits_dual_3 (Primal_point_3 o =
                                                    Primal_point_3(0, 0, 0)) : origin(o)
                {}

                // Exact traits is based on the exact kernel.
                typedef Convex_hull_traits_dual_3<typename R_::Exact_kernel>
                    Exact_traits;

                // Filtering traits is based on the filtering kernel.
                typedef Convex_hull_traits_dual_3<typename R_::Approximate_kernel>
                    Filtering_traits;

                // Converters
                typedef Cartesian_converter_dual<R_, typename R_::Exact_kernel>
                    Converter_exact_dual;
                typedef Cartesian_converter_dual<R_, typename R_::Approximate_kernel>
                    Converter_approx_dual;

                // Filtered predicates
                typedef Filtered_predicate<
                    typename Exact_traits::Equal_3,
                    typename Filtering_traits::Equal_3,
                    Converter_exact_dual ,
                    Converter_approx_dual > Equal_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Collinear_3,
                    typename Filtering_traits::Collinear_3,
                    Converter_exact_dual,
                    Converter_approx_dual > Collinear_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Coplanar_3,
                    typename Filtering_traits::Coplanar_3,
                    Converter_exact_dual,
                    Converter_approx_dual > Coplanar_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Less_distance_to_point_3,
                    typename Filtering_traits::Less_distance_to_point_3,
                    Converter_exact_dual,
                    Converter_approx_dual > Less_distance_to_point_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Has_on_positive_side_3,
                    typename Filtering_traits::Has_on_positive_side_3,
                    Converter_exact_dual,
                    Converter_approx_dual > Has_on_positive_side_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Less_signed_distance_to_plane_3,
                    typename Filtering_traits::Less_signed_distance_to_plane_3,
                    Converter_exact_dual,
                    Converter_approx_dual > Less_signed_distance_to_plane_3;

                Collinear_3 collinear_3_object() const
                { return Collinear_3(origin); }

                Coplanar_3 coplanar_3_object() const
                { return Coplanar_3(origin); }

                Less_distance_to_point_3 less_distance_to_point_3_object() const
                { return Less_distance_to_point_3(origin); }

                Equal_3 equal_3_object() const
                { return Equal_3(origin); }

                Has_on_positive_side_3 has_on_positive_side_3_object() const
                { return Has_on_positive_side_3(origin); }

                Less_signed_distance_to_plane_3
                less_signed_distance_to_plane_3_object() const
                { return Less_signed_distance_to_plane_3(origin); }

                // Constructions are inherited
        };

    // Traits specialization
    template <typename R_>
        class Convex_hull_traits_dual_3<R_, true>
        : public Convex_hull_filtered_traits_dual_3<R_>
        {
            private:
                typedef typename R_::Point_3 Primal_point_3;

            public:
                Convex_hull_traits_dual_3 (Primal_point_3 o =
                                           Primal_point_3(0, 0, 0)) : Convex_hull_filtered_traits_dual_3<R_>(o)
                {}
        } ;
  } // namespace Convex_hull_3
} // namespace CGAL

#endif // CGAL_CONVEX_HULL_TRAITS_DUAL_3_H

