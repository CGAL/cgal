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

#ifndef CGAL_CH3_DUAL_PREDICATES_H
#define CGAL_CH3_DUAL_PREDICATES_H

#include <CGAL/license/Convex_hull_3.h>


#include <CGAL/predicates/sign_of_determinant.h>

// Predicates used during the computation of the dual convex hull
// The dual point associated to the plane which  does not contain the origin
// and whose equation is
// a x + b y + c z + d = 0
// is the point (- a / d, - b / d, - c / d)

// The predicates are the ones required by the ConvexHullTraits_3 concept

namespace CGAL
{
    namespace Convex_hull_3
    {
        // Plane in the dual space : 3 dual points = 3 planes
        template <typename K>
            struct Plane_dual
            {
                typedef typename K::Plane_3 Plane_3;

                Plane_3 p1;
                Plane_3 p2;
                Plane_3 p3;

                Plane_dual (Plane_3 p1, Plane_3 p2, Plane_3 p3) :
                    p1(p1), p2(p2), p3(p3)
                {}

                Plane_dual ()
                {}
            };

        // Segment in the dual space : 2 dual points = 2 planes
        template <typename K>
            struct Segment_dual
            {
                typedef typename K::Plane_3 Plane_3;

                Plane_3 p;
                Plane_3 q;

                Segment_dual (Plane_3 p, Plane_3 q) : p(p), q(q)
                {}
            };

        // Predicates for dual points
        // Equal
        template < typename K >
            struct Equal_3_dual_point
            {
                typedef typename K::RT        RT;
                typedef typename K::Plane_3   Plane_3;
                typedef typename K::Point_3   Point_3;
                typedef bool                  result_type;

                Equal_3_dual_point (Point_3 const& o = Point_3(0, 0, 0)) :
                    origin(o)
                {}

                result_type
                    operator()(const Plane_3 &p, const Plane_3 &q) const
                    {
                        RT dp = p.d() + origin.x() * p.a()
                            + origin.y() * p.b() + origin.z() * p.c();
                        RT dq = q.d() + origin.x() * q.a()
                            + origin.y() * q.b() + origin.z() * q.c();

                        RT diffa = p.a() * dq - dp * q.a();
                        RT diffb = p.b() * dq - dp * q.b();
                        RT diffc = p.c() * dq - dp * q.c();

                        return CGAL_AND_3(CGAL::is_zero(diffa),
                                          CGAL::is_zero(diffb),
                                          CGAL::is_zero(diffc));
                    }

                private:
                    Point_3 origin;
            };

        // Collinear
        template < typename K >
            struct Collinear_3_dual_point
            {
                typedef typename K::RT        RT;
                typedef typename K::Plane_3   Plane_3;
                typedef typename K::Point_3   Point_3;
                typedef bool                  result_type;

                Collinear_3_dual_point (Point_3 const& o = Point_3(0, 0, 0)) :
                    origin(o)
                {}

                result_type
                    operator()(const Plane_3 &p,
                               const Plane_3 &q,
                               const Plane_3 &r) const
                    {
                        RT dp = p.d() + origin.x() * p.a()
                            + origin.y() * p.b() + origin.z() * p.c();
                        RT dq = q.d() + origin.x() * q.a()
                            + origin.y() * q.b() + origin.z() * q.c();
                        RT dr = r.d() + origin.x() * r.a()
                            + origin.y() * r.b() + origin.z() * r.c();

                        RT diffapq = dp * q.a() - dq * p.a();
                        RT diffbpq = dp * q.b() - dq * p.b();
                        RT diffcpq = dp * q.c() - dq * p.c();

                        RT diffapr = dp * r.a() - dr * p.a();
                        RT diffbpr = dp * r.b() - dr * p.b();
                        RT diffcpr = dp * r.c() - dr * p.c();

                        // Cross product
                        RT cross1 = diffbpq * diffcpr - diffcpq * diffbpr;
                        RT cross2 = diffcpq * diffapr - diffapq * diffcpr;
                        RT cross3 = diffapq * diffbpr - diffbpq * diffapr;

                        return CGAL_AND_3(CGAL::is_zero(cross1),
                                          CGAL::is_zero(cross2),
                                          CGAL::is_zero(cross3));
                    }

                private:
                    Point_3 origin;
            };

        // Coplanar
        template < typename K >
            struct Coplanar_3_dual_point
            {
                typedef typename K::RT        RT;
                typedef typename K::Plane_3   Plane_3;
                typedef typename K::Point_3   Point_3;
                typedef bool                  result_type;

                Coplanar_3_dual_point (Point_3 const& o = Point_3(0, 0, 0)) :
                    origin(o)
                {}

                result_type
                    operator()(const Plane_3 &p,
                               const Plane_3 &q,
                               const Plane_3 &r,
                               const Plane_3 &s) const
                    {
                        RT dp = p.d() + origin.x() * p.a()
                            + origin.y() * p.b() + origin.z() * p.c();
                        RT dq = q.d() + origin.x() * q.a()
                            + origin.y() * q.b() + origin.z() * q.c();
                        RT dr = r.d() + origin.x() * r.a()
                            + origin.y() * r.b() + origin.z() * r.c();
                        RT ds = s.d() + origin.x() * s.a()
                            + origin.y() * s.b() + origin.z() * s.c();

                        RT diffapq = dp * q.a() - dq * p.a();
                        RT diffbpq = dp * q.b() - dq * p.b();
                        RT diffcpq = dp * q.c() - dq * p.c();

                        RT diffapr = dp * r.a() - dr * p.a();
                        RT diffbpr = dp * r.b() - dr * p.b();
                        RT diffcpr = dp * r.c() - dr * p.c();

                        RT diffaps = dp * s.a() - ds * p.a();
                        RT diffbps = dp * s.b() - ds * p.b();
                        RT diffcps = dp * s.c() - ds * p.c();

                        return (CGAL::sign_of_determinant(diffapq, diffapr, diffaps,
                                                          diffbpq, diffbpr, diffbps,
                                                          diffcpq, diffcpr, diffcps)
                            == CGAL::ZERO);
                    }

                private:
                    Point_3 origin;
            };

        // Has on positive side
        template < typename K >
            struct Has_on_positive_side_3_dual_point
            {
                typedef typename K::RT         RT;
                typedef typename K::Plane_3    Plane_3;
                typedef typename K::Point_3   Point_3;
                typedef bool                   result_type;

                Has_on_positive_side_3_dual_point (Point_3 const& o =
                                                   Point_3(0, 0, 0)) :
                    origin(o)
                {}

                result_type
                    operator()(const Plane_dual<K> p, const Plane_3 &q) const
                    {
                        Plane_3 p1 = p.p1;
                        Plane_3 p2 = p.p2;
                        Plane_3 p3 = p.p3;

                        RT dp1 = p1.d() + origin.x() * p1.a()
                            + origin.y() * p1.b() + origin.z() * p1.c();
                        RT dp2 = p2.d() + origin.x() * p2.a()
                            + origin.y() * p2.b() + origin.z() * p2.c();
                        RT dp3 = p3.d() + origin.x() * p3.a()
                            + origin.y() * p3.b() + origin.z() * p3.c();
                        RT dq = q.d() + origin.x() * q.a()
                            + origin.y() * q.b() + origin.z() * q.c();

                        // Compute the normal to the plane
                        RT alpha = (dp1 * p2.b() - dp2 * p1.b()) *
                            (dp1 * p3.c() - dp3 * p1.c()) -
                            (dp1 * p2.c() - dp2 * p1.c()) *
                            (dp1 * p3.b() - dp3 * p1.b());

                        RT beta  = (dp1 * p2.c() - dp2 * p1.c()) *
                            (dp1 * p3.a() - dp3 * p1.a()) -
                            (dp1 * p2.a() - dp2 * p1.a()) *
                            (dp1 * p3.c() - dp3 * p1.c());

                        RT gamma = (dp1 * p2.a() - dp2 * p1.a()) *
                            (dp1 * p3.b() - dp3 * p1.b()) -
                            (dp1 * p2.b() - dp2 * p1.b()) *
                            (dp1 * p3.a() - dp3 * p1.a());

                        // Test if q is on the positive side of p
                        RT prod = alpha * (p1.a() * dq - q.a() * dp1) +
                            beta * (p1.b() * dq - q.b() * dp1) +
                            gamma * (p1.c() * dq - q.c() * dp1);

                        if (CGAL::is_positive(dp1 * dq)) {
                            return CGAL::is_positive(prod);
                        } else {
                            return CGAL::is_negative(prod);
                        }
                    }

                private:
                    Point_3 origin;
            };

        // Less distance to point
        template < typename K >
            struct Less_distance_to_point_3_dual_point
            {
                typedef typename K::RT        RT;
                typedef typename K::Plane_3   Plane_3;
                typedef typename K::Point_3   Point_3;
                typedef bool                  result_type;

                Less_distance_to_point_3_dual_point (Point_3 const& o =
                                                     Point_3(0, 0, 0)) :
                    origin(o)
                {}

                result_type
                    operator()(const Plane_3 &p,
                               const Plane_3 &q,
                               const Plane_3 &r) const
                    {
                        RT dp = p.d() + origin.x() * p.a()
                            + origin.y() * p.b() + origin.z() * p.c();
                        RT dq = q.d() + origin.x() * q.a()
                            + origin.y() * q.b() + origin.z() * q.c();
                        RT dr = r.d() + origin.x() * r.a()
                            + origin.y() * r.b() + origin.z() * r.c();

                        RT diffapq = p.a() * dq - q.a() * dp;
                        RT diffbpq = p.b() * dq - q.b() * dp;
                        RT diffcpq = p.c() * dq - q.c() * dp;

                        RT diffapr = p.a() * dr - r.a() * dp;
                        RT diffbpr = p.b() * dr - r.b() * dp;
                        RT diffcpr = p.c() * dr - r.c() * dp;

                        RT distpq = diffapq * diffapq +
                            diffbpq * diffbpq +
                            diffcpq * diffcpq;

                        RT distpr = diffapr * diffapr +
                            diffbpr * diffbpr +
                            diffcpr * diffcpr;

                        return CGAL::is_positive(dq * dq *
                                                 distpr - dr * dr * distpq);
                    }

                private:
                    Point_3 origin;
            };

        // Less signed distance to plane
        template < typename K >
            struct Less_signed_distance_to_plane_3_dual_point
            {
                typedef typename K::RT         RT;
                typedef typename K::Plane_3    Plane_3;
                typedef typename K::Point_3   Point_3;
                typedef bool                   result_type;

                Less_signed_distance_to_plane_3_dual_point (Point_3 const& o =
                                                            Point_3(0, 0, 0)) :
                    origin(o)
                {}

                result_type
                    operator()(const Plane_dual<K> &p,
                               const Plane_3 &q,
                               const Plane_3 &r) const
                    {
                        Plane_3 p1 = p.p1;
                        Plane_3 p2 = p.p2;
                        Plane_3 p3 = p.p3;

                        RT dp1 = p1.d() + origin.x() * p1.a()
                            + origin.y() * p1.b() + origin.z() * p1.c();
                        RT dp2 = p2.d() + origin.x() * p2.a()
                            + origin.y() * p2.b() + origin.z() * p2.c();
                        RT dp3 = p3.d() + origin.x() * p3.a()
                            + origin.y() * p3.b() + origin.z() * p3.c();
                        RT dq = q.d() + origin.x() * q.a()
                            + origin.y() * q.b() + origin.z() * q.c();
                        RT dr = r.d() + origin.x() * r.a()
                            + origin.y() * r.b() + origin.z() * r.c();

                        // Compute the normal to the plane
                        RT alpha = (dp1 * p2.b() - dp2 * p1.b()) *
                            (dp1 * p3.c() - dp3 * p1.c()) -
                            (dp1 * p2.c() - dp2 * p1.c()) *
                            (dp1 * p3.b() - dp3 * p1.b());

                        RT beta  = (dp1 * p2.c() - dp2 * p1.c()) *
                            (dp1 * p3.a() - dp3 * p1.a()) -
                            (dp1 * p2.a() - dp2 * p1.a()) *
                            (dp1 * p3.c() - dp3 * p1.c());

                        RT gamma = (dp1 * p2.a() - dp2 * p1.a()) *
                            (dp1 * p3.b() - dp3 * p1.b()) -
                            (dp1 * p2.b() - dp2 * p1.b()) *
                            (dp1 * p3.a() - dp3 * p1.a());

                        RT sumpq = alpha * q.a() + beta * q.b() + gamma * q.c();
                        RT sumpr = alpha * r.a() + beta * r.b() + gamma * r.c();
                        RT diff = dr * sumpq - dq * sumpr;

                        if (CGAL::is_positive(dq * dr)) {
                            return CGAL::is_positive(diff);
                        } else {
                            return CGAL::is_negative(diff);
                        }
                    }

                private:
                    Point_3 origin;
            };

        // Oriented side
        template < typename K >
            struct Oriented_side_3_dual_point
            {
                typedef typename K::RT         RT;
                typedef typename K::Plane_3    Plane_3;
                typedef typename K::Point_3   Point_3;
                typedef CGAL::Oriented_side    result_type;

                Oriented_side_3_dual_point (Point_3 const& o =
                                                   Point_3(0, 0, 0)) :
                    origin(o)
                {}

                result_type
                    operator()(const Plane_dual<K> p, const Plane_3 &q) const
                    {
                        Plane_3 p1 = p.p1;
                        Plane_3 p2 = p.p2;
                        Plane_3 p3 = p.p3;

                        RT dp1 = p1.d() + origin.x() * p1.a()
                            + origin.y() * p1.b() + origin.z() * p1.c();
                        RT dp2 = p2.d() + origin.x() * p2.a()
                            + origin.y() * p2.b() + origin.z() * p2.c();
                        RT dp3 = p3.d() + origin.x() * p3.a()
                            + origin.y() * p3.b() + origin.z() * p3.c();
                        RT dq = q.d() + origin.x() * q.a()
                            + origin.y() * q.b() + origin.z() * q.c();

                        // Compute the normal to the plane
                        RT alpha = (dp1 * p2.b() - dp2 * p1.b()) *
                            (dp1 * p3.c() - dp3 * p1.c()) -
                            (dp1 * p2.c() - dp2 * p1.c()) *
                            (dp1 * p3.b() - dp3 * p1.b());

                        RT beta  = (dp1 * p2.c() - dp2 * p1.c()) *
                            (dp1 * p3.a() - dp3 * p1.a()) -
                            (dp1 * p2.a() - dp2 * p1.a()) *
                            (dp1 * p3.c() - dp3 * p1.c());

                        RT gamma = (dp1 * p2.a() - dp2 * p1.a()) *
                            (dp1 * p3.b() - dp3 * p1.b()) -
                            (dp1 * p2.b() - dp2 * p1.b()) *
                            (dp1 * p3.a() - dp3 * p1.a());

                        // Test if q is on the positive side of p
                        RT prod = alpha * (p1.a() * dq - q.a() * dp1) +
                            beta * (p1.b() * dq - q.b() * dp1) +
                            gamma * (p1.c() * dq - q.c() * dp1);

                        if (CGAL::is_positive(dp1 * dq)) {
                            if (CGAL::is_positive(prod)) {
                                return CGAL::ON_POSITIVE_SIDE;
                            } else if (CGAL::is_negative(prod)) {
                                return CGAL::ON_NEGATIVE_SIDE;
                            } else {
                                return CGAL::ON_ORIENTED_BOUNDARY;
                            }
                        } else {
                            if (CGAL::is_negative(prod)) {
                                return CGAL::ON_POSITIVE_SIDE;
                            } else if (CGAL::is_positive(prod)) {
                                return CGAL::ON_NEGATIVE_SIDE;
                            } else {
                                return CGAL::ON_ORIENTED_BOUNDARY;
                            }
                        }
                    }

                private:
                    Point_3 origin;
            };

    } // namespace Convex_hull_3
} // namespace CGAL

#endif // CGAL_CH3_DUAL_PREDICATES_H

