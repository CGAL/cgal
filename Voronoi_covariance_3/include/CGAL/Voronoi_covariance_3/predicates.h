#ifndef PREDICATS_H
#define PREDICATS_H

#include <CGAL/predicates/sign_of_determinant.h>

namespace CGAL
{
    namespace Voronoi_covariance_3
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
                typedef bool                  result_type;

                result_type
                    operator()(const Plane_3 &p, const Plane_3 &q) const
                    {
                        RT diffa = p.a() * q.d() - p.d() * q.a();
                        RT diffb = p.b() * q.d() - p.d() * q.b();
                        RT diffc = p.c() * q.d() - p.d() * q.c();

                        return CGAL_AND_3(CGAL::is_zero(diffa),
                                          CGAL::is_zero(diffb),
                                          CGAL::is_zero(diffc));
                    }
            };

        // Collinear
        template < typename K >
            struct Collinear_3_dual_point
            {
                typedef typename K::RT        RT;
                typedef typename K::Plane_3   Plane_3;
                typedef typename K::Vector_3  Vector_3;
                typedef bool                  result_type;

                result_type
                    operator()(const Plane_3 &p,
                               const Plane_3 &q,
                               const Plane_3 &r) const
                    {
                        RT diffapq = p.d() * q.a() - q.d() * p.a();
                        RT diffbpq = p.d() * q.b() - q.d() * p.b();
                        RT diffcpq = p.d() * q.c() - q.d() * p.c();

                        RT diffapr = p.d() * r.a() - r.d() * p.a();
                        RT diffbpr = p.d() * r.b() - r.d() * p.b();
                        RT diffcpr = p.d() * r.c() - r.d() * p.c();

                        // Cross product
                        RT cross1 = diffbpq * diffcpr - diffcpq * diffbpr;
                        RT cross2 = diffcpq * diffapr - diffapq * diffcpr;
                        RT cross3 = diffapq * diffbpr - diffbpq * diffapr;

                        return CGAL_AND_3(CGAL::is_zero(cross1),
                                          CGAL::is_zero(cross2),
                                          CGAL::is_zero(cross3));
                    }
            };

        // Coplanar
        template < typename K >
            struct Coplanar_3_dual_point
            {
                typedef typename K::RT        RT;
                typedef typename K::Plane_3   Plane_3;
                typedef bool                  result_type;

                result_type
                    operator()(const Plane_3 &p,
                               const Plane_3 &q,
                               const Plane_3 &r,
                               const Plane_3 &s) const
                    {
                        RT diffapq = p.d() * q.a() - q.d() * p.a();
                        RT diffbpq = p.d() * q.b() - q.d() * p.b();
                        RT diffcpq = p.d() * q.c() - q.d() * p.c();

                        RT diffapr = p.d() * r.a() - r.d() * p.a();
                        RT diffbpr = p.d() * r.b() - r.d() * p.b();
                        RT diffcpr = p.d() * r.c() - r.d() * p.c();

                        RT diffaps = p.d() * s.a() - s.d() * p.a();
                        RT diffbps = p.d() * s.b() - s.d() * p.b();
                        RT diffcps = p.d() * s.c() - s.d() * p.c();

                        return (CGAL::sign_of_determinant(diffapq, diffapr, diffaps,
                                                          diffbpq, diffbpr, diffbps,
                                                          diffcpq, diffcpr, diffcps)
                            == CGAL::ZERO);
                    }
            };

        // Has on positive side
        template < typename K >
            struct Has_on_positive_side_3_dual_point
            {
                typedef typename K::RT         RT;
                typedef typename K::Plane_3    Plane_3;
                typedef bool                   result_type;

                result_type
                    operator()(const Plane_dual<K> p, const Plane_3 &q) const
                    {
                        Plane_3 p1 = p.p1;
                        Plane_3 p2 = p.p2;
                        Plane_3 p3 = p.p3;

                        // Compute the normal to the plane
                        RT alpha = (p1.d() * p2.b() - p2.d() * p1.b()) *
                            (p1.d() * p3.c() - p3.d() * p1.c()) -
                            (p1.d() * p2.d() - p2.d() * p1.c()) *
                            (p1.d() * p3.b() - p3.d() * p1.b());

                        RT beta  = (p1.d() * p2.c() - p2.d() * p1.c()) *
                            (p1.d() * p3.a() - p3.d() * p1.a()) -
                            (p1.d() * p2.a() - p2.d() * p1.a()) *
                            (p1.d() * p3.c() - p3.d() * p1.c());

                        RT gamma = (p1.d() * p2.a() - p2.d() * p1.a()) *
                            (p1.d() * p3.b() - p3.d() * p1.b()) -
                            (p1.d() * p2.b() - p2.d() * p1.b()) *
                            (p1.d() * p3.a() - p3.d() * p1.a());

                        // Test if q is on the positive side of p
                        RT prod = alpha * (p1.a() * q.d() - q.a() * p1.d()) +
                            beta * (p1.b() * q.d() - q.b() * p1.d()) +
                            gamma * (p1.c() * q.d() - q.c() * p1.d());

                        if (CGAL::is_positive(p1.d() * q.d())) {
                            return CGAL::is_positive(prod);
                        } else {
                            return CGAL::is_negative(prod);
                        }
                    }
            };

        // Less distance to point
        template < typename K >
            struct Less_distance_to_point_3_dual_point
            {
                typedef typename K::RT        RT;
                typedef typename K::Plane_3   Plane_3;
                typedef bool                  result_type;

                result_type
                    operator()(const Plane_3 &p,
                               const Plane_3 &q,
                               const Plane_3 &r) const
                    {
                        RT diffapq = p.a() * q.d() - q.a() * p.d();
                        RT diffbpq = p.b() * q.d() - q.b() * p.d();
                        RT diffcpq = p.c() * q.d() - q.c() * p.d();

                        RT diffapr = p.a() * r.d() - r.a() * r.d();
                        RT diffbpr = p.b() * r.d() - r.b() * r.d();
                        RT diffcpr = p.c() * r.d() - r.c() * r.d();

                        RT distpq = diffapq * diffapq +
                            diffbpq * diffbpq +
                            diffcpq * diffcpq;

                        RT distpr = diffapr * diffapr +
                            diffbpr * diffbpr +
                            diffcpr * diffcpr;

                        return CGAL::is_positive(q.d() * q.d() *
                                                 distpr - r.d() * r.d() * distpq);
                    }
            };

        // Less signed distance to plane
        template < typename K >
            struct Less_signed_distance_to_plane_3_dual_point
            {
                typedef typename K::RT         RT;
                typedef typename K::Plane_3    Plane_3;
                typedef bool                   result_type;

                result_type
                    operator()(const Plane_dual<K> &p,
                               const Plane_3 &q,
                               const Plane_3 &r) const
                    {
                        Plane_3 p1 = p.p1;
                        Plane_3 p2 = p.p2;
                        Plane_3 p3 = p.p3;

                        // Compute the normal to the plane
                        RT alpha = (p1.d() * p2.b() - p2.d() * p1.b()) *
                            (p1.d() * p3.c() - p3.d() * p1.c()) -
                            (p1.d() * p2.d() - p2.d() * p1.c()) *
                            (p1.d() * p3.b() - p3.d() * p1.b());

                        RT beta  = (p1.d() * p2.c() - p2.d() * p1.c()) *
                            (p1.d() * p3.a() - p3.d() * p1.a()) -
                            (p1.d() * p2.a() - p2.d() * p1.a()) *
                            (p1.d() * p3.c() - p3.d() * p1.c());

                        RT gamma = (p1.d() * p2.a() - p2.d() * p1.a()) *
                            (p1.d() * p3.b() - p3.d() * p1.b()) -
                            (p1.d() * p2.b() - p2.d() * p1.b()) *
                            (p1.d() * p3.a() - p3.d() * p1.a());

                        RT sumpq = alpha * q.a() + beta * q.b() + gamma * q.c();
                        RT sumpr = alpha * r.a() + beta * r.b() + gamma * r.c();
                        RT diff = r.d() * sumpq - q.d() * sumpr;

                        if (CGAL::is_positive(q.d() * r.d())) {
                            return CGAL::is_positive(diff);
                        } else {
                            return CGAL::is_negative(diff);
                        }
                    }
            };

        // Oriented side
        template < typename K >
            struct Oriented_side_3_dual_point
            {
                typedef typename K::RT         RT;
                typedef typename K::Plane_3    Plane_3;
                typedef CGAL::Oriented_side    result_type;

                result_type
                    operator()(const Plane_dual<K> p, const Plane_3 &q) const
                    {
                        Plane_3 p1 = p.p1;
                        Plane_3 p2 = p.p2;
                        Plane_3 p3 = p.p3;

                        // Compute the normal to the plane
                        RT alpha = (p1.d() * p2.b() - p2.d() * p1.b()) *
                            (p1.d() * p3.c() - p3.d() * p1.c()) -
                            (p1.d() * p2.d() - p2.d() * p1.c()) *
                            (p1.d() * p3.b() - p3.d() * p1.b());

                        RT beta  = (p1.d() * p2.c() - p2.d() * p1.c()) *
                            (p1.d() * p3.a() - p3.d() * p1.a()) -
                            (p1.d() * p2.a() - p2.d() * p1.a()) *
                            (p1.d() * p3.c() - p3.d() * p1.c());

                        RT gamma = (p1.d() * p2.a() - p2.d() * p1.a()) *
                            (p1.d() * p3.b() - p3.d() * p1.b()) -
                            (p1.d() * p2.b() - p2.d() * p1.b()) *
                            (p1.d() * p3.a() - p3.d() * p1.a());

                        // Test if q is on the positive side of p
                        RT prod = alpha * (p1.a() * q.d() - q.a() * p1.d()) +
                            beta * (p1.b() * q.d() - q.b() * p1.d()) +
                            gamma * (p1.c() * q.d() - q.c() * p1.d());

                        if (CGAL::is_positive(p1.d() * q.d())) {
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
            };
    } // namespace Voronoi_covariance_3
} // namespace CGAL

#endif

