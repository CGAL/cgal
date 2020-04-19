// Copyright (c) 1999-2005
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra, Sylvain Pion, Michael Hoffmann

#ifndef CGAL_CARTESIAN_FUNCTION_OBJECTS_H
#define CGAL_CARTESIAN_FUNCTION_OBJECTS_H

#include <CGAL/Kernel/function_objects.h>
#include <CGAL/predicates/kernel_ftC2.h>
#include <CGAL/predicates/kernel_ftC3.h>
#include <CGAL/constructions/kernel_ftC2.h>
#include <CGAL/constructions/kernel_ftC3.h>
#include <CGAL/Cartesian/solve_3.h>

namespace CGAL {

namespace CartesianKernelFunctors {

  using namespace CommonKernelFunctors;

  template <typename K>
  class Angle_2
  {
    typedef typename K::Point_2  Point_2;
    typedef typename K::Vector_2 Vector_2;
  public:
    typedef typename K::Angle   result_type;

    result_type
    operator()(const Vector_2& u, const Vector_2& v) const
    { return angleC2(u.x(), u.y(), v.x(), v.y()); }

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return angleC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y()); }

    result_type
    operator()(const Point_2& p, const Point_2& q,
               const Point_2& r, const Point_2& s) const
    {
      return angleC2(p.x(), p.y(),
                     q.x(), q.y(),
                     r.x(), r.y(),
                     s.x(), s.y());
    }
  };

  template <typename K>
  class Angle_3
  {
    typedef typename K::Point_3  Point_3;
    typedef typename K::Vector_3 Vector_3;
  public:
    typedef typename K::Angle    result_type;

    result_type
    operator()(const Vector_3& u, const Vector_3& v) const
    {
      return angleC3(u.x(), u.y(), u.z(),
                     v.x(), v.y(), v.z());
    }
    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      return angleC3(p.x(), p.y(), p.z(),
                     q.x(), q.y(), q.z(),
                     r.x(), r.y(), r.z());
    }

    result_type
    operator()(const Point_3& p, const Point_3& q,
               const Point_3& r, const Point_3& s) const
    {
      return angleC3(p.x(), p.y(), p.z(),
                     q.x(), q.y(), q.z(),
                     r.x(), r.y(), r.z(),
                     s.x(), s.y(), s.z());
    }

    result_type
    operator()(const Point_3& p, const Point_3& q,
               const Point_3& r, const Vector_3& n) const
    {
      return enum_cast<Angle>(orientation(p,q,r,r+n));
    }
  };

  template <typename K>
  class Are_parallel_2
  {
    typedef typename K::Line_2          Line_2;
    typedef typename K::Segment_2       Segment_2;
    typedef typename K::Ray_2           Ray_2;

  public:
    typedef typename K::Boolean         result_type;

    result_type
    operator()(const Line_2& l1, const Line_2& l2) const
    { return parallelC2(l1.a(), l1.b(), l2.a(), l2.b()); }

    result_type
    operator()(const Segment_2& s1, const Segment_2& s2) const
    { return parallelC2(s1.source().x(), s1.source().y(),
                        s1.target().x(), s1.target().y(),
                        s2.source().x(), s2.source().y(),
                        s2.target().x(), s2.target().y());
    }

    result_type
    operator()(const Ray_2& r1, const Ray_2& r2) const
    { return parallelC2(r1.source().x(), r1.source().y(),
                        r1.second_point().x(), r1.second_point().y(),
                        r2.source().x(), r2.source().y(),
                        r2.second_point().x(), r2.second_point().y());
    }
  };

  template <typename K>
  class Are_parallel_3
  {
    typedef typename K::Line_3          Line_3;
    typedef typename K::Segment_3       Segment_3;
    typedef typename K::Ray_3           Ray_3;
    typedef typename K::Plane_3         Plane_3;

  public:
    typedef typename K::Boolean         result_type;

    result_type
    operator()(const Line_3& l1, const Line_3& l2) const
    { return parallelC3(
                l1.to_vector().x(), l1.to_vector().y(), l1.to_vector().z(),
                l2.to_vector().x(), l2.to_vector().y(), l2.to_vector().z());
    }

    result_type
    operator()(const Plane_3& h1, const Plane_3& h2) const
    { return parallelC3(h1.a(), h1.b(), h1.c(),
                        h2.a(), h2.b(), h2.c());
    }

    result_type
    operator()(const Segment_3& s1, const Segment_3& s2) const
    { return parallelC3(s1.source().x(), s1.source().y(), s1.source().z(),
                        s1.target().x(), s1.target().y(), s1.target().z(),
                        s2.source().x(), s2.source().y(), s2.source().z(),
                        s2.target().x(), s2.target().y(), s2.target().z());
    }

    result_type
    operator()(const Ray_3& r1, const Ray_3& r2) const
    { return parallelC3(r1.source().x(), r1.source().y(), r1.source().z(),
        r1.second_point().x(), r1.second_point().y(), r1.second_point().z(),
                        r2.source().x(), r2.source().y(), r2.source().z(),
        r2.second_point().x(), r2.second_point().y(), r2.second_point().z());
    }
  };

  template <typename K>
  class Bounded_side_2
  {
    typedef typename K::Point_2         Point_2;
    typedef typename K::Circle_2        Circle_2;
    typedef typename K::Triangle_2      Triangle_2;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;
  public:
    typedef typename K::Bounded_side    result_type;

    result_type
    operator()( const Circle_2& c, const Point_2& p) const
    {
      typename K::Compute_squared_distance_2 squared_distance;
      return enum_cast<Bounded_side>(CGAL::compare(c.squared_radius(),
                                                   squared_distance(c.center(),p)));
    }

    result_type
    operator()( const Triangle_2& t, const Point_2& p) const
    {
      typename K::Collinear_are_ordered_along_line_2
        collinear_are_ordered_along_line;
      typename K::Orientation_2 orientation;
      typename K::Orientation o1 = orientation(t.vertex(0), t.vertex(1), p),
                              o2 = orientation(t.vertex(1), t.vertex(2), p),
                              o3 = orientation(t.vertex(2), t.vertex(3), p);

      if (o2 == o1 && o3 == o1)
        return ON_BOUNDED_SIDE;
      return
        (o1 == COLLINEAR
         && collinear_are_ordered_along_line(t.vertex(0), p, t.vertex(1))) ||
        (o2 == COLLINEAR
         && collinear_are_ordered_along_line(t.vertex(1), p, t.vertex(2))) ||
        (o3 == COLLINEAR
         && collinear_are_ordered_along_line(t.vertex(2), p, t.vertex(3)))
        ? ON_BOUNDARY
        : ON_UNBOUNDED_SIDE;
    }

    result_type
    operator()( const Iso_rectangle_2& r, const Point_2& p) const
    {
      bool x_incr = (r.xmin() < p.x()) && (p.x() < r.xmax()),
           y_incr = (r.ymin() < p.y()) && (p.y() < r.ymax());
      if (x_incr)
        {
          if (y_incr)
            return ON_BOUNDED_SIDE;
          if ( (p.y() == r.ymin()) || (r.ymax() == p.y()) )
            return ON_BOUNDARY;
        }
      if ( (p.x() == r.xmin()) || (r.xmax() == p.x()) )
        if ( y_incr || (p.y() == r.ymin()) || (r.ymax() == p.y()) )
          return ON_BOUNDARY;

      return ON_UNBOUNDED_SIDE;
    }
  };

  template <typename K>
  class Bounded_side_3
  {
    typedef typename K::FT              FT;
    typedef typename K::Point_3         Point_3;
    typedef typename K::Sphere_3        Sphere_3;
    typedef typename K::Circle_3        Circle_3;
    typedef typename K::Tetrahedron_3   Tetrahedron_3;
    typedef typename K::Iso_cuboid_3    Iso_cuboid_3;
  public:
    typedef typename K::Bounded_side    result_type;

    result_type
    operator()( const Sphere_3& s, const Point_3& p) const
    { return s.rep().bounded_side(p); }

    result_type
    operator()( const Circle_3& s, const Point_3& p) const
    { return s.rep().bounded_side(p); }

    result_type
    operator()( const Tetrahedron_3& t, const Point_3& p) const
    {
      FT alpha, beta, gamma;

      Cartesian_internal::solve(t.vertex(1)-t.vertex(0),
                                t.vertex(2)-t.vertex(0),
                                t.vertex(3)-t.vertex(0),
                                p - t.vertex(0), alpha, beta, gamma);
      if (   (alpha < 0) || (beta < 0) || (gamma < 0)
          || (alpha + beta + gamma > 1) )
          return ON_UNBOUNDED_SIDE;

      if (   (alpha == 0) || (beta == 0) || (gamma == 0)
          || (alpha+beta+gamma == 1) )
        return ON_BOUNDARY;

      return ON_BOUNDED_SIDE;
    }

    result_type
    operator()( const Iso_cuboid_3& c, const Point_3& p) const
    {
      return c.rep().bounded_side(p);
    }

  };

  template <typename K>
  class Collinear_are_ordered_along_line_2
  {
    typedef typename K::Point_2         Point_2;
  public:
    typedef typename K::Boolean         result_type;

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      CGAL_kernel_exactness_precondition( collinear(p, q, r) );
      return collinear_are_ordered_along_lineC2
        (p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
    }
  };

  template <typename K>
  class Collinear_are_ordered_along_line_3
  {
    typedef typename K::Point_3         Point_3;
  public:
    typedef typename K::Boolean         result_type;

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      CGAL_kernel_exactness_precondition( collinear(p, q, r) );
      return collinear_are_ordered_along_lineC3(p.x(), p.y(), p.z(),
                                                q.x(), q.y(), q.z(),
                                                r.x(), r.y(), r.z());
    }
  };

  template <typename K>
  class Collinear_are_strictly_ordered_along_line_2
  {
    typedef typename K::Point_2         Point_2;
  public:
    typedef typename K::Boolean         result_type;

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      CGAL_kernel_exactness_precondition( collinear(p, q, r) );
      return collinear_are_strictly_ordered_along_lineC2
        (p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
    }
  };

  template <typename K>
  class Collinear_are_strictly_ordered_along_line_3
  {
    typedef typename K::Point_3         Point_3;
  public:
    typedef typename K::Boolean         result_type;

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      CGAL_kernel_exactness_precondition( collinear(p, q, r) );
      return collinear_are_strictly_ordered_along_lineC3(p.x(), p.y(), p.z(),
                                                         q.x(), q.y(), q.z(),
                                                         r.x(), r.y(), r.z());
    }
  };

  template <typename K>
  class Collinear_has_on_2
  {
    typedef typename K::Point_2               Point_2;
    typedef typename K::Ray_2                 Ray_2;
    typedef typename K::Segment_2             Segment_2;
  public:
    typedef typename K::Boolean               result_type;

    result_type
    operator()( const Ray_2& r, const Point_2& p) const
    {
      const Point_2 & source = r.source();
      const Point_2 & second = r.second_point();
      switch(make_certain(compare_x(source, second))) {
      case SMALLER:
        return compare_x(source, p) != LARGER;
      case LARGER:
        return compare_x(p, source) != LARGER;
      default:
        switch(make_certain(compare_y(source, second))){
        case SMALLER:
          return compare_y(source, p) != LARGER;
        case LARGER:
          return compare_y(p, source) != LARGER;
        default:
          return true; // p == source
        }
      } // switch
    }

    result_type
    operator()( const Segment_2& s, const Point_2& p) const
    {
      return collinear_are_ordered_along_line(s.source(), p, s.target());
    }
  };

  template <typename K>
  class Collinear_2
  {
    typedef typename K::Point_2        Point_2;
    typedef typename K::Orientation_2  Orientation_2;
    Orientation_2 o;
  public:
    typedef typename K::Boolean        result_type;

    Collinear_2() {}
    Collinear_2(const Orientation_2 o_) : o(o_) {}

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    { return o(p, q, r) == COLLINEAR; }
  };

  template <typename K>
  class Collinear_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef typename K::Boolean    result_type;

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      return collinearC3(p.x(), p.y(), p.z(),
                         q.x(), q.y(), q.z(),
                         r.x(), r.y(), r.z());
    }
  };

  template <typename K>
  class Compare_angle_with_x_axis_2
  {
    typedef typename K::Direction_2        Direction_2;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type
    operator()(const Direction_2& d1, const Direction_2& d2) const
    {
      return compare_angle_with_x_axisC2(d1.dx(), d1.dy(), d2.dx(), d2.dy());
    }
  };

  template <typename K>
  class Compare_distance_2
  {
    typedef typename K::Point_2            Point_2;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      return cmp_dist_to_pointC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
    }

    template <class T1, class T2, class T3>
    result_type
    operator()(const T1& p, const T2& q, const T3& r) const
    {
      return CGAL::compare(squared_distance(p, q), squared_distance(p, r));
    }

    template <class T1, class T2, class T3, class T4>
    result_type
    operator()(const T1& p, const T2& q, const T3& r, const T4& s) const
    {
      return CGAL::compare(squared_distance(p, q), squared_distance(r, s));
    }
  };

  template <typename K>
  class Compare_distance_3
  {
    typedef typename K::Point_3            Point_3;
    typedef typename K::Segment_3          Segment_3;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      return cmp_dist_to_pointC3(p.x(), p.y(), p.z(),
                                 q.x(), q.y(), q.z(),
                                 r.x(), r.y(), r.z());
    }

    result_type
    operator()(const Point_3& p1, const Segment_3& s1, const Segment_3& s2) const
    {
      return CGAL::internal::compare_distance_pssC3(p1,s1,s2, K());
    }

    result_type
    operator()(const Point_3& p1, const Point_3& p2, const Segment_3& s2) const
    {
      return CGAL::internal::compare_distance_ppsC3(p1,p2,s2, K());
    }

    result_type
    operator()(const Point_3& p1, const Segment_3& s2, const Point_3& p2) const
    {
      return opposite(CGAL::internal::compare_distance_ppsC3(p1,p2,s2, K()));
    }

    template <class T1, class T2, class T3>
    result_type
    operator()(const T1& p, const T2& q, const T3& r) const
    {
      return CGAL::compare(squared_distance(p, q), squared_distance(p, r));
    }

    template <class T1, class T2, class T3, class T4>
    result_type
    operator()(const T1& p, const T2& q, const T3& r, const T4& s) const
    {
      return CGAL::compare(squared_distance(p, q), squared_distance(r, s));
    }
  };

  template < typename K >
  class Compare_power_distance_2
  {
  public:
    typedef typename K::Weighted_point_2         Weighted_point_2;
    typedef typename K::Point_2                  Point_2;
    typedef typename K::Comparison_result        Comparison_result;

    typedef Comparison_result                    result_type;

    Comparison_result operator()(const Point_2& r,
                                 const Weighted_point_2& p,
                                 const Weighted_point_2& q) const
    {
      return CGAL::compare_power_distanceC2(p.x(), p.y(), p.weight(),
                                            q.x(), q.y(), q.weight(),
                                            r.x(), r.y());
    }
  };

  template <typename K>
  class Compare_signed_distance_to_line_2
  {
    typedef typename K::Point_2   Point_2;
    typedef typename K::Line_2    Line_2;
    typedef typename K::Equal_2   Equal_2;
    typedef typename K::Less_signed_distance_to_line_2 Less_signed_distance_to_line_2;

  public:
    typedef typename K::Comparison_result     result_type;

    result_type
    operator()(const Point_2& a, const Point_2& b,
               const Point_2& c, const Point_2& d) const
    {
      CGAL_kernel_precondition_code(Equal_2 equal;)
          CGAL_kernel_precondition(! equal(a,b));
      return cmp_signed_dist_to_lineC2( a.x(), a.y(),
                                        b.x(), b.y(),
                                        c.x(), c.y(),
                                        d.x(), d.y());
    }

    result_type
    operator()(const Line_2& l, const Point_2& p, const Point_2& q) const
    {
      Less_signed_distance_to_line_2 less = K().less_signed_distance_to_line_2_object();
      if (less(l, p, q)) return SMALLER;
      if (less(l, q, p)) return LARGER;
      return EQUAL;
    }
  };

  template <typename K>
  class Compare_squared_radius_3
  {
    typedef typename K::Point_3            Point_3;
    typedef typename K::FT                 FT;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r, const Point_3& s, const FT& ft) const
    {
      return CGAL::compare(squared_radiusC3(p.x(), p.y(), p.z(),
                                            q.x(), q.y(), q.z(),
                                            r.x(), r.y(), r.z(),
                                            s.x(), s.y(), s.z() ),
                           ft);
    }

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r, const FT& ft) const
    {
      return CGAL::compare(squared_radiusC3(p.x(), p.y(), p.z(),
                                            q.x(), q.y(), q.z(),
                                            r.x(), r.y(), r.z()),
                           ft);
    }

    result_type
    operator()(const Point_3& p, const Point_3& q, const FT& ft) const
    {
      return CGAL::compare(squared_radiusC3(p.x(), p.y(), p.z(),
                                            q.x(), q.y(), q.z() ),
                           ft);
    }

    result_type
    operator()(const Point_3&, const FT& ft) const
    {
      return - CGAL_NTS sign(ft);
    }
  };



  template <typename K>
  class Compare_slope_2
  {
    typedef typename K::Line_2             Line_2;
    typedef typename K::Segment_2          Segment_2;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type
    operator()(const Line_2& l1, const Line_2& l2) const
    {
      return compare_slopesC2(l1.a(), l1.b(), l2.a(), l2.b());
    }

    result_type
    operator()(const Segment_2& s1, const Segment_2& s2) const
    {
      return compare_slopesC2(s1.source().x(), s1.source().y(),
                              s1.target().x(), s1.target().y(),
                              s2.source().x(), s2.source().y(),
                              s2.target().x(), s2.target().y());
    }
  };

  template <typename K>
  class Compare_x_at_y_2
  {
    typedef typename K::Point_2             Point_2;
    typedef typename K::Line_2              Line_2;
  public:
    typedef typename K::Comparison_result   result_type;

    result_type
    operator()( const Point_2& p, const Line_2& h) const
    { return compare_y_at_xC2(p.y(), p.x(), h.b(), h.a(), h.c()); }

    result_type
    operator()( const Point_2& p, const Line_2& h1, const Line_2& h2) const
    {
      return compare_y_at_xC2(p.y(), h1.b(), h1.a(), h1.c(),
                              h2.b(), h2.a(), h2.c());
    }

    result_type
    operator()( const Line_2& l1, const Line_2& l2, const Line_2& h) const
    {
      return compare_y_at_xC2(l1.b(), l1.a(), l1.c(), l2.b(), l2.a(), l2.c(),
                              h.b(), h.a(), h.c());
    }

    result_type
    operator()( const Line_2& l1, const Line_2& l2,
                const Line_2& h1, const Line_2& h2) const
    {
      return compare_y_at_xC2(l1.b(), l1.a(), l1.c(), l2.b(), l2.a(), l2.c(),
                              h1.b(), h1.a(), h1.c(), h2.b(), h2.a(), h2.c());
    }
  };

  template <typename K>
  class Compare_xyz_3
  {
    typedef typename K::Point_3             Point_3;
  public:
    typedef typename K::Comparison_result   result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    {
      return compare_lexicographically_xyzC3(p.x(), p.y(), p.z(),
                                             q.x(), q.y(), q.z());
    }
  };

  template <typename K>
  class Compare_xy_2
  {
    typedef typename K::Point_2            Point_2;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    { return compare_lexicographically_xyC2(p.x(), p.y(), q.x(), q.y()); }
  };

  template <typename K>
  class Compare_xy_3
  {
    typedef typename K::Point_3            Point_3;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return compare_lexicographically_xyC2(p.x(), p.y(), q.x(), q.y()); }
  };

  template <typename K>
  class Compare_x_2
  {
    typedef typename K::Point_2             Point_2;
    typedef typename K::Line_2              Line_2;
  public:
    typedef typename K::Comparison_result   result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    { return CGAL::compare(p.x(), q.x()); }

    result_type
    operator()( const Point_2& p, const Line_2& l, const Line_2& h) const
    { return compare_xC2(p.x(), l.a(), l.b(), l.c(), h.a(), h.b(), h.c()); }

    result_type
    operator()( const Line_2& l, const Line_2& h1, const Line_2& h2) const
    {
      return compare_xC2(l.a(), l.b(), l.c(), h1.a(), h1.b(), h1.c(),
                         h2.a(), h2.b(), h2.c());
    }

    result_type
    operator()( const Line_2& l1, const Line_2& l2,
                const Line_2& h1, const Line_2& h2) const
    {
      return compare_xC2(l1.a(), l1.b(), l1.c(), l2.a(), l2.b(), l2.c(),
                         h1.a(), h1.b(), h1.c(), h2.a(), h2.b(), h2.c());
    }
  };

  template <typename K>
  class Compare_x_3
  {
    typedef typename K::Point_3             Point_3;
  public:
    typedef typename K::Comparison_result   result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return CGAL::compare(p.x(), q.x()); }
  };

  template <typename K>
  class Compare_yx_2
  {
    typedef typename K::Point_2            Point_2;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    { return compare_lexicographically_xyC2(p.y(), p.x(), q.y(), q.x()); }
  };

  template <typename K>
  class Compare_y_at_x_2
  {
    typedef typename K::Point_2             Point_2;
    typedef typename K::Line_2              Line_2;
    typedef typename K::Segment_2           Segment_2;
  public:
    typedef typename K::Comparison_result   result_type;

    result_type
    operator()( const Point_2& p, const Line_2& h) const
    { return compare_y_at_xC2(p.x(), p.y(), h.a(), h.b(), h.c()); }

    result_type
    operator()( const Point_2& p, const Line_2& h1, const Line_2& h2) const
    {
      return compare_y_at_xC2(p.x(), h1.a(), h1.b(), h1.c(),
                              h2.a(), h2.b(), h2.c());
    }

    result_type
    operator()( const Line_2& l1, const Line_2& l2, const Line_2& h) const
    {
      return compare_y_at_xC2(l1.a(), l1.b(), l1.c(), l2.a(), l2.b(), l2.c(),
                              h.a(), h.b(), h.c());
    }

    result_type
    operator()( const Line_2& l1, const Line_2& l2,
                const Line_2& h1, const Line_2& h2) const
    {
      return compare_y_at_xC2(l1.a(), l1.b(), l1.c(), l2.a(), l2.b(), l2.c(),
                              h1.a(), h1.b(), h1.c(), h2.a(), h2.b(), h2.c());
    }

    result_type
    operator()( const Point_2& p, const Segment_2& s) const
    {
      return compare_y_at_xC2(p.x(), p.y(),
                              s.source().x(), s.source().y(),
                              s.target().x(), s.target().y());
    }

    result_type
    operator()( const Point_2& p,
                const Segment_2& s1, const Segment_2& s2) const
    {
      return compare_y_at_x_segment_C2(p.x(),
                                       s1.source().x(), s1.source().y(),
                                       s1.target().x(), s1.target().y(),
                                       s2.source().x(), s2.source().y(),
                                       s2.target().x(), s2.target().y());
    }
  };

  template <typename K>
  class Compare_y_2
  {
    typedef typename K::Point_2             Point_2;
    typedef typename K::Line_2              Line_2;
  public:
    typedef typename K::Comparison_result   result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    { return CGAL::compare(p.y(), q.y()); }

    result_type
    operator()( const Point_2& p, const Line_2& l1, const Line_2& l2) const
    {
      return compare_xC2(p.y(),
                         l1.b(), l1.a(), l1.c(),
                         l2.b(), l2.a(), l2.c());
    }

    result_type
    operator()( const Line_2& l, const Line_2& h1, const Line_2& h2) const
    {
      return compare_xC2(l.b(), l.a(), l.c(), h1.b(), h1.a(), h1.c(),
                         l.b(), l.a(), l.c(), h2.b(), h2.a(), h2.c());
    }

    result_type
    operator()( const Line_2& l1, const Line_2& l2,
                const Line_2& h1, const Line_2& h2) const
    {
      return compare_xC2(l1.b(), l1.a(), l1.c(), l2.b(), l2.a(), l2.c(),
                         h1.b(), h1.a(), h1.c(), h2.b(), h2.a(), h2.c());
    }
  };

  template <typename K>
  class Compare_y_3
  {
    typedef typename K::Point_3             Point_3;
  public:
    typedef typename K::Comparison_result   result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return CGAL::compare(p.y(), q.y()); }
  };

  template <typename K>
  class Compare_z_3
  {
    typedef typename K::Point_3            Point_3;
  public:
    typedef typename K::Comparison_result  result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return CGAL::compare(p.z(), q.z()); }
  };

  template <class K>
  class Compute_approximate_area_3
  {
    typedef typename K::Circle_3                  Circle_3;
    typedef typename K::FT                        FT;

  public:

    typedef double result_type;

    result_type
    operator() (const Circle_3 & c) const
    // { return c.rep().approximate_area(); }
    { return CGAL_PI * to_double(c.squared_radius()); }
  };

  template <class K>
  class Compute_approximate_squared_length_3
  {
    typedef typename K::Circle_3                  Circle_3;
    typedef typename K::FT                        FT;

  public:

    typedef double result_type;

    result_type
    operator() (const Circle_3 & c) const
    // { return c.rep().approximate_squared_length(); }
    { return CGAL_PI * CGAL_PI * 4.0 * to_double(c.squared_radius()); }
  };


  template <typename K>
  class Compute_area_2
  {
    typedef typename K::FT                FT;
    typedef typename K::Iso_rectangle_2   Iso_rectangle_2;
    typedef typename K::Triangle_2        Triangle_2;
    typedef typename K::Point_2           Point_2;
  public:
    typedef FT               result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q, const Point_2& r ) const
    {
      FT v1x = q.x() - p.x();
      FT v1y = q.y() - p.y();
      FT v2x = r.x() - p.x();
      FT v2y = r.y() - p.y();
      return determinant(v1x, v1y, v2x, v2y)/2;
    }

    result_type
    operator()( const Iso_rectangle_2& r ) const
    { return (r.xmax()-r.xmin()) * (r.ymax()-r.ymin()); }

    result_type
    operator()( const Triangle_2& t ) const
    { return t.area(); }
  };

  template <typename K>
  class Compute_area_divided_by_pi_3
  {
    typedef typename K::Circle_3                  Circle_3;
    typedef typename K::FT                        FT;

  public:

    typedef FT result_type;

    result_type
    operator()(const Circle_3 & c) const
    { return c.rep().area_divided_by_pi(); }

  };

  template <typename K>
  class Compute_determinant_2
  {
    typedef typename K::FT                FT;
    typedef typename K::Vector_2          Vector_2;
  public:
    typedef FT               result_type;

    result_type
    operator()(const Vector_2& v, const Vector_2& w) const
    {
        return determinant(v.x(), v.y(), w.x(), w.y());
    }
  };

  template <typename K>
  class Compute_determinant_3
  {
    typedef typename K::FT                FT;
    typedef typename K::Vector_3          Vector_3;
  public:
    typedef FT               result_type;

    result_type
    operator()(const Vector_3& v, const Vector_3& w, const Vector_3& t) const
    {
        return determinant(v.x(), v.y(), v.z(),
                                 w.x(), w.y(), w.z(),
                                 t.x(), t.y(), t.z());
    }
  };

  template <typename K>
  class Compute_scalar_product_2
  {
    typedef typename K::FT                FT;
    typedef typename K::Vector_2          Vector_2;
  public:
    typedef FT               result_type;

    result_type
    operator()(const Vector_2& v, const Vector_2& w) const
    {
        return v.x() * w.x() + v.y() * w.y();
    }
  };

  template <typename K>
  class Compute_scalar_product_3
  {
    typedef typename K::FT                FT;
    typedef typename K::Vector_3          Vector_3;
  public:
    typedef FT               result_type;

    result_type
    operator()(const Vector_3& v, const Vector_3& w) const
    {
        return v.x() * w.x() + v.y() * w.y() + v.z() * w.z();
    }
  };

  template <typename K>
  class Compute_squared_area_3
  {
    typedef typename K::FT                FT;
    typedef typename K::Point_3           Point_3;
    typedef typename K::Triangle_3        Triangle_3;
  public:
    typedef FT               result_type;

    result_type
    operator()( const Triangle_3& t ) const
    {
        return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }

    result_type
    operator()( const Point_3& p, const Point_3& q, const Point_3& r ) const
    {
        return squared_areaC3(p.x(), p.y(), p.z(),
                              q.x(), q.y(), q.z(),
                              r.x(), r.y(), r.z());
    }
  };

  // FIXME
  template <typename K>
  class Compute_squared_distance_Point_Point_2
  {
    typedef typename K::FT       FT;
    typedef typename K::Point_2  Point_2;
  public:
    typedef FT               result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    {
      return squared_distanceC2(p.x(), p.y(), q.x(), q.y());
    }
  };

  template <class K>
  class Compute_squared_length_divided_by_pi_square_3
  {
    typedef typename K::Circle_3                  Circle_3;
    typedef typename K::FT                        FT;

  public:

    typedef FT result_type;

    result_type
    operator() (const Circle_3 & c) const
    { return c.rep().squared_length_divided_by_pi_square(); }

  };

  template <typename K>
  class Compute_squared_radius_2
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Circle_2    Circle_2;
  public:
    typedef FT                      result_type;

    result_type
    operator()( const Circle_2& c) const
    { return c.rep().squared_radius(); }

    result_type
    operator()( const Point_2& /*p*/) const
    { return FT(0); }

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    { return squared_radiusC2(p.x(), p.y(), q.x(), q.y()); }

    result_type
    operator()( const Point_2& p, const Point_2& q, const Point_2& r) const
    { return squared_radiusC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y()); }
  };

} //namespace CartesianKernelFunctors

// For the non specialized template will do the right thing,
// namely return a copy of an FT

namespace CartesianKernelFunctors {

  template <typename K>
  class Compute_squared_radius_3
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_3     Point_3;
    typedef typename K::Sphere_3    Sphere_3;
    typedef typename K::Circle_3    Circle_3;
  public:
    typedef FT               result_type;

    result_type
    operator()( const Sphere_3& s) const
    { return s.rep().squared_radius(); }

    result_type
    operator()( const Circle_3& c) const
    { return c.rep().squared_radius(); }

    result_type
    operator()( const Point_3& /*p*/) const
    { return FT(0); }

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    {
      return squared_radiusC3(p.x(), p.y(), p.z(),
                              q.x(), q.y(), q.z());
    }

    result_type
    operator()( const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      return squared_radiusC3(p.x(), p.y(), p.z(),
                              q.x(), q.y(), q.z(),
                              r.x(), r.y(), r.z());
    }

    result_type
    operator()( const Point_3& p, const Point_3& q,
                const Point_3& r, const Point_3& s) const
    {
      return squared_radiusC3(p.x(), p.y(), p.z(),
                              q.x(), q.y(), q.z(),
                              r.x(), r.y(), r.z(),
                              s.x(), s.y(), s.z());
    }
  };

  template <typename K>
  class Compute_volume_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
    typedef typename K::Iso_cuboid_3   Iso_cuboid_3;
  public:
    typedef FT               result_type;

    result_type
    operator()(const Point_3& p0, const Point_3& p1,
               const Point_3& p2, const Point_3& p3) const
    {
      return determinant<FT>(p1.x()-p0.x(), p1.y()-p0.y(), p1.z()-p0.z(),
                             p2.x()-p0.x(), p2.y()-p0.y(), p2.z()-p0.z(),
                             p3.x()-p0.x(), p3.y()-p0.y(), p3.z()-p0.z())/6;
    }

    result_type
    operator()( const Tetrahedron_3& t ) const
    {
      return this->operator()(t.vertex(0), t.vertex(1),
                              t.vertex(2), t.vertex(3));
    }

    result_type
    operator()( const Iso_cuboid_3& c ) const
    { return c.rep().volume(); }
  };


  template <typename K>
  class Compute_x_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_2        Point_2;
    typedef typename K::Vector_2       Vector_2;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_2& p) const
    {
      return p.rep().x();
    }

    result_type
    operator()(const Vector_2& v) const
    {
      return v.rep().x();
    }
  };

  template <typename K>
  class Compute_x_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Vector_3       Vector_3;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_3& p) const
    {
      return p.rep().x();
    }

    result_type
    operator()(const Vector_3& v) const
    {
      return v.rep().x();
    }
  };


  template <typename K>
  class Compute_y_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_2        Point_2;
    typedef typename K::Vector_2       Vector_2;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_2& p) const
    {
      return p.rep().y();
    }

    result_type
    operator()(const Vector_2& v) const
    {
      return v.rep().y();
    }
  };


  template <typename K>
  class Compute_y_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Vector_3       Vector_3;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_3& p) const
    {
      return p.rep().y();
    }

    result_type
    operator()(const Vector_3& v) const
    {
      return v.rep().y();
    }
  };

  template <typename K>
  class Compute_z_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Vector_3       Vector_3;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_3& p) const
    {
      return p.rep().z();
    }

    result_type
    operator()(const Vector_3& v) const
    {
      return v.rep().z();
    }
  };



  template <typename K>
  class Compute_dx_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Direction_2    Direction_2;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Direction_2& d) const
    {
      return d.rep().dx();
    }
  };

  template <typename K>
  class Compute_dx_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Direction_3    Direction_3;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Direction_3& d) const
    {
      return d.rep().dx();
    }
  };

  template <typename K>
  class Compute_dy_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Direction_2    Direction_2;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Direction_2& d) const
    {
      return d.rep().dy();
    }
  };

  template <typename K>
  class Compute_dy_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Direction_3    Direction_3;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Direction_3& d) const
    {
      return d.rep().dy();
    }
  };

  template <typename K>
  class Compute_dz_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Direction_3        Direction_3;

  public:
    typedef const FT&               result_type;

    result_type
    operator()(const Direction_3& d) const
    {
      return d.rep().dz();
    }
  };

  template <typename K>
  class Compute_hx_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_2        Point_2;
    typedef typename K::Vector_2       Vector_2;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_2& p) const
    {
      return p.rep().hx();
    }

    result_type
    operator()(const Vector_2& v) const
    {
      return v.rep().hx();
    }
  };

  template <typename K>
  class Compute_hx_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Vector_3       Vector_3;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_3& p) const
    {
      return p.rep().hx();
    }

    result_type
    operator()(const Vector_3& v) const
    {
      return v.rep().hx();
    }
  };

  template <typename K>
  class Compute_hy_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_2        Point_2;
    typedef typename K::Vector_2       Vector_2;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_2& p) const
    {
      return p.rep().hy();
    }

    result_type
    operator()(const Vector_2& v) const
    {
      return v.rep().hy();
    }
  };

  template <typename K>
  class Compute_hy_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Vector_3       Vector_3;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_3& p) const
    {
      return p.rep().hy();
    }

    result_type
    operator()(const Vector_3& v) const
    {
      return v.rep().hy();
    }
  };

  template <typename K>
  class Compute_hz_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Vector_3       Vector_3;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_3& p) const
    {
      return p.rep().hz();
    }

    result_type
    operator()(const Vector_3& v) const
    {
      return v.rep().hz();
    }
  };

  template <typename K>
  class Compute_hw_2
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_2        Point_2;
    typedef typename K::Vector_2       Vector_2;

  public:
    typedef const FT&                  result_type;

    result_type
    operator()(const Point_2& p) const
    {
      return p.rep().hw();
    }

    result_type
    operator()(const Vector_2& v) const
    {
      return v.rep().hw();
    }
  };

  template <typename K>
  class Compute_hw_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Vector_3       Vector_3;

  public:
    typedef const FT&                 result_type;

    result_type
    operator()(const Point_3& p) const
    {
      return p.rep().hw();
    }

    result_type
    operator()(const Vector_3& v) const
    {
      return v.rep().hw();
    }
  };


  template <typename K>
  class Compute_xmin_2
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;

  public:
    typedef const FT&                   result_type;

    result_type
    operator()(const Iso_rectangle_2& r) const
    {
      return (r.min)().x();
    }
  };

  template <typename K>
  class Compute_xmax_2
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;

  public:
    typedef const FT&                   result_type;

    result_type
    operator()(const Iso_rectangle_2& r) const
    {
      return (r.max)().x();
    }
  };

  template <typename K>
  class Compute_ymin_2
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;

  public:
    typedef const FT&                   result_type;

    result_type
    operator()(const Iso_rectangle_2& r) const
    {
      return (r.min)().y();
    }
  };

  template <typename K>
  class Compute_ymax_2
  {
    typedef typename K::FT              FT;
    typedef typename K::Iso_rectangle_2 Iso_rectangle_2;

  public:
    typedef const FT&                   result_type;

    result_type
    operator()(const Iso_rectangle_2& r) const
    {
      return (r.max)().y();
    }
  };


  template <typename K>
  class Construct_barycenter_2
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
  public:
    typedef Point_2                 result_type;

    result_type
    operator()(const Point_2& p1, const FT&w1, const Point_2& p2) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      barycenterC2(p1.x(), p1.y(), w1, p2.x(), p2.y(), x, y);
      return construct_point_2(x, y);
    }

    result_type
    operator()(const Point_2& p1, const FT& w1, const Point_2& p2, const FT& w2) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      barycenterC2(p1.x(), p1.y(), w1, p2.x(), p2.y(), w2, x, y);
      return construct_point_2(x, y);
    }

    result_type
    operator()(const Point_2& p1, const FT& w1, const Point_2& p2, const FT& w2,
               const Point_2& p3) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      barycenterC2(p1.x(), p1.y(), w1, p2.x(), p2.y(), w2, p3.x(), p3.y(), x, y);
      return construct_point_2(x, y);
    }

    result_type
    operator()(const Point_2& p1, const FT& w1, const Point_2& p2, const FT& w2,
               const Point_2& p3, const FT& w3) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      barycenterC2(p1.x(), p1.y(), w1, p2.x(), p2.y(), w2, p3.x(), p3.y(), w3, x, y);
      return construct_point_2(x, y);
    }

    result_type
    operator()(const Point_2& p1, const FT& w1, const Point_2& p2, const FT& w2,
               const Point_2& p3, const FT& w3, const Point_2& p4) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      barycenterC2(p1.x(), p1.y(), w1, p2.x(), p2.y(), w2, p3.x(), p3.y(), w3, p4.x(), p4.y(), x, y);
      return construct_point_2(x, y);
    }

    result_type
    operator()(const Point_2& p1, const FT& w1, const Point_2& p2, const FT& w2,
               const Point_2& p3, const FT& w3, const Point_2& p4, const FT& w4) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      barycenterC2(p1.x(), p1.y(), w1, p2.x(), p2.y(), w2, p3.x(), p3.y(), w3, p4.x(), p4.y(), w4, x, y);
      return construct_point_2(x, y);
    }

  };

  template <typename K>
  class Construct_barycenter_3
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_3     Point_3;
  public:
    typedef Point_3                 result_type;

    result_type
    operator()(const Point_3& p1, const FT&w1, const Point_3& p2) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      barycenterC3(p1.x(), p1.y(), p1.z(), w1, p2.x(), p2.y(), p2.z(), x, y, z);
      return construct_point_3(x, y, z);
    }

    result_type
    operator()(const Point_3& p1, const FT& w1, const Point_3& p2, const FT& w2) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      barycenterC3(p1.x(), p1.y(), p1.z(), w1, p2.x(), p2.y(), p2.z(), w2, x, y, z);
      return construct_point_3(x, y, z);
    }

    result_type
    operator()(const Point_3& p1, const FT& w1, const Point_3& p2, const FT& w2,
               const Point_3& p3) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      barycenterC3(p1.x(), p1.y(), p1.z(), w1, p2.x(), p2.y(), p2.z(), w2, p3.x(), p3.y(), p3.z(), x, y, z);
      return construct_point_3(x, y, z);
    }

    result_type
    operator()(const Point_3& p1, const FT& w1, const Point_3& p2, const FT& w2,
               const Point_3& p3, const FT& w3) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      barycenterC3(p1.x(), p1.y(), p1.z(), w1, p2.x(), p2.y(), p2.z(), w2,
                   p3.x(), p3.y(), p3.z(), w3, x, y, z);
      return construct_point_3(x, y, z);
    }

    result_type
    operator()(const Point_3& p1, const FT& w1, const Point_3& p2, const FT& w2,
               const Point_3& p3, const FT& w3, const Point_3& p4) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      barycenterC3(p1.x(), p1.y(), p1.z(), w1, p2.x(), p2.y(), p2.z(), w2,
                   p3.x(), p3.y(), p3.z(), w3, p4.x(), p4.y(), p4.z(), x, y, z);
      return construct_point_3(x, y, z);
    }

    result_type
    operator()(const Point_3& p1, const FT& w1, const Point_3& p2, const FT& w2,
               const Point_3& p3, const FT& w3, const Point_3& p4, const FT& w4) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      barycenterC3(p1.x(), p1.y(), p1.z(), w1, p2.x(), p2.y(), p2.z(), w2,
                   p3.x(), p3.y(), p3.z(), w3, p4.x(), p4.y(), p4.z(), w4, x, y, z);
      return construct_point_3(x, y, z);
    }

  };

  template <typename K>
  class Construct_base_vector_3
  {
    typedef typename K::Vector_3   Vector_3;
    typedef typename K::Plane_3    Plane_3;
    typedef typename K::FT         FT;
    typedef typename K::Construct_cross_product_vector_3
    Construct_cross_product_vector_3;
    typedef typename K::Construct_orthogonal_vector_3
    Construct_orthogonal_vector_3;
    Construct_cross_product_vector_3 cp;
    Construct_orthogonal_vector_3 co;
  public:
    typedef Vector_3         result_type;

    Construct_base_vector_3() {}
    Construct_base_vector_3(const Construct_cross_product_vector_3& cp_,
                            const Construct_orthogonal_vector_3& co_)
      : cp(cp_), co(co_)
    {}

    result_type
    operator()( const Plane_3& h, int index ) const
    {
      if (index == 1) {
        if ( CGAL_NTS is_zero(h.a()) )  // parallel to x-axis
          return Vector_3(FT(1), FT(0), FT(0));

        if ( CGAL_NTS is_zero(h.b()) )  // parallel to y-axis
          return Vector_3(FT(0), FT(1), FT(0));

        if ( CGAL_NTS is_zero(h.c()) )  // parallel to z-axis
          return Vector_3(FT(0), FT(0), FT(1));

        FT a = CGAL::abs(h.a()),
          b = CGAL::abs(h.b()),
          c = CGAL::abs(h.c());

        // to avoid badly defined vectors with coordinates all close
        // to 0 when the plane is almost horizontal, we ignore the
        // smallest coordinate instead of always ignoring Z
        if (a <= b && a <= c)
          return Vector_3(FT(0), -h.c(), h.b());

        if (b <= a && b <= c)
          return Vector_3(-h.c(), FT(0), h.a());

        return Vector_3(-h.b(), h.a(), FT(0));
      } else {
        return cp(co(h), this->operator()(h,1));
      }
    }
  };


  template <typename K>
  class Construct_bbox_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Segment_2        Segment_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Triangle_2       Triangle_2;
    typedef typename K::Circle_2         Circle_2;
  public:
    typedef Bbox_2                       result_type;

    result_type
    operator()(const Point_2& p) const
    {
      std::pair<double,double> xp = CGAL_NTS to_interval(p.x());
      std::pair<double,double> yp = CGAL_NTS to_interval(p.y());
      return Bbox_2(xp.first, yp.first, xp.second, yp.second);
    }

    result_type
    operator()(const Segment_2& s) const
    { return s.source().bbox() + s.target().bbox(); }

    result_type
    operator()(const Triangle_2& t) const
    {
      Bbox_2 bb = this->operator()(t.vertex(0));
      bb += this->operator()(t.vertex(1));
      bb += this->operator()(t.vertex(2));
      return bb;
      /*
          Microsoft (R) C/C++ Optimizing Compiler Version 18.00.40629.0 for x64
          produces a segfault of this functor for Simple_cartesian<Interval_nt<0>>
          with the original version of the code below
          Note that it also worked for 18.00.21005.1

      typename K::Construct_bbox_2 construct_bbox_2;
      return construct_bbox_2(t.vertex(0))
           + construct_bbox_2(t.vertex(1))
           + construct_bbox_2(t.vertex(2));
      */
    }

    result_type
    operator()(const Iso_rectangle_2& r) const
    {
      typename K::Construct_bbox_2 construct_bbox_2;
      return construct_bbox_2((r.min)()) + construct_bbox_2((r.max)());
    }

    result_type
    operator()(const Circle_2& c) const
    {
      typename K::Construct_bbox_2 construct_bbox_2;
      Bbox_2 b = construct_bbox_2(c.center());

      Interval_nt<> x (b.xmin(), b.xmax());
      Interval_nt<> y (b.ymin(), b.ymax());

      Interval_nt<> sqr = CGAL_NTS to_interval(c.squared_radius());
      Interval_nt<> r = CGAL::sqrt(sqr);
      Interval_nt<> minx = x-r;
      Interval_nt<> maxx = x+r;
      Interval_nt<> miny = y-r;
      Interval_nt<> maxy = y+r;

      return Bbox_2(minx.inf(), miny.inf(), maxx.sup(), maxy.sup());
    }
  };


  template <typename K>
  class Construct_bbox_3
  {
    typedef typename K::Point_3          Point_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Iso_cuboid_3     Iso_cuboid_3;
    typedef typename K::Triangle_3       Triangle_3;
    typedef typename K::Tetrahedron_3    Tetrahedron_3;
    typedef typename K::Sphere_3         Sphere_3;
    typedef typename K::Circle_3         Circle_3;
  public:
    typedef Bbox_3          result_type;

    Bbox_3
    operator()(const Point_3& p) const
    {
      std::pair<double,double> xp = CGAL_NTS to_interval(p.x());
      std::pair<double,double> yp = CGAL_NTS to_interval(p.y());
      std::pair<double,double> zp = CGAL_NTS to_interval(p.z());
      return Bbox_3(xp.first, yp.first, zp.first,
                    xp.second, yp.second, zp.second);
    }

    Bbox_3
    operator()(const Segment_3& s) const
    { return s.source().bbox() + s.target().bbox(); }

    Bbox_3
    operator()(const Triangle_3& t) const
    {
      typename K::Construct_bbox_3 construct_bbox_3;
      return construct_bbox_3(t.vertex(0))
           + construct_bbox_3(t.vertex(1))
           + construct_bbox_3(t.vertex(2));
    }

    Bbox_3
    operator()(const Iso_cuboid_3& r) const
    {
      typename K::Construct_bbox_3 construct_bbox_3;
      return construct_bbox_3((r.min)()) + construct_bbox_3((r.max)());
    }

    Bbox_3
    operator()(const Tetrahedron_3& t) const
    {
      typename K::Construct_bbox_3 construct_bbox_3;
      return construct_bbox_3(t.vertex(0)) + construct_bbox_3(t.vertex(1))
           + construct_bbox_3(t.vertex(2)) + construct_bbox_3(t.vertex(3));
    }

    Bbox_3
    operator()(const Sphere_3& s) const
    {
      typename K::Construct_bbox_3 construct_bbox_3;
      Bbox_3 b = construct_bbox_3(s.center());

      Interval_nt<> x (b.xmin(), b.xmax());
      Interval_nt<> y (b.ymin(), b.ymax());
      Interval_nt<> z (b.zmin(), b.zmax());

      Interval_nt<> sqr = CGAL_NTS to_interval(s.squared_radius());
      Interval_nt<> r = CGAL::sqrt(sqr);
      Interval_nt<> minx = x-r;
      Interval_nt<> maxx = x+r;
      Interval_nt<> miny = y-r;
      Interval_nt<> maxy = y+r;
      Interval_nt<> minz = z-r;
      Interval_nt<> maxz = z+r;

      return Bbox_3(minx.inf(), miny.inf(), minz.inf(),
                    maxx.sup(), maxy.sup(), maxz.sup());
    }

    Bbox_3
    operator()(const Circle_3& c) const
    { return c.rep().bbox(); }

  };


  template <typename K>
  class Construct_bisector_2
  {
    typedef typename K::FT      FT;
    typedef typename K::Point_2 Point_2;
    typedef typename K::Line_2  Line_2;
  public:
    typedef Line_2              result_type;

    result_type
    operator()(const Point_2& p, const Point_2& q) const
    {
      FT a, b, c;
      bisector_of_pointsC2(p.x(), p.y(), q.x(), q.y(), a, b, c);
      return Line_2(a, b, c);
    }

    result_type
    operator()(const Line_2& p, const Line_2& q) const
    {
      FT a, b, c;
      bisector_of_linesC2(p.a(), p.b(), p.c(),
                          q.a(), q.b(), q.c(),
                          a, b, c);
      return Line_2(a, b, c);
    }
  };

  template <typename K>
  class Construct_bisector_3
  {
    typedef typename K::FT        FT;
    typedef typename K::Point_3   Point_3;
    typedef typename K::Plane_3   Plane_3;
  public:
    typedef Plane_3               result_type;

    result_type
    operator()(const Point_3& p, const Point_3& q) const
    {
      FT a, b, c, d;
      bisector_of_pointsC3(p.x(), p.y(), p.z(),
                           q.x(), q.y(), q.z(),
                           a, b, c, d);
      return Plane_3(a, b, c, d);
    }

    result_type
    operator()(const Plane_3& p, const Plane_3& q) const
    {
      FT a, b, c, d;
      bisector_of_planesC3(p.a(), p.b(), p.c(), p.d(),
                           q.a(), q.b(), q.c(), q.d(),
                           a, b, c, d);
      return Plane_3(a, b, c, d);
    }
  };

  template <typename K>
  class Construct_centroid_2
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_2     Point_2;
    typedef typename K::Triangle_2  Triangle_2;
  public:
    typedef Point_2                 result_type;

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      centroidC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(), x, y);
      return construct_point_2(x, y);
    }

    result_type
    operator()(const Triangle_2& t) const
    {
      return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }

    result_type
    operator()(const Point_2& p, const Point_2& q,
               const Point_2& r, const Point_2& s) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      centroidC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(), s.x(), s.y(), x, y);
      return construct_point_2(x, y);
    }
  };

  template <typename K>
  class Construct_centroid_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Point_3        Point_3;
    typedef typename K::Triangle_3     Triangle_3;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
  public:
    typedef Point_3                    result_type;

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      centroidC3(p.x(), p.y(), p.z(),
                 q.x(), q.y(), q.z(),
                 r.x(), r.y(), r.z(),
                 x, y, z);
      return construct_point_3(x, y, z);
    }

    result_type
    operator()(const Point_3& p, const Point_3& q,
               const Point_3& r, const Point_3& s) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      centroidC3(p.x(), p.y(), p.z(),
                 q.x(), q.y(), q.z(),
                 r.x(), r.y(), r.z(),
                 s.x(), s.y(), s.z(),
                 x, y, z);
      return construct_point_3(x, y, z);
    }

    result_type
    operator()(const Triangle_3& t) const
    {
      return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }

    result_type
    operator()(const Tetrahedron_3& t) const
    {
      return this->operator()(t.vertex(0), t.vertex(1),
                              t.vertex(2), t.vertex(3));
    }
  };

  template <typename K>
  class Construct_circumcenter_2
  {
    typedef typename K::Point_2     Point_2;
    typedef typename K::Triangle_2  Triangle_2;
  public:
    typedef Point_2                 result_type;

    Point_2
    operator()(const Point_2& p, const Point_2& q) const
    {
      typename K::Construct_midpoint_2 construct_midpoint_2;
      return construct_midpoint_2(p, q);
    }

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      typename K::Construct_point_2 construct_point_2;
      typedef typename K::FT        FT;
      FT x, y;
      circumcenterC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(), x, y);
      return construct_point_2(x, y);
    }

    result_type
    operator()(const Triangle_2& t) const
    {
      return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }
  };

  template <typename K>
  class Construct_circumcenter_3
  {
    typedef typename K::FT             FT;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
    typedef typename K::Triangle_3     Triangle_3;
    typedef typename K::Point_3        Point_3;
  public:
    typedef Point_3                    result_type;

    Point_3
    operator()(const Point_3& p, const Point_3& q) const
    {
      typename K::Construct_midpoint_3 construct_midpoint_3;
      return construct_midpoint_3(p, q);
    }

    Point_3
    operator()(const Point_3& p, const Point_3& q, const Point_3& s) const
    {
      typename K::Construct_point_3 construct_point_3;
      // Translate s to origin to simplify the expression.
      FT psx = p.x()-s.x();
      FT psy = p.y()-s.y();
      FT psz = p.z()-s.z();
      FT ps2 = CGAL_NTS square(psx) + CGAL_NTS square(psy) + CGAL_NTS square(psz);
      FT qsx = q.x()-s.x();
      FT qsy = q.y()-s.y();
      FT qsz = q.z()-s.z();
      FT qs2 = CGAL_NTS square(qsx) + CGAL_NTS square(qsy) + CGAL_NTS square(qsz);
      FT rsx = psy*qsz-psz*qsy;
      FT rsy = psz*qsx-psx*qsz;
      FT rsz = psx*qsy-psy*qsx;

      // The following determinants can be developped and simplified.
      //
      // FT num_x = determinant(psy,psz,ps2,
      //                              qsy,qsz,qs2,
      //                              rsy,rsz,0);
      // FT num_y = determinant(psx,psz,ps2,
      //                              qsx,qsz,qs2,
      //                              rsx,rsz,0);
      // FT num_z = determinant(psx,psy,ps2,
      //                              qsx,qsy,qs2,
      //                              rsx,rsy,0);

      FT num_x = ps2 * determinant(qsy,qsz,rsy,rsz)
               - qs2 * determinant(psy,psz,rsy,rsz);
      FT num_y = ps2 * determinant(qsx,qsz,rsx,rsz)
               - qs2 * determinant(psx,psz,rsx,rsz);
      FT num_z = ps2 * determinant(qsx,qsy,rsx,rsy)
               - qs2 * determinant(psx,psy,rsx,rsy);

      FT den   = determinant(psx,psy,psz,
                             qsx,qsy,qsz,
                             rsx,rsy,rsz);

      CGAL_kernel_assertion( den != 0 );
      FT inv = 1 / (2 * den);

      FT x = s.x() + num_x*inv;
      FT y = s.y() - num_y*inv;
      FT z = s.z() + num_z*inv;
      return construct_point_3(x, y, z);
    }

    Point_3
    operator()(const Triangle_3& t) const
    {
      return this->operator()(t.vertex(0), t.vertex(1), t.vertex(2));
    }

    Point_3
    operator()(const Point_3& p, const Point_3& q,
               const Point_3& r, const Point_3& s) const
    {
      typename K::Construct_point_3 construct_point_3;
      // Translate p to origin to simplify the expression.
      FT qpx = q.x()-p.x();
      FT qpy = q.y()-p.y();
      FT qpz = q.z()-p.z();
      FT qp2 = CGAL_NTS square(qpx) + CGAL_NTS square(qpy) + CGAL_NTS square(qpz);
      FT rpx = r.x()-p.x();
      FT rpy = r.y()-p.y();
      FT rpz = r.z()-p.z();
      FT rp2 = CGAL_NTS square(rpx) + CGAL_NTS square(rpy) + CGAL_NTS square(rpz);
      FT spx = s.x()-p.x();
      FT spy = s.y()-p.y();
      FT spz = s.z()-p.z();
      FT sp2 = CGAL_NTS square(spx) + CGAL_NTS square(spy) + CGAL_NTS square(spz);

      FT num_x = determinant(qpy,qpz,qp2,
                             rpy,rpz,rp2,
                             spy,spz,sp2);
      FT num_y = determinant(qpx,qpz,qp2,
                             rpx,rpz,rp2,
                             spx,spz,sp2);
      FT num_z = determinant(qpx,qpy,qp2,
                             rpx,rpy,rp2,
                             spx,spy,sp2);
      FT den   = determinant(qpx,qpy,qpz,
                             rpx,rpy,rpz,
                             spx,spy,spz);
      CGAL_kernel_assertion( ! CGAL_NTS is_zero(den) );
      FT inv = 1 / (2 * den);

      FT x = p.x() + num_x*inv;
      FT y = p.y() - num_y*inv;
      FT z = p.z() + num_z*inv;
      return construct_point_3(x, y, z);
    }

    Point_3
    operator()(const Tetrahedron_3& t) const
    {
      return this->operator()(t.vertex(0), t.vertex(1),
                              t.vertex(2), t.vertex(3));
    }
  };




  template <typename K>
  class Construct_cross_product_vector_3
  {
    typedef typename K::Vector_3  Vector_3;
  public:
    typedef Vector_3              result_type;

    Vector_3
    operator()(const Vector_3& v, const Vector_3& w) const
    {
      return Vector_3(v.y() * w.z() - v.z() * w.y(),
                      v.z() * w.x() - v.x() * w.z(),
                      v.x() * w.y() - v.y() * w.x());
    }
  };

  template <typename K>
  class Construct_lifted_point_3
  {
    typedef typename K::Point_2                    Point_2;
    typedef typename K::Point_3                    Point_3;
    typedef typename K::Plane_3                    Plane_3;
    typedef typename K::Construct_base_vector_3    Construct_base_vector_3;
    typedef typename K::Construct_point_on_3       Construct_point_on_3;
    typedef typename K::Construct_scaled_vector_3  Construct_scaled_vector_3;
    typedef typename K::Construct_translated_point_3
    Construct_translated_point_3;
    Construct_base_vector_3 cb;
    Construct_point_on_3 cp;
    Construct_scaled_vector_3 cs;
    Construct_translated_point_3 ct;
  public:
    typedef Point_3          result_type;

    Construct_lifted_point_3() {}
    Construct_lifted_point_3(const Construct_base_vector_3& cb_,
                             const Construct_point_on_3& cp_,
                             const Construct_scaled_vector_3& cs_,
                             const Construct_translated_point_3& ct_)
      : cb(cb_), cp(cp_), cs(cs_), ct(ct_)
    {}

    Point_3
    operator()(const Plane_3& h, const Point_2& p) const
    {
      return ct(ct(cp(h), cs(cb(h,1), p.x())), cs(cb(h,2), p.y()));
    }
  };

  template <typename K>
  class Construct_direction_2
  {
    typedef typename K::Direction_2     Direction_2;
    typedef typename Direction_2::Rep   Rep;
    typedef typename K::Point_2         Point_2;
    typedef typename K::Vector_2        Vector_2;
    typedef typename K::Line_2          Line_2;
    typedef typename K::Ray_2           Ray_2;
    typedef typename K::Segment_2       Segment_2;
    typedef typename K::RT              RT;

  public:
    typedef Direction_2                 result_type;

    Rep // Direction_2
    operator()(Return_base_tag, const RT& x, const RT& y) const
    { return Rep(x, y); }

    Rep // Direction_2
    operator()(Return_base_tag, const Vector_2& v) const
    {
      return Rep(v.x(),v.y()); }

    Rep // Direction_2
    operator()(Return_base_tag, const Line_2& l) const
    { return Rep(l.b(), -l.a()); }

    Rep // Direction_2
    operator()(Return_base_tag, const Point_2& p, const Point_2& q) const
    {
      return Rep(q.x() - p.x(), q.y() - p.y());
    }

    Rep // Direction_2
    operator()(Return_base_tag, const Ray_2& r) const
    {
      return this->operator()(Return_base_tag(), r.source(), r.second_point());
    }

    Rep // Direction_2
    operator()(Return_base_tag, const Segment_2& s) const
    {
      return this->operator()(Return_base_tag(), s.source(), s.target());
    }


    Direction_2
    operator()(const RT& x, const RT& y) const
    { return this->operator()(Return_base_tag(), x, y); }

    Direction_2
    operator()(const Vector_2& v) const
    {
      return this->operator()(Return_base_tag(), v); }

    Direction_2
    operator()(const Line_2& l) const
    { return this->operator()(Return_base_tag(), l); }

    Direction_2
    operator()(const Point_2& p, const Point_2& q) const
    {
      return this->operator()(Return_base_tag(), p, q);
    }

    Direction_2
    operator()(const Ray_2& r) const
    {
      return this->operator()(Return_base_tag(), r);
    }

    Direction_2
    operator()(const Segment_2& s) const
    {
      return this->operator()(Return_base_tag(), s);
    }
  };

  template <typename K>
  class Construct_direction_3
  {
    typedef typename K::Direction_3     Direction_3;
    typedef typename K::Vector_3        Vector_3;
    typedef typename K::Line_3          Line_3;
    typedef typename K::Ray_3           Ray_3;
    typedef typename K::Segment_3       Segment_3;
    typedef typename K::RT              RT;
    typedef typename Direction_3::Rep   Rep;
  public:
    typedef Direction_3       result_type;

    Rep // Direction_3
    operator()(Return_base_tag, const RT& x, const RT& y, const RT& z) const
    { return Rep(x, y, z); }

    Rep // Direction_3
    operator()(Return_base_tag, const Vector_3& v) const
    { return Rep(v); }

    Rep // Direction_3
    operator()(Return_base_tag, const Line_3& l) const
    { return Rep(l); }

    Rep // Direction_3
    operator()(Return_base_tag, const Ray_3& r) const
    { return Rep(r); }

    Rep // Direction_3
    operator()(Return_base_tag, const Segment_3& s) const
    { return Rep(s); }


    Direction_3
    operator()(const RT& x, const RT& y, const RT& z) const
    { return this->operator()(Return_base_tag(), x, y, z); }

    Direction_3
    operator()(const Vector_3& v) const
    { return this->operator()(Return_base_tag(), v); }

    Direction_3
    operator()(const Line_3& l) const
    { return this->operator()(Return_base_tag(), l); }

    Direction_3
    operator()(const Ray_3& r) const
    { return this->operator()(Return_base_tag(), r); }

    Direction_3
    operator()(const Segment_3& s) const
    { return this->operator()(Return_base_tag(), s); }
  };

  template <typename K>
  class Construct_equidistant_line_3
  {
    typedef typename K::FT          FT;
    typedef typename K::Point_3     Point_3;
    typedef typename K::Vector_3    Vector_3;
    typedef typename K::Line_3      Line_3;
    typedef typename Line_3::Rep    Rep;
  public:
    typedef Line_3           result_type;

    Line_3
    operator()( const Point_3& p, const Point_3& q, const Point_3& s) const
    {
      CGAL_kernel_precondition(! collinear(p, q, s));

      // Translate s to origin to simplify the expression.
      FT psx = p.x()-s.x();
      FT psy = p.y()-s.y();
      FT psz = p.z()-s.z();
      FT ps2 = CGAL_NTS square(psx) + CGAL_NTS square(psy) + CGAL_NTS square(psz);
      FT qsx = q.x()-s.x();
      FT qsy = q.y()-s.y();
      FT qsz = q.z()-s.z();
      FT qs2 = CGAL_NTS square(qsx) + CGAL_NTS square(qsy) + CGAL_NTS square(qsz);
      FT rsx = psy*qsz-psz*qsy;
      FT rsy = psz*qsx-psx*qsz;
      FT rsz = psx*qsy-psy*qsx;

      // The following determinants can be developped and simplified.
      //
      // FT num_x = determinant(psy,psz,ps2,
      //                              qsy,qsz,qs2,
      //                              rsy,rsz,0);
      // FT num_y = determinant(psx,psz,ps2,
      //                              qsx,qsz,qs2,
      //                              rsx,rsz,0);
      // FT num_z = determinant(psx,psy,ps2,
      //                              qsx,qsy,qs2,
      //                              rsx,rsy,0);

      FT num_x = ps2 * determinant(qsy,qsz,rsy,rsz)
               - qs2 * determinant(psy,psz,rsy,rsz);
      FT num_y = ps2 * determinant(qsx,qsz,rsx,rsz)
               - qs2 * determinant(psx,psz,rsx,rsz);
      FT num_z = ps2 * determinant(qsx,qsy,rsx,rsy)
               - qs2 * determinant(psx,psy,rsx,rsy);

      FT den   = determinant(psx,psy,psz,
                             qsx,qsy,qsz,
                             rsx,rsy,rsz);

      CGAL_kernel_assertion( den != 0 );
      FT inv = 1 / (2 * den);

      FT x = s.x() + num_x*inv;
      FT y = s.y() - num_y*inv;
      FT z = s.z() + num_z*inv;
      return Rep(Point_3(x, y, z), Vector_3(rsx, rsy, rsz));
    }

  };

  template <typename K>
  class Construct_iso_rectangle_2
  {
    typedef typename K::RT               RT;
    typedef typename K::FT               FT;
    typedef typename K::Point_2          Point_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename Iso_rectangle_2::Rep     Rep;

  public:
    typedef Iso_rectangle_2              result_type;

    Rep // Iso_rectangle_2
    operator()(Return_base_tag, const Point_2& p, const Point_2& q, int) const
    {
      // I have to remove the assertions, because of Cartesian_converter.
      // CGAL_kernel_assertion(p.x()<=q.x());
      // CGAL_kernel_assertion(p.y()<=q.y());
      return Rep(p, q, 0);
    }

    Rep // Iso_rectangle_2
    operator()(Return_base_tag, const Point_2& p, const Point_2& q) const
    {
      FT minx, maxx, miny, maxy;
      if (p.x() < q.x()) { minx = p.x(); maxx = q.x(); }
      else               { minx = q.x(); maxx = p.x(); }
      if (p.y() < q.y()) { miny = p.y(); maxy = q.y(); }
      else               { miny = q.y(); maxy = p.y(); }

      return Rep(Point_2(minx, miny),
                 Point_2(maxx, maxy), 0);
    }

    Rep // Iso_rectangle_2
    operator()(Return_base_tag, const Point_2 &left,   const Point_2 &right,
               const Point_2 &bottom, const Point_2 &top) const
    {
      CGAL_kernel_assertion_code(typename K::Less_x_2 less_x;)
      CGAL_kernel_assertion_code(typename K::Less_y_2 less_y;)
      CGAL_kernel_assertion(!less_x(right, left));
      CGAL_kernel_assertion(!less_y(top, bottom));
      return Rep(Point_2(left.x(), bottom.y()),
                 Point_2(right.x(), top.y()), 0);
    }

    Rep // Iso_rectangle_2
    operator()(Return_base_tag, const RT& min_hx, const RT& min_hy,
               const RT& max_hx, const RT& max_hy) const
    {
      CGAL_kernel_precondition(min_hx <= max_hx);
      CGAL_kernel_precondition(min_hy <= max_hy);
      return Rep(Point_2(min_hx, min_hy),
                 Point_2(max_hx, max_hy), 0);
    }

    Rep // Iso_rectangle_2
    operator()(Return_base_tag, const RT& min_hx, const RT& min_hy,
               const RT& max_hx, const RT& max_hy, const RT& hw) const
    {
      if (hw == 1)
        return Rep(Point_2(min_hx, min_hy),
                   Point_2(max_hx, max_hy), 0);
      return Rep(Point_2(min_hx/hw, min_hy/hw),
                 Point_2(max_hx/hw, max_hy/hw), 0);
    }


    Iso_rectangle_2
    operator()(const Point_2& p, const Point_2& q, int i) const
    {
      return this->operator()(Return_base_tag(), p, q, i);
    }

    Iso_rectangle_2
    operator()(const Point_2& p, const Point_2& q) const
    {
      return this->operator()(Return_base_tag(), p, q);
    }

    Iso_rectangle_2
    operator()(const Point_2 &left,   const Point_2 &right,
               const Point_2 &bottom, const Point_2 &top) const
    {
      return this->operator()(Return_base_tag(), left, right, bottom, top);
    }

    Iso_rectangle_2
    operator()(const RT& min_hx, const RT& min_hy,
               const RT& max_hx, const RT& max_hy) const
    {
      return this->operator()(Return_base_tag(), min_hx, min_hy, max_hx, max_hy);
    }

    Iso_rectangle_2
    operator()(const RT& min_hx, const RT& min_hy,
               const RT& max_hx, const RT& max_hy, const RT& hw) const
    {
      return this->operator()(Return_base_tag(), min_hx, min_hy, max_hx, max_hy, hw);
    }
  };

  template <typename K>
  class Construct_line_2
  {
    typedef typename K::RT                        RT;
    typedef typename K::FT                        FT;
    typedef typename K::Point_2                   Point_2;
    typedef typename K::Direction_2               Direction_2;
    typedef typename K::Vector_2                  Vector_2;
    typedef typename K::Segment_2                 Segment_2;
    typedef typename K::Ray_2                     Ray_2;
    typedef typename K::Line_2                    Line_2;
    typedef typename Line_2::Rep                  Rep;
    typedef typename K::Construct_point_on_2      Construct_point_on_2;
    Construct_point_on_2 c;
  public:
    typedef Line_2            result_type;

    Construct_line_2() {}
    Construct_line_2(const Construct_point_on_2& c_) : c(c_) {}

    Rep // Line_2
    operator()(Return_base_tag, const RT& a, const RT& b, const RT& cc) const
    { return Rep(a, b, cc); }

    Rep // Line_2
    operator()(Return_base_tag, const Point_2& p, const Point_2& q) const
    {
      FT a, b, cc;
      line_from_pointsC2(p.x(), p.y(), q.x(), q.y(), a, b, cc);
      return Rep(a, b, cc);
    }

    Rep // Line_2
    operator()(Return_base_tag, const Point_2& p, const Direction_2& d) const
    {
      FT a, b, cc;
      line_from_point_directionC2(p.x(), p.y(), d.dx(), d.dy(), a, b, cc);
      return Rep(a, b, cc);
    }

    Rep // Line_2
    operator()(Return_base_tag, const Point_2& p, const Vector_2& v) const
    {
      FT a, b, cc;
      line_from_point_directionC2(p.x(), p.y(), v.x(), v.y(), a, b, cc);
      return Rep(a, b, cc);
    }

    Rep // Line_2
    operator()(Return_base_tag, const Segment_2& s) const
    { return this->operator()(Return_base_tag(), c(s, 0), c(s, 1)); }

    Rep // Line_2
    operator()(Return_base_tag, const Ray_2& r) const
    { return this->operator()(Return_base_tag(), c(r, 0), c(r, 1)); }


    Line_2
    operator()(const RT& a, const RT& b, const RT& cc) const
    { return this->operator()(Return_base_tag(), a, b, cc); }

    Line_2
    operator()(const Point_2& p, const Point_2& q) const
    { return this->operator()(Return_base_tag(), p, q); }

    Line_2
    operator()(const Point_2& p, const Direction_2& d) const
    { return this->operator()(Return_base_tag(), p, d); }

    Line_2
    operator()(const Point_2& p, const Vector_2& v) const
    { return this->operator()(Return_base_tag(), p, v); }

    Line_2
    operator()(const Segment_2& s) const
    { return this->operator()(Return_base_tag(), s); }

    Line_2
    operator()(const Ray_2& r) const
    { return this->operator()(Return_base_tag(), r); }
  };

  template <typename K>
  class Construct_line_3
  {
    typedef typename K::Point_3                   Point_3;
    typedef typename K::Direction_3               Direction_3;
    typedef typename K::Segment_3                 Segment_3;
    typedef typename K::Ray_3                     Ray_3;
    typedef typename K::Line_3                    Line_3;
    typedef typename K::Vector_3                  Vector_3;
    typedef typename Line_3::Rep                  Rep;
  public:
    typedef Line_3            result_type;

    Rep // Line_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q) const
    { return Rep(p, Vector_3(p, q)); }

    Rep // Line_3
    operator()(Return_base_tag, const Point_3& p, const Direction_3& d) const
    { return operator()(Return_base_tag(), p, Vector_3(d.dx(), d.dy(), d.dz())); }

    Rep // Line_3
    operator()(Return_base_tag, const Point_3& p, const Vector_3& v) const
    { return Rep(p, v); }

    Rep // Line_3
    operator()(Return_base_tag, const Segment_3& s) const
    { return Rep(s.source(), Vector_3(s.source(), s.target())); }

    Rep // Line_3
    operator()(Return_base_tag, const Ray_3& r) const
    { return Rep(r.source(), Vector_3(r.source(), r.second_point())); }


    Line_3
    operator()(const Point_3& p, const Point_3& q) const
    { return this->operator()(Return_base_tag(), p, q); }

    Line_3
    operator()(const Point_3& p, const Direction_3& d) const
    { return this->operator()(Return_base_tag(), p, d); }

    Line_3
    operator()(const Point_3& p, const Vector_3& v) const
    { return this->operator()(Return_base_tag(), p, v); }

    Line_3
    operator()(const Segment_3& s) const
    { return this->operator()(Return_base_tag(), s); }

    Line_3
    operator()(const Ray_3& r) const
    { return this->operator()(Return_base_tag(), r); }
  };

  template <typename K>
  class Construct_midpoint_2
  {
    typedef typename K::FT        FT;
    typedef typename K::Point_2   Point_2;
  public:
    typedef Point_2          result_type;

    Point_2
    operator()(const Point_2& p, const Point_2& q) const
    {
      typename K::Construct_point_2 construct_point_2;
      FT x, y;
      midpointC2(p.x(), p.y(), q.x(), q.y(), x, y);
      return construct_point_2(x, y);
    }
  };

  template <typename K>
  class Construct_midpoint_3
  {
    typedef typename K::FT        FT;
    typedef typename K::Point_3   Point_3;
  public:
    typedef Point_3               result_type;

    Point_3
    operator()(const Point_3& p, const Point_3& q) const
    {
      typename K::Construct_point_3 construct_point_3;
      FT x, y, z;
      midpointC3(p.x(), p.y(), p.z(), q.x(), q.y(), q.z(), x, y, z);
      return construct_point_3(x, y, z);
    }
  };

  template <typename K>
  class Construct_opposite_vector_2
  {
    typedef typename K::Vector_2    Vector_2;
  public:
    typedef Vector_2                result_type;

    Vector_2
    operator()( const Vector_2& v) const
    { return Vector_2(-v.x(), -v.y()); }
  };

  template <typename K>
  class Construct_difference_of_vectors_2
  {
    typedef typename K::Vector_2    Vector_2;
  public:
    typedef Vector_2                result_type;

    Vector_2
    operator()( const Vector_2& v, const Vector_2& w) const
    { return Vector_2(v.x()-w.x(), v.y()-w.y()); }
  };

  template <typename K>
  class Construct_difference_of_vectors_3
  {
    typedef typename K::Vector_3    Vector_3;
  public:
    typedef Vector_3                result_type;

    Vector_3
    operator()( const Vector_3& v, const Vector_3& w) const
    { return Vector_3(v.x()-w.x(), v.y()-w.y(), v.z()-w.z()); }
  };


  template < typename K >
  class Construct_radical_axis_2
  {
  public:
    typedef typename K::Weighted_point_2                Weighted_point_2;
    typedef typename K::Line_2                          Line_2;

    typedef Line_2           result_type;

    Line_2
    operator() ( const Weighted_point_2 & p, const Weighted_point_2 & q) const
    {
      typedef typename K::RT RT;
      RT a,b,c;
      radical_axisC2(p.x(),p.y(),p.weight(),q.x(),q.y(),q.weight(),a,b,c);
      return Line_2(a,b,c);
    }
  };




  template <typename K>
  class Construct_sum_of_vectors_2
  {
    typedef typename K::Vector_2    Vector_2;
  public:
    typedef Vector_2                result_type;

    Vector_2
    operator()( const Vector_2& v, const Vector_2& w) const
    { return Vector_2(v.x()+w.x(), v.y()+w.y()); }
  };

  template <typename K>
  class Construct_sum_of_vectors_3
  {
    typedef typename K::Vector_3    Vector_3;
  public:
    typedef Vector_3                result_type;

    Vector_3
    operator()( const Vector_3& v, const Vector_3& w) const
    { return Vector_3(v.x()+w.x(), v.y()+w.y(), v.z()+w.z()); }
  };

  template <typename K>
  class Construct_opposite_vector_3
  {
    typedef typename K::Vector_3    Vector_3;
  public:
    typedef Vector_3                result_type;

    Vector_3
    operator()( const Vector_3& v) const
    { return Vector_3(-v.x(), -v.y(), -v.z()); }
  };

  template <typename K>
  class Construct_orthogonal_vector_3
  {
    typedef typename K::FT FT;
    typedef typename K::Point_3     Point_3;
    typedef typename K::Vector_3    Vector_3;
    typedef typename K::Plane_3     Plane_3;
  public:
    typedef Vector_3                result_type;

    Vector_3
    operator()( const Plane_3& p ) const
    { return Vector_3(p.a(), p.b(), p.c()); }

    Vector_3
    operator()( const Point_3& p, const Point_3& q, const Point_3& r ) const
    {
      FT rpx = p.x()-r.x();
      FT rpy = p.y()-r.y();
      FT rpz = p.z()-r.z();
      FT rqx = q.x()-r.x();
      FT rqy = q.y()-r.y();
      FT rqz = q.z()-r.z();
      // Cross product rp * rq
      FT vx = rpy*rqz - rqy*rpz;
      FT vy = rpz*rqx - rqz*rpx;
      FT vz = rpx*rqy - rqx*rpy;
      typename K::Construct_vector_3 construct_vector;

      return construct_vector(vx, vy, vz);
    }
  };

  template <typename K>
  class Construct_perpendicular_vector_2
  {
    typedef typename K::Vector_2   Vector_2;
  public:
    typedef Vector_2               result_type;

    Vector_2
    operator()( const Vector_2& v, Orientation o) const
    {
      CGAL_kernel_precondition( o != COLLINEAR );
      if (o == COUNTERCLOCKWISE)
        return K().construct_vector_2_object()(-v.y(), v.x());
      else
        return K().construct_vector_2_object()(v.y(), -v.x());
    }
  };

  template <typename K>
  class Construct_perpendicular_direction_2
  {
    typedef typename K::Direction_2   Direction_2;
  public:
    typedef Direction_2               result_type;

    Direction_2
    operator()( const Direction_2& d, Orientation o) const
    {
      CGAL_kernel_precondition( o != COLLINEAR );
      if (o == COUNTERCLOCKWISE)
        return K().construct_direction_2_object()(-d.dy(), d.dx());
      else
        return K().construct_direction_2_object()(d.dy(), -d.dx());
    }
  };


  template <typename K>
  class Construct_perpendicular_line_2
  {
    typedef typename K::Line_2    Line_2;
    typedef typename K::Point_2   Point_2;
  public:
    typedef Line_2                result_type;

    Line_2
    operator()( const Line_2& l, const Point_2& p) const
    {
      typename K::FT fta, ftb, ftc;
      perpendicular_through_pointC2(l.a(), l.b(), p.x(), p.y(), fta, ftb, ftc);
      return Line_2(fta, ftb, ftc);
    }
  };


  template <typename K>
  class Construct_point_2
  {
    typedef typename K::RT         RT;
    typedef typename K::FT         FT;
    typedef typename K::Point_2    Point_2;
    typedef typename K::Weighted_point_2 Weighted_point_2;
    typedef typename K::Line_2     Line_2;
    typedef typename Point_2::Rep  Rep;
  public:

    template<typename>
    struct result {
      typedef Point_2 type;
    };

    template<typename F>
    struct result<F(Weighted_point_2)> {
      typedef const Point_2& type;
    };

    template<typename F>
    struct result<F(Point_2)> {
      typedef const Point_2& type;
    };

    Rep // Point_2
    operator()(Return_base_tag, Origin o) const
    { return Rep(o); }

    Rep // Point_2
    operator()(Return_base_tag, const RT& x, const RT& y) const
    { return Rep(x, y); }

    Rep // Point_2
    operator()(Return_base_tag, const RT& x, const RT& y, const RT& w) const
    { return Rep(x, y, w); }

    Point_2
    operator()(const Line_2& l) const
    {
      typename K::Construct_point_2 construct_point_2;
      typename K::FT x, y;
      line_get_pointC2(l.a(), l.b(), l.c(), FT(0), x, y);
      return construct_point_2(x,y);
    }

    Point_2
    operator()(const Line_2& l, const FT i) const
    {
      typename K::Construct_point_2 construct_point_2;
      typename K::FT x, y;
      line_get_pointC2(l.a(), l.b(), l.c(), i, x, y);
      return construct_point_2(x,y);
    }

    const Point_2&
    operator()(const Point_2 & p) const
    { return p; }

    const Point_2&
    operator()(const Weighted_point_2 & p) const
    { return p.rep().point(); }

    Point_2
    operator()(Origin o) const
    { return Point_2(o); }

    Point_2
    operator()(const RT& x, const RT& y) const
    { return Point_2(x, y); }

    Point_2
    operator()(const RT& x, const RT& y, const RT& w) const
    { return Point_2(x, y, w); }
  };

  template <typename K>
  class Construct_point_3
  {
    typedef typename K::RT               RT;
    typedef typename K::Point_3          Point_3;
    typedef typename K::Weighted_point_3 Weighted_point_3;
    typedef typename Point_3::Rep        Rep;

  public:

    template<typename>
    struct result {
      typedef Point_3 type;
    };

    template<typename F>
    struct result<F(Weighted_point_3)> {
      typedef const Point_3& type;
    };

    template<typename F>
    struct result<F(Point_3)> {
      typedef const Point_3& type;
    };


    Rep // Point_3
    operator()(Return_base_tag, Origin o) const
    { return Rep(o); }

    Rep // Point_3
    operator()(Return_base_tag, const RT& x, const RT& y, const RT& z) const
    { return Rep(x, y, z); }

    Rep // Point_3
    operator()(Return_base_tag, const RT& x, const RT& y, const RT& z, const RT& w) const
    { return Rep(x, y, z, w); }

    const Point_3&
    operator()(const Point_3 & p) const
    { return p; }

    const Point_3&
    operator()(const Weighted_point_3 & p) const
    { return p.rep().point(); }

    Point_3
    operator()(Origin o) const
    { return Point_3(o); }

    Point_3
    operator()(const RT& x, const RT& y, const RT& z) const
    { return Point_3(x, y, z); }

    Point_3
    operator()(const RT& x, const RT& y, const RT& z, const RT& w) const
    { return Point_3(x, y, z, w); }
  };


  template <typename K>
  class Construct_weighted_point_2
  {
    typedef typename K::FT                  FT;
    typedef typename K::Point_2             Point_2;
    typedef typename K::Weighted_point_2    Weighted_point_2;
    typedef typename Weighted_point_2::Rep  Rep;
  public:
    typedef Weighted_point_2                result_type;

    Rep
    operator()(Return_base_tag, Origin o) const
    { return Rep(o); }

    Rep
    operator()(Return_base_tag, const Point_2& p, const FT& w) const
    { return Rep(p,w); }

    Rep
    operator()(Return_base_tag, const FT& x, const FT& y) const
    { return Rep(x,y); }

    Weighted_point_2
    operator()(Origin o) const
    { return Weighted_point_2(o); }

    Weighted_point_2
    operator()(const Point_2& p, const FT& w) const
    { return Weighted_point_2(p,w); }

    Weighted_point_2
    operator()(const FT& x, const FT& y) const
    { return Weighted_point_2(x, y); }

    Weighted_point_2
    operator()(const Point_2& p) const
    { return Weighted_point_2(p,0); }

    const Weighted_point_2&
    operator()(const Weighted_point_2 & wp) const
    { return wp; }
  };

  template <typename K>
  class Construct_weighted_point_3
  {
    typedef typename K::FT                  FT;
    typedef typename K::Point_3             Point_3;
    typedef typename K::Weighted_point_3    Weighted_point_3;
    typedef typename Weighted_point_3::Rep  Rep;
  public:
    typedef Weighted_point_3                result_type;

    Rep
    operator()(Return_base_tag, Origin o) const
    { return Rep(o); }

    Rep
    operator()(Return_base_tag, const Point_3& p, const FT& w) const
    { return Rep(p,w); }

    Rep
    operator()(Return_base_tag, const FT& x, const FT& y, const FT& z) const
    { return Rep(x,y,z); }

    Weighted_point_3
    operator()(Origin o) const
    { return Weighted_point_3(o); }

    Weighted_point_3
    operator()(const Point_3& p, const FT& w) const
    { return Rep(p,w); }

    Weighted_point_3
    operator()(const FT& x, const FT& y, const FT& z) const
    { return Weighted_point_3(x,y,z); }

    Weighted_point_3
    operator()(const Point_3& p) const
    { return Weighted_point_3(p,0); }

    const Weighted_point_3&
    operator()(const Weighted_point_3& wp) const
    { return wp; }
  };


  template <typename K>
  class Construct_projected_point_2
  {
    typedef typename K::Point_2    Point_2;
    typedef typename K::Line_2     Line_2;
  public:
    typedef Point_2                result_type;

    Point_2
    operator()( const Line_2& l, const Point_2& p ) const
    {
      typename K::FT x, y;
      typename K::Construct_point_2 construct_point_2;
      line_project_pointC2(l.a(), l.b(), l.c(), p.x(), p.y(), x, y);
      return construct_point_2(x, y);
    }
  };


  template <typename K>
  class Construct_projected_point_3
  {
    typedef typename K::Point_3    Point_3;
    typedef typename K::Plane_3    Plane_3;
    typedef typename K::Line_3     Line_3;
    typedef typename K::Triangle_3 Triangle_3;
    typedef typename K::Segment_3  Segment_3;
    typedef typename K::Ray_3      Ray_3;
    typedef typename K::FT         FT;
  public:
    typedef Point_3                result_type;

    Point_3
    operator()( const Line_3& l, const Point_3& p ) const
    {
      // projects p on the line l
      FT lpx = l.point().x();
      FT lpy = l.point().y();
      FT lpz = l.point().z();
      FT ldx = l.direction().dx();
      FT ldy = l.direction().dy();
      FT ldz = l.direction().dz();
      FT dpx = p.x()-lpx;
      FT dpy = p.y()-lpy;
      FT dpz = p.z()-lpz;
      FT lambda = (ldx*dpx+ldy*dpy+ldz*dpz) / (ldx*ldx+ldy*ldy+ldz*ldz);
      return Point_3(lpx + lambda * ldx,
                     lpy + lambda * ldy,
                     lpz + lambda * ldz);
    }

    Point_3
    operator()( const Plane_3& h, const Point_3& p ) const
    { return h.rep().projection(p); }

    Point_3
    operator()( const Triangle_3& t, const Point_3& p ) const
    { return CommonKernelFunctors::Construct_projected_point_3<K>()(p,t,K()); }

    Point_3
    operator()( const Segment_3& s, const Point_3& p ) const
    { return CommonKernelFunctors::Construct_projected_point_3<K>()(p,s,K()); }

    Point_3
    operator()( const Ray_3& r, const Point_3& p ) const
    { return CommonKernelFunctors::Construct_projected_point_3<K>()(p,r,K()); }
  };

  template <class K>
  class Construct_radical_line_2
  {
    typedef typename K::Line_2            Line_2;
    typedef typename K::Circle_2          Circle_2;
    typedef typename K::FT                 FT;

  public:

    typedef Line_2 result_type;

    result_type
    operator() (const Circle_2 & c1, const Circle_2 & c2) const
          {
      // Concentric Circles don't have radical line
      CGAL_kernel_precondition (c1.center() != c2.center());
      const FT a = 2*(c2.center().x() - c1.center().x());
      const FT b = 2*(c2.center().y() - c1.center().y());
      const FT c = CGAL::square(c1.center().x()) +
        CGAL::square(c1.center().y()) - c1.squared_radius() -
        CGAL::square(c2.center().x()) -
        CGAL::square(c2.center().y()) + c2.squared_radius();
      return Line_2(a, b, c);
    }
  };

  template <class K>
  class Construct_radical_plane_3
  {
    typedef typename K::Plane_3            Plane_3;
    typedef typename K::Sphere_3           Sphere_3;
    typedef typename K::FT                 FT;

  public:

    typedef Plane_3 result_type;

    result_type
    operator() (const Sphere_3 & s1, const Sphere_3 & s2) const
    {
      // Concentric Spheres don't have radical plane
      CGAL_kernel_precondition (s1.center() != s2.center());
      const FT a = 2*(s2.center().x() - s1.center().x());
      const FT b = 2*(s2.center().y() - s1.center().y());
      const FT c = 2*(s2.center().z() - s1.center().z());
      const FT d = CGAL::square(s1.center().x()) +
        CGAL::square(s1.center().y()) +
        CGAL::square(s1.center().z()) - s1.squared_radius() -
        CGAL::square(s2.center().x()) -
        CGAL::square(s2.center().y()) -
        CGAL::square(s2.center().z()) + s2.squared_radius();
      return Plane_3(a, b, c, d);
    }
  };


  template <typename K>
  class Construct_scaled_vector_2
  {
    typedef typename K::FT         FT;
    typedef typename K::Vector_2   Vector_2;
  public:
    typedef Vector_2               result_type;

    Vector_2
    operator()( const Vector_2& v, const FT& c) const
    {
      return Vector_2(c * v.x(), c * v.y());
    }
  };

  template <typename K>
  class Construct_divided_vector_2
  {
    typedef typename K::FT         FT;
    typedef typename K::Vector_2   Vector_2;
  public:
    typedef Vector_2               result_type;

    Vector_2
    operator()( const Vector_2& v, const FT& c) const
    {
      return Vector_2(v.x()/c, v.y()/c);
    }
  };

  template <typename K>
  class Construct_divided_vector_3
  {
    typedef typename K::FT         FT;
    typedef typename K::Vector_3   Vector_3;
  public:
    typedef Vector_3               result_type;

    Vector_3
    operator()( const Vector_3& v, const FT& c) const
    {
      return Vector_3(v.x()/c, v.y()/c, v.z()/c);
    }
  };

  template <typename K>
  class Construct_scaled_vector_3
  {
    typedef typename K::FT         FT;
    typedef typename K::Vector_3   Vector_3;
  public:
    typedef Vector_3               result_type;

    Vector_3
    operator()( const Vector_3& w, const FT& c) const
    {
      return Vector_3(c * w.x(), c * w.y(), c * w.z());
    }
  };

  template <typename K>
  class Construct_translated_point_2
  {
    typedef typename K::Point_2   Point_2;
    typedef typename K::Vector_2  Vector_2;
  public:
    typedef Point_2               result_type;

    Point_2
    operator()( const Point_2& p, const Vector_2& v) const
    {
      typename K::Construct_point_2 construct_point_2;
      return construct_point_2(p.x() + v.x(), p.y() + v.y());
    }

    Point_2
    operator()( const Origin& , const Vector_2& v) const
    {
      typename K::Construct_point_2 construct_point_2;
      return construct_point_2(v.x(), v.y());
    }
  };

  template <typename K>
  class Construct_translated_point_3
  {
    typedef typename K::Point_3   Point_3;
    typedef typename K::Vector_3  Vector_3;
  public:
    typedef Point_3               result_type;

    Point_3
    operator()( const Point_3& p, const Vector_3& v) const
    {
      typename K::Construct_point_3 construct_point_3;
      return construct_point_3(p.x() + v.x(), p.y() + v.y(), p.z() + v.z());
    }

    Point_3
    operator()( const Origin& , const Vector_3& v) const
    {
      typename K::Construct_point_3 construct_point_3;
      return construct_point_3(v.x(), v.y(), v.z());
    }
  };

  template <typename K>
  class Construct_vector_2
  {
    typedef typename K::RT           RT;
    typedef typename K::FT           FT;
    typedef typename K::Segment_2    Segment_2;
    typedef typename K::Ray_2        Ray_2;
    typedef typename K::Line_2       Line_2;
    typedef typename K::Vector_2     Vector_2;
    typedef typename K::Point_2      Point_2;
    typedef typename K::Direction_2  Direction_2;
    typedef typename Vector_2::Rep   Rep;
  public:
    typedef Vector_2                 result_type;

    Rep // Vector_2
    operator()(Return_base_tag, const Point_2& p, const Point_2& q) const
    { return Rep(q.x() - p.x(), q.y() - p.y()); }

    Rep // Vector_2
    operator()(Return_base_tag, const Origin&, const Point_2& q) const
    { return Rep(q.x(), q.y()); }

    Rep // Vector_2
    operator()(Return_base_tag, const Point_2& p, const Origin& ) const
    { return Rep(-p.x(), -p.y()); }

    Rep // Vector_2
    operator()(Return_base_tag, const Direction_2& d ) const
    { return Rep(d.dx(), d.dy()); }

    Vector_2
    operator()(Return_base_tag, const Segment_2& s) const
    { return s.to_vector(); }

    Vector_2
    operator()(Return_base_tag, const Ray_2& r) const
    { return r.to_vector(); }

    Rep // Vector_2
    operator()(Return_base_tag, const Line_2& l) const
    { return Rep(l.b(), -l.a()); }

    Rep // Vector_2
    operator()(Return_base_tag, Null_vector) const
    { return Rep(FT(0), FT(0)); }

    Rep // Vector_2
    operator()(Return_base_tag, const RT& x, const RT& y) const
    { return Rep(x, y); }

    Rep // Vector_2
    operator()(Return_base_tag, const RT& x, const RT& y, const RT& w) const
    { return Rep(x, y, w); }



    Vector_2
    operator()( const Point_2& p, const Point_2& q) const
    { return this->operator()(Return_base_tag(), p, q); }

    Vector_2
    operator()( const Origin& o, const Point_2& q) const
    { return this->operator()(Return_base_tag(), o, q); }

    Vector_2
    operator()( const Point_2& p, const Origin& o) const
    { return this->operator()(Return_base_tag(), p, o); }

    Vector_2
    operator()( const Direction_2& d ) const
    { return this->operator()(Return_base_tag(), d); }

    Vector_2
    operator()( const Segment_2& s) const
    { return this->operator()(Return_base_tag(), s); }

    Vector_2
    operator()( const Ray_2& r) const
    { return this->operator()(Return_base_tag(), r); }

    Vector_2
    operator()( const Line_2& l) const
    { return this->operator()(Return_base_tag(), l); }

    Vector_2
    operator()( Null_vector n) const
    { return this->operator()(Return_base_tag(), n); }

    Vector_2
    operator()( const RT& x, const RT& y) const
    { return this->operator()(Return_base_tag(), x, y); }

    Vector_2
    operator()( const RT& x, const RT& y, const RT& w) const
    { return this->operator()(Return_base_tag(), x, y, w); }
  };

  template <typename K>
  class Construct_vector_3
  {
    typedef typename K::RT           RT;
    typedef typename K::FT           FT;
    typedef typename K::Segment_3    Segment_3;
    typedef typename K::Direction_3  Direction_3;
    typedef typename K::Ray_3        Ray_3;
    typedef typename K::Line_3       Line_3;
    typedef typename K::Vector_3     Vector_3;
    typedef typename K::Point_3      Point_3;
    typedef typename Vector_3::Rep   Rep;
  public:
    typedef Vector_3                 result_type;

    Rep // Vector_3
    operator()(Return_base_tag, const Point_3& p, const Point_3& q) const
    {
      return Rep(q.x() - p.x(), q.y() - p.y(), q.z() - p.z());
    }

    Rep // Vector_3
    operator()(Return_base_tag, const Origin&, const Point_3& q) const
    {
      return Rep(q.x(), q.y(), q.z());
    }

    Rep // Vector_3
    operator()(Return_base_tag, const Point_3& p, const Origin&) const
    {
      return Rep(- p.x(), - p.y(), - p.z());
    }

    Rep // Vector_3
    operator()(Return_base_tag, const Direction_3& d) const
    { return d.rep().to_vector(); }

    Rep // Vector_3
    operator()(Return_base_tag, const Segment_3& s) const
    { return s.rep().to_vector(); }

    Rep // Vector_3
    operator()(Return_base_tag, const Ray_3& r) const
    { return r.rep().to_vector(); }

    Rep // Vector_3
    operator()(Return_base_tag, const Line_3& l) const
    { return l.rep().to_vector(); }

    Rep // Vector_3
    operator()(Return_base_tag, const Null_vector&) const
    { return Rep(FT(0), FT(0), FT(0)); }

    Rep // Vector_3
    operator()(Return_base_tag, const RT& x, const RT& y, const RT& z) const
    { return Rep(x, y, z); }

    Rep // Vector_3
    operator()(Return_base_tag, const RT& x, const RT& y, const RT& z, const RT& w) const
    { return Rep(x, y, z, w); }


    Vector_3
    operator()( const Point_3& p, const Point_3& q) const
    { return this->operator()(Return_base_tag(), p, q); }

    Vector_3
    operator()( const Origin& o, const Point_3& q) const
    { return this->operator()(Return_base_tag(), o, q); }

    Vector_3
    operator()( const Point_3& p, const Origin& q) const
    { return this->operator()(Return_base_tag(), p, q); }

    Vector_3
    operator()( const Direction_3& d) const
    { return this->operator()(Return_base_tag(), d); }

    Vector_3
    operator()( const Segment_3& s) const
    { return this->operator()(Return_base_tag(), s); }

    Vector_3
    operator()( const Ray_3& r) const
    { return this->operator()(Return_base_tag(), r); }

    Vector_3
    operator()( const Line_3& l) const
    { return this->operator()(Return_base_tag(), l); }

    Vector_3
    operator()( const Null_vector& n) const
    { return this->operator()(Return_base_tag(), n); }

    Vector_3
    operator()( int x, int y, int z) const
    { return this->operator()(Return_base_tag(), x, y, z); }

    Vector_3
    operator()( const RT& x, const RT& y, const RT& z) const
    { return this->operator()(Return_base_tag(), x, y, z); }

    Vector_3
    operator()( const RT& x, const RT& y, const RT& z, const RT& w) const
    { return this->operator()(Return_base_tag(), x, y, z, w); }
  };

  template <typename K>
  class Construct_vertex_2
  {
    typedef typename K::Point_2          Point_2;
    typedef typename K::Segment_2        Segment_2;
    typedef typename K::Iso_rectangle_2  Iso_rectangle_2;
    typedef typename K::Triangle_2       Triangle_2;
  public:
    template<class>
    struct result {
      typedef const Point_2& type;
    };

    template<typename F>
    struct result<F(Iso_rectangle_2, int)> {
      typedef Point_2 type;
    };

    const Point_2 &
    operator()( const Segment_2& s, int i) const
    { return s.vertex(i); }

    const Point_2 &
    operator()( const Triangle_2& t, int i) const
    { return t.rep().vertex(i); }

    Point_2
    operator()( const Iso_rectangle_2& r, int i) const
    {
      switch (i%4) {
      case 0: return (r.min)();
      case 1: return Point_2(r.xmax(), r.ymin());
      case 2: return (r.max)();
      default: return Point_2(r.xmin(), r.ymax());
      }
    }
  };

} //namespace CartesianKernelFunctors

namespace CartesianKernelFunctors {

  template <typename K>
  class Coplanar_orientation_3
  {
    typedef typename K::Point_3      Point_3;
#ifdef CGAL_kernel_exactness_preconditions
    typedef typename K::Coplanar_3   Coplanar_3;
    typedef typename K::Collinear_3  Collinear_3;
    Coplanar_3 cp;
    Collinear_3 cl;
#endif // CGAL_kernel_exactness_preconditions
  public:
    typedef typename K::Orientation  result_type;

#ifdef CGAL_kernel_exactness_preconditions
    Coplanar_orientation_3() {}
    Coplanar_orientation_3(const Coplanar_3& cp_, const Collinear_3& cl_)
      : cp(cp_), cl(cl_)
    {}
#endif // CGAL_kernel_exactness_preconditions

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      return coplanar_orientationC3(p.x(), p.y(), p.z(),
                                    q.x(), q.y(), q.z(),
                                    r.x(), r.y(), r.z());
    }

    result_type
    operator()( const Point_3& p, const Point_3& q,
                const Point_3& r, const Point_3& s) const
    {
      // p,q,r,s supposed to be coplanar
      // p,q,r supposed to be non collinear
      // tests whether s is on the same side of p,q as r
      // returns :
      // COLLINEAR if pqr collinear
      // POSITIVE if qrp and qrs have the same orientation
      // NEGATIVE if qrp and qrs have opposite orientations
      CGAL_kernel_exactness_precondition( ! cl(p, q, r) );
      CGAL_kernel_exactness_precondition( cp(p, q, r, s) );
      return coplanar_orientationC3(p.x(), p.y(), p.z(),
                                    q.x(), q.y(), q.z(),
                                    r.x(), r.y(), r.z(),
                                    s.x(), s.y(), s.z());
    }
  };

  template <typename K>
  class Coplanar_side_of_bounded_circle_3
  {
    typedef typename K::Point_3   Point_3;
#ifdef CGAL_kernel_exactness_preconditions
    typedef typename K::Coplanar_3   Coplanar_3;
    typedef typename K::Collinear_3  Collinear_3;
    Coplanar_3 cp;
    Collinear_3 cl;
#endif // CGAL_kernel_exactness_preconditions
  public:
    typedef typename K::Bounded_side     result_type;

#ifdef CGAL_kernel_exactness_preconditions
    Coplanar_side_of_bounded_circle_3() {}
    Coplanar_side_of_bounded_circle_3(const Coplanar_3& cp_,
                                      const Collinear_3& cl_)
      : cp(cp_), cl(cl_)
    {}
#endif // CGAL_kernel_exactness_preconditions

    result_type
    operator()( const Point_3& p, const Point_3& q,
                const Point_3& r, const Point_3& t) const
    {
      // p,q,r,t are supposed to be coplanar.
      // p,q,r determine an orientation of this plane (not collinear).
      // returns the equivalent of side_of_bounded_circle(p,q,r,t)
      // in this plane
      CGAL_kernel_exactness_precondition( cp(p,q,r,t) );
      CGAL_kernel_exactness_precondition( !cl(p,q,r) );
      return coplanar_side_of_bounded_circleC3(p.x(), p.y(), p.z(),
                                               q.x(), q.y(), q.z(),
                                               r.x(), r.y(), r.z(),
                                               t.x(), t.y(), t.z());
    }
  };

  template <typename K>
  class Equal_xy_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef typename K::Boolean    result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    {
      return CGAL_AND( p.x() == q.x() , p.y() == q.y() );
    }
  };

  template <typename K>
  class Equal_x_2
  {
    typedef typename K::Point_2    Point_2;
  public:
    typedef typename K::Boolean    result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    { return p.x() == q.x(); }
  };

  template <typename K>
  class Equal_x_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef typename K::Boolean    result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return p.x() == q.x(); }
  };

  template <typename K>
  class Equal_y_2
  {
    typedef typename K::Point_2    Point_2;
  public:
    typedef typename K::Boolean    result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    { return p.y() == q.y(); }
  };

  template <typename K>
  class Equal_y_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef typename K::Boolean    result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return p.y() == q.y(); }
  };

  template <typename K>
  class Equal_z_3
  {
    typedef typename K::Point_3    Point_3;
  public:
    typedef typename K::Boolean    result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return p.z() == q.z(); }
  };

  template <typename K>
  class Has_on_3
  {
    typedef typename K::FT               FT;
    typedef typename K::Point_3          Point_3;
    typedef typename K::Vector_3         Vector_3;
    typedef typename K::Line_3           Line_3;
    typedef typename K::Ray_3            Ray_3;
    typedef typename K::Segment_3        Segment_3;
    typedef typename K::Plane_3          Plane_3;
    typedef typename K::Triangle_3       Triangle_3;
    typedef typename K::Circle_3         Circle_3;
    typedef typename K::Sphere_3         Sphere_3;
  public:
    typedef typename K::Boolean          result_type;

    result_type
    operator()( const Line_3& l, const Point_3& p) const
    { return l.rep().has_on(p); }

    result_type
    operator()( const Ray_3& r, const Point_3& p) const
    { return r.rep().has_on(p); }

    result_type
    operator()( const Segment_3& s, const Point_3& p) const
    { return s.has_on(p); }

    result_type
    operator()( const Plane_3& pl, const Point_3& p) const
    { return pl.rep().has_on(p); }

    result_type
    operator()( const Plane_3& pl, const Line_3& l) const
    { return pl.rep().has_on(l); }

    result_type
    operator()( const Triangle_3& t, const Point_3& p) const
    {
      Point_3  o  = t.vertex(0) + t.supporting_plane().orthogonal_vector();
      Vector_3 v0 = t.vertex(0)-o,
               v1 = t.vertex(1)-o,
               v2 = t.vertex(2)-o;

      FT alpha, beta, gamma;
      Cartesian_internal::solve(v0, v1, v2, p-o, alpha, beta, gamma);
      return (alpha >= FT(0)) && (beta >= FT(0)) && (gamma >= FT(0))
          && ((alpha+beta+gamma == FT(1)));
    }

    result_type
    operator()(const Circle_3 &a, const Point_3 &p) const
    { return a.rep().has_on(p); }

    result_type
    operator()(const Sphere_3 &a, const Circle_3 &p) const
    { return a.rep().has_on(p); }

    result_type
    operator()(const Sphere_3 &a, const Point_3 &p) const
    { return a.rep().has_on(p); }

    result_type
    operator()(const Plane_3 &a, const Circle_3 &p) const
    { return a.rep().has_on(p); }


  };

  template <typename K>
  class Less_distance_to_point_2
  {
    typedef typename K::Point_2   Point_2;
  public:
    typedef typename K::Boolean   result_type;

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      return has_smaller_dist_to_pointC2(p.x(), p.y(),
                                         q.x(), q.y(),
                                         r.x(), r.y());
    }
  };

  template <typename K>
  class Less_distance_to_point_3
  {
    typedef typename K::Point_3   Point_3;
  public:
    typedef typename K::Boolean   result_type;

    result_type
    operator()(const Point_3& p, const Point_3& q, const Point_3& r) const
    {
      return has_smaller_dist_to_pointC3(p.x(), p.y(), p.z(),
                                         q.x(), q.y(), q.z(),
                                         r.x(), r.y(), r.z());
    }
  };


  template <typename K>
  class Less_signed_distance_to_line_2
  {
    typedef typename K::Point_2   Point_2;
    typedef typename K::Line_2    Line_2;
    typedef typename K::Equal_2   Equal_2;
  public:
    typedef typename K::Boolean   result_type;

    result_type
    operator()(const Point_2& a, const Point_2& b,
               const Point_2& c, const Point_2& d) const
    {
      CGAL_kernel_precondition_code(Equal_2 equal;)
      CGAL_kernel_precondition(! equal(a,b));
      return cmp_signed_dist_to_lineC2( a.x(), a.y(),
                                        b.x(), b.y(),
                                        c.x(), c.y(),
                                        d.x(), d.y()) == SMALLER;
    }

    result_type
    operator()(const Line_2& l, const Point_2& p, const Point_2& q) const
    {
      return has_smaller_signed_dist_to_directionC2(l.a(), l.b(),
                                                    p.x(), p.y(),
                                                    q.x(), q.y());
    }
  };

  template <typename K>
  class Less_signed_distance_to_plane_3
  {
    typedef typename K::Point_3       Point_3;
    typedef typename K::Plane_3       Plane_3;
    typedef typename K::Collinear_3   Collinear_3;
  public:
    typedef typename K::Boolean       result_type;

    result_type
    operator()( const Plane_3& h, const Point_3& p, const Point_3& q) const
    {
      return has_smaller_signed_dist_to_directionC3(h.a(), h.b(), h.c(),
                                                    p.x(), p.y(), p.z(),
                                                    q.x(), q.y(), q.z());
    }

    result_type
    operator()( const Point_3& hp, const Point_3& hq,  const Point_3& hr,
                const Point_3& p, const Point_3& q) const
    {
      CGAL_kernel_precondition_code(Collinear_3 collinear_3;)
      CGAL_kernel_precondition(! collinear_3(hp, hq, hr));
      return has_smaller_signed_dist_to_planeC3(hp.x(), hp.y(), hp.z(),
                                                hq.x(), hq.y(), hq.z(),
                                                hr.x(), hr.y(), hr.z(),
                                                p.x(),  p.y(),  p.z(),
                                                q.x(),  q.y(),  q.z());;
    }
  };

  template <typename K>
  class Less_xyz_3
  {
    typedef typename K::Point_3         Point_3;
    typedef typename K::Compare_xyz_3   Compare_xyz_3;
    Compare_xyz_3 c;
  public:
    typedef typename K::Boolean         result_type;

    Less_xyz_3() {}
    Less_xyz_3(const Compare_xyz_3& c_) : c(c_) {}

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return c(p, q) == SMALLER; }
  };

  template <typename K>
  class Less_xy_2
  {
    typedef typename K::Point_2        Point_2;
    typedef typename K::Compare_xy_2   Compare_xy_2;
    Compare_xy_2 c;
  public:
    typedef typename K::Boolean        result_type;

    Less_xy_2() {}
    Less_xy_2(const Compare_xy_2& c_) : c(c_) {}

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    { return c(p, q) == SMALLER; }
  };

  template <typename K>
  class Less_xy_3
  {
    typedef typename K::Point_3        Point_3;
    typedef typename K::Compare_xy_3   Compare_xy_3;
    Compare_xy_3 c;
  public:
    typedef typename K::Boolean        result_type;

    Less_xy_3() {}
    Less_xy_3(const Compare_xy_3& c_) : c(c_) {}

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return c(p, q) == SMALLER; }
  };

  template <typename K>
  class Less_x_2
  {
    typedef typename K::Point_2        Point_2;
  public:
    typedef typename K::Boolean        result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    { return p.x() < q.x(); }
  };

  template <typename K>
  class Less_x_3
  {
    typedef typename K::Point_3        Point_3;
  public:
    typedef typename K::Boolean        result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return p.x() < q.x(); }
  };

  template <typename K>
  class Less_yx_2
  {
    typedef typename K::Point_2        Point_2;
  public:
    typedef typename K::Boolean        result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    {
      return compare_lexicographically_xyC2(p.y(), p.x(),
                                            q.y(), q.x()) == SMALLER;
    }
  };

  template <typename K>
  class Less_y_2
  {
    typedef typename K::Point_2        Point_2;
  public:
    typedef typename K::Boolean        result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q) const
    { return p.y() < q.y(); }
  };

  template <typename K>
  class Less_y_3
  {
    typedef typename K::Point_3        Point_3;
  public:
    typedef typename K::Boolean        result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return p.y() < q.y(); }
  };

  template <typename K>
  class Less_z_3
  {
    typedef typename K::Point_3        Point_3;
  public:
    typedef typename K::Boolean        result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q) const
    { return p.z() < q.z(); }
  };

  template <typename K>
  class Orientation_2
  {
    typedef typename K::Point_2       Point_2;
    typedef typename K::Vector_2      Vector_2;
    typedef typename K::Circle_2      Circle_2;
  public:
    typedef typename K::Orientation   result_type;

    result_type
    operator()(const Point_2& p, const Point_2& q, const Point_2& r) const
    {
      return orientationC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
    }

    result_type
    operator()(const Vector_2& u, const Vector_2& v) const
    {
      return orientationC2(u.x(), u.y(), v.x(), v.y());
    }

    result_type
    operator()(const Circle_2& c) const
    {
      return c.rep().orientation();
    }
  };

  template <typename K>
  class Orientation_3
  {
    typedef typename K::Point_3        Point_3;
    typedef typename K::Vector_3       Vector_3;
    typedef typename K::Tetrahedron_3  Tetrahedron_3;
    typedef typename K::Sphere_3       Sphere_3;
  public:
    typedef typename K::Orientation    result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q,
                const Point_3& r, const Point_3& s) const
    {
      return orientationC3(p.x(), p.y(), p.z(),
                           q.x(), q.y(), q.z(),
                           r.x(), r.y(), r.z(),
                           s.x(), s.y(), s.z());
    }

    result_type
    operator()( const Vector_3& u, const Vector_3& v, const Vector_3& w) const
    {
      return orientationC3(u.x(), u.y(), u.z(),
                           v.x(), v.y(), v.z(),
                           w.x(), w.y(), w.z());
    }

    result_type
    operator()( const Tetrahedron_3& t) const
    {
      return t.rep().orientation();
    }

    result_type
    operator()(const Sphere_3& s) const
    {
      return s.rep().orientation();
    }
  };

  template < typename K >
  class Power_side_of_oriented_power_circle_2
  {
  public:
    typedef typename K::Weighted_point_2         Weighted_point_2;
    typedef typename K::Oriented_side            Oriented_side;

    typedef Oriented_side                        result_type;

    Oriented_side operator()(const Weighted_point_2& p,
                             const Weighted_point_2& q,
                             const Weighted_point_2& r,
                             const Weighted_point_2& t) const
    {
      //CGAL_kernel_precondition( ! collinear(p, q, r) );
      return power_side_of_oriented_power_circleC2(p.x(), p.y(), p.weight(),
                                                   q.x(), q.y(), q.weight(),
                                                   r.x(), r.y(), r.weight(),
                                                   t.x(), t.y(), t.weight());
    }

    // The methods below are currently undocumented because the definition of
    // orientation is unclear for 2 and 1 point configurations in a 2D space.

    // One should be (very) careful with the order of vertices when using them,
    // as swapping points will change the result and one must therefore have a
    // precise idea of what is the positive orientation in the full space.
    // For example, these functions are (currently) used safely in the regular
    // triangulations classes because we always call them on vertices of
    // triangulation cells, which are always positively oriented.

    Oriented_side operator()(const Weighted_point_2& p,
                             const Weighted_point_2& q,
                             const Weighted_point_2& t) const
    {
      //CGAL_kernel_precondition( collinear(p, q, r) );
      //CGAL_kernel_precondition( p.point() != q.point() );
      return power_side_of_oriented_power_circleC2(p.point().x(), p.y(), p.weight(),
                                                   q.x(), q.y(), q.weight(),
                                                   t.x(), t.y(), t.weight());
    }

    Oriented_side operator()(const Weighted_point_2& p,
                             const Weighted_point_2& t) const
    {
      //CGAL_kernel_precondition( p.point() == r.point() );
      Comparison_result r = CGAL::compare(p.weight(), t.weight());
      if(r == LARGER)    return ON_NEGATIVE_SIDE;
      else if (r == SMALLER) return ON_POSITIVE_SIDE;
      return ON_ORIENTED_BOUNDARY;
    }
  };

  template <typename K>
  class Oriented_side_2
  {
    typedef typename K::Point_2        Point_2;
    typedef typename K::Circle_2       Circle_2;
    typedef typename K::Line_2         Line_2;
    typedef typename K::Triangle_2     Triangle_2;
  public:
    typedef typename K::Oriented_side  result_type;

    result_type
    operator()( const Circle_2& c, const Point_2& p) const
    { return enum_cast<Oriented_side>(c.bounded_side(p)) * c.orientation(); }

    result_type
    operator()( const Line_2& l, const Point_2& p) const
    { return side_of_oriented_lineC2(l.a(), l.b(), l.c(), p.x(), p.y()); }

    result_type
    operator()( const Triangle_2& t, const Point_2& p) const
    {
      typename K::Collinear_are_ordered_along_line_2
        collinear_are_ordered_along_line;
      typename K::Orientation_2 orientation;
      // depends on the orientation of the vertices
      typename K::Orientation
                  o1 = orientation(t.vertex(0), t.vertex(1), p),
                  o2 = orientation(t.vertex(1), t.vertex(2), p),
                  o3 = orientation(t.vertex(2), t.vertex(3), p),
                  ot = orientation(t.vertex(0), t.vertex(1), t.vertex(2));

      if (o1 == ot && o2 == ot && o3 == ot) // ot cannot be COLLINEAR
        return ot;
      return
        (o1 == COLLINEAR
         && collinear_are_ordered_along_line(t.vertex(0), p, t.vertex(1))) ||
        (o2 == COLLINEAR
         && collinear_are_ordered_along_line(t.vertex(1), p, t.vertex(2))) ||
        (o3 == COLLINEAR
         && collinear_are_ordered_along_line(t.vertex(2), p, t.vertex(3)))
        ? result_type(ON_ORIENTED_BOUNDARY)
        : opposite(ot);
    }
  };

  template <typename K>
  class Side_of_bounded_circle_2
  {
    typedef typename K::Point_2        Point_2;
  public:
    typedef typename K::Bounded_side   result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q, const Point_2& t) const
    {
      return side_of_bounded_circleC2(p.x(), p.y(),
                                      q.x(), q.y(),
                                      t.x(), t.y());
    }

    result_type
    operator()( const Point_2& p, const Point_2& q,
                const Point_2& r, const Point_2& t) const
    {
      return side_of_bounded_circleC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y(),
                                      t.x(), t.y());
    }
  };

  template <typename K>
  class Side_of_bounded_sphere_3
  {
    typedef typename K::Point_3        Point_3;
  public:
    typedef typename K::Bounded_side   result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q, const Point_3& test) const
    {
      return side_of_bounded_sphereC3(p.x(), p.y(), p.z(),
                                      q.x(), q.y(), q.z(),
                                      test.x(), test.y(), test.z());
    }

    result_type
    operator()( const Point_3& p, const Point_3& q,
                const Point_3& r, const Point_3& test) const
    {
      return side_of_bounded_sphereC3(p.x(), p.y(), p.z(),
                                      q.x(), q.y(), q.z(),
                                      r.x(), r.y(), r.z(),
                                      test.x(), test.y(), test.z());
    }

    result_type
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
                const Point_3& s, const Point_3& test) const
    {
      return side_of_bounded_sphereC3(p.x(), p.y(), p.z(),
                                      q.x(), q.y(), q.z(),
                                      r.x(), r.y(), r.z(),
                                      s.x(), s.y(), s.z(),
                                      test.x(), test.y(), test.z());
    }
  };

  template <typename K>
  class Side_of_oriented_circle_2
  {
    typedef typename K::Point_2        Point_2;
  public:
    typedef typename K::Oriented_side  result_type;

    result_type
    operator()( const Point_2& p, const Point_2& q,
                const Point_2& r, const Point_2& t) const
    {
      return side_of_oriented_circleC2(p.x(), p.y(),
                                       q.x(), q.y(),
                                       r.x(), r.y(),
                                       t.x(), t.y());
    }
  };

  template <typename K>
  class Side_of_oriented_sphere_3
  {
    typedef typename K::Point_3        Point_3;
  public:
    typedef typename K::Oriented_side  result_type;

    result_type
    operator()( const Point_3& p, const Point_3& q, const Point_3& r,
                const Point_3& s, const Point_3& test) const
    {
      return side_of_oriented_sphereC3(p.x(), p.y(), p.z(),
                                       q.x(), q.y(), q.z(),
                                       r.x(), r.y(), r.z(),
                                       s.x(), s.y(), s.z(),
                                       test.x(), test.y(), test.z());
    }
  };


} // namespace CartesianKernelFunctors

} //namespace CGAL

#endif // CGAL_CARTESIAN_FUNCTION_OBJECTS_H
