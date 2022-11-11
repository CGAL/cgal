// Copyright (c) 2016-2018  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     :   Iordan Iordanov
//                   Monique Teillaud
//

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H

#include <CGAL/license/Periodic_4_hyperbolic_triangulation_2.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/Cartesian.h>

#include <boost/tuple/tuple.hpp>
#include <boost/variant.hpp>

#include <utility>

namespace CGAL {

namespace internal {


template <typename K, typename Predicate_>
class Hyperbolic_traits_with_translations_2_adaptor
  : Predicate_
{
  typedef K                                               Kernel;
  typedef Predicate_                                      Predicate;

  typedef typename Kernel::FT                             FT;
  typedef typename Kernel::Point_2                        Point;
  typedef typename Kernel::Hyperbolic_translation         Hyperbolic_translation;

  // Use the construct_point_2 predicate from the kernel to convert the periodic points to Euclidean points
  typedef typename Kernel::Construct_hyperbolic_point_2   Construct_hyperbolic_point_2;

public:
  typedef typename Predicate::result_type                 result_type;

  using Predicate::operator();

  Hyperbolic_traits_with_translations_2_adaptor(const Predicate_ pred = Predicate_()) : Predicate_(pred) {}

  result_type operator()(const Point& p0, const Point& p1,
                         const Hyperbolic_translation& o0, const Hyperbolic_translation& o1) const
  {
    return operator()(pp(p0, o0), pp(p1, o1));
  }

  result_type operator()(const Point& p0, const Point& p1, const Point& p2,
                         const Hyperbolic_translation& o0, const Hyperbolic_translation& o1, const Hyperbolic_translation& o2) const
  {
    return operator()(pp(p0, o0), pp(p1, o1), pp(p2, o2));
  }

  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3,
                         const Hyperbolic_translation& o0, const Hyperbolic_translation& o1,
                         const Hyperbolic_translation& o2, const Hyperbolic_translation& o3) const
  {
    return operator()(pp(p0, o0), pp(p1, o1), pp(p2, o2), pp(p3, o3));
  }

  // Added for Side_of_hyperbolic_triangle_2
  template <typename T>
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3, T& t) const
  {
    return operator()(p0, p1, p2, p3, t);
  }

  template <typename T>
  result_type operator()(const Point& p0, const Point& p1,
                         const Point& p2, const Point& p3,
                         const Hyperbolic_translation& o0, const Hyperbolic_translation& o1,
                         const Hyperbolic_translation& o2, const Hyperbolic_translation& o3, T& t) const
  {
    return operator()(pp(p0, o0), pp(p1, o1), pp(p2, o2), pp(p3, o3), t);
  }
private:
  Point pp(const Point &p, const Hyperbolic_translation &o) const
  {
    return Construct_hyperbolic_point_2()(p, o);
  }
};

template < typename K, typename Construct_point_base>
class Periodic_4_construct_hyperbolic_point_2
  : public Construct_point_base
{

private:
  typedef K                                               Kernel;
  typedef typename Kernel::FT                             NT;
  typedef typename Kernel::Point_2                        Point;
  typedef typename Kernel::Hyperbolic_translation         Hyperbolic_translation;

public:
  typedef Point                                           result_type;

  using Construct_point_base::operator();

  Periodic_4_construct_hyperbolic_point_2() { }

  Point operator()(const Point& pt, const Hyperbolic_translation& tr) const
  {
    if(tr.is_identity())
      return operator()(pt);

    Point p = operator()(pt);
    std::pair<NT, NT> a, b;
    a = tr.alpha();
    b = tr.beta();

    // Prepare variables
    NT ax(a.first);
    NT bx(b.first);
    NT zx(p.x());
    NT ay(a.second);
    NT by(b.second);
    NT zy(p.y());

    // Compute parts of fraction
    NT rn(ax*zx - ay*zy + bx); // real part of numerator
    NT in(ay*zx + ax*zy + by); // imaginary part of numerator
    NT rd(bx*zx + by*zy + ax); // real part of denominator
    NT id(bx*zy - by*zx - ay); // imaginary part of denominator

    // The denominator cannot be zero
    NT den(rd*rd + id*id);
    CGAL_assertion(den != NT(0));

    // Compute real and imaginary part of result
    NT rx((rn*rd + in*id));
    NT ry((in*rd - rn*id));

    return operator()(rx/den, ry/den);
  }

  Point operator() (const Point& p) const
  {
    return p;
  }

};


template <typename Traits>
class Compute_approximate_hyperbolic_diameter
  {

    typedef typename Traits:: Euclidean_line_2     Euclidean_line_2;
    typedef typename Traits::Circle_2              Circle_2;
    typedef typename Traits::Hyperbolic_point_2    Hyperbolic_point_2;
  public:
    typedef double result_type;

    Compute_approximate_hyperbolic_diameter(const Traits& gt = Traits()) : _gt(gt) {}

    result_type operator()(const Hyperbolic_point_2& p1,
                           const Hyperbolic_point_2& p2,
                           const Hyperbolic_point_2& p3)
    {

      Circle_2 c = _gt.construct_circle_2_object()(p1, p2, p3);

      Hyperbolic_point_2 p0 = _gt.construct_point_2_object()(0, 0);
      Circle_2           c0 = _gt.construct_circle_2_object()(p0, 1);
      Euclidean_line_2  ell = _gt.construct_line_2_object()(p0, c.center());

      if(ell.is_degenerate())
        return 5.; // for the Bolza surface, this just needs to be larger than half the systole (which is ~1.5)

      std::pair<Hyperbolic_point_2, Hyperbolic_point_2> res1 = _gt.construct_inexact_intersection_2_object()(c0, ell);
      std::pair<Hyperbolic_point_2, Hyperbolic_point_2> res2 = _gt.construct_inexact_intersection_2_object()(c , ell);

      Hyperbolic_point_2 a = res1.first;
      Hyperbolic_point_2 b = res1.second;

      Hyperbolic_point_2 p = res2.first;
      Hyperbolic_point_2 q = res2.second;

      double aq = CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(a, q)));
      double pb = CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(p, b)));
      double ap = CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(a, p)));
      double qb = CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(q, b)));

      double hyperdist = std::fabs(std::log(CGAL::to_double((aq*pb)/(ap*qb))));

      return hyperdist;
    }

  private:
    const Traits& _gt;
  };


  template <typename Traits>
  class Construct_hyperbolic_line_2
  {
    typedef typename Traits::Weighted_point_2                    Weighted_point_2;
    typedef typename Traits::Hyperbolic_point_2                  Hyperbolic_point_2;
    typedef typename Traits::Hyperbolic_segment_2                Hyperbolic_segment_2;
    typedef typename Traits::Point_2                             Bare_point;
    typedef typename Traits::Circle_2                            Circle_2;
    typedef typename Traits::Circular_arc_2                      Circular_arc_2;
    typedef typename Traits::FT                                  FT;

  public:
    Construct_hyperbolic_line_2(const Traits gt = Traits()) : _gt(gt) { }

    Hyperbolic_segment_2 operator()(const Hyperbolic_point_2& p, const Hyperbolic_point_2& q)
    {
      Origin o;
      if(_gt.collinear_2_object()(p, q, _gt.construct_hyperbolic_point_2_object(o)))
      {
        std::pair<Hyperbolic_point_2,Hyperbolic_point_2> inters = _gt.construct_intersection_2_object()(_gt.construct_line_2_object()(p,q), _gt.construct_circle_2_object()(_gt.construct_hyperbolic_point_2_object()(0,0),1));
        return Euclidean_segment_2(inters.first, inters.second);
      }

      Weighted_point_2 wp(p);
      Weighted_point_2 wq(q);
      Weighted_point_2 wo(_gt.construct_hyperbolic_point_2_object()(o), FT(1)); // Poincaré circle

      Bare_point center = _gt.construct_weighted_circumcenter_2_object()(wp, wo, wq);
      FT sq_radius = _gt.compute_squared_Euclidean_distance_2_object()(p, center);

      Circle_2 circle = _gt.construct_circle_2_object()(center, sq_radius);
      std::pair<Hyperbolic_point_2,Hyperbolic_point_2> inters = _gt.construct_intersection_2_object()(circle, _gt.construct_circle_2_object()(_gt.construct_hyperbolic_point_2_object()(0,0),1));

      if(_gt.orientation_2_object()(circle.center(), inters.first, inters.second) == POSITIVE) {
        return Circular_arc_2(circle, inters.first, inters.second);
      } else {
        return Circular_arc_2(circle, inters.second, inters.first);
      }
    }

  private:
    const Traits& _gt;
  };



  template <typename Traits>
  class Construct_inexact_intersection_2
  {
    typedef typename Traits::FT                                 FT;
    typedef typename Traits::Circle_2                           Circle_2;
    typedef typename Traits::Euclidean_line_2                   Euclidean_line_2;
    typedef typename Traits::Hyperbolic_point_2                 Hyperbolic_point_2;
    typedef typename Traits::Hyperbolic_segment_2               Hyperbolic_segment_2;
    typedef typename Traits::Euclidean_segment_2                Euclidean_segment_2;
    typedef typename Traits::Circular_arc_2                     Circular_arc_2;

  public:
    Construct_inexact_intersection_2(const Traits& gt = Traits()) : _gt(gt) {}

    Hyperbolic_point_2 operator()(const Euclidean_line_2& ell1, const Euclidean_line_2& ell2)
    {
      if(std::fabs(CGAL::to_double(ell1.b())) < 1e-16)
        std::swap(ell1, ell2);

      double a1 = CGAL::to_double(ell1.a()), b1 = CGAL::to_double(ell1.b()), c1 = CGAL::to_double(ell1.c());
      double a2 = CGAL::to_double(ell2.a()), b2 = CGAL::to_double(ell2.b()), c2 = CGAL::to_double(ell2.c());

      CGAL_assertion(std::fabs(b1) > 1e-16);
      if(std::fabs(b2) > 1e-16)
        CGAL_assertion(std::fabs(a1/b1 - a2/b2) > 1e-16);

      double lambda1 = -a1/b1;
      double mu1 = -c1/b1;
      double x = (-c2 - mu1*b2)/(a2 + lambda1*b2);
      double y = lambda1*x + mu1;
      return Hyperbolic_point_2(x, y);
    }

    std::pair<Hyperbolic_point_2, Hyperbolic_point_2> operator()(const Euclidean_line_2& ell,
                                                                 const Circle_2& cc)
    {
      double a = CGAL::to_double(ell.a()),
             b = CGAL::to_double(ell.b()),
             c = CGAL::to_double(ell.c());
      double p = CGAL::to_double(cc.center().x()),
             q = CGAL::to_double(cc.center().y()),
             r2 = CGAL::to_double(cc.squared_radius());

      double A, B, C, D;
      double x1, y1, x2, y2;
      if(std::fabs(a) < 1e-16)
      {
        y1 = -c/b;  y2 = -c/b;
        A = b*p;
        D = -b*b*q*q + b*b*r2 - 2.*b*c*q - c*c;
        x1 = (A + CGAL::sqrt(D))/b;
        x2 = (A - CGAL::sqrt(D))/b;
      }
      else if(std::fabs(b) < 1e-16)
      {
        x1 = -c/a;  x2 = -c/a;
        A = q*a;
        D = -a*a*p*p + r2*a*a - 2.*a*c*p - c*c;
        y1 = (A + CGAL::sqrt(D))/a;
        y2 = (A - CGAL::sqrt(D))/a;
      }
      else
      {
        A = a*a*q - a*b*p-b*c;
        C = (-b*q - c)*a*a + b*b*p*a;
        D = -a*a*(b*b*q*q + 2.*q*(p*a + c)*b - b*b*r2 + (p*p - r2)*a*a + 2.*a*c*p + c*c);
        B = a*a + b*b;

        y1 = (A + CGAL::sqrt(D))/B;
        y2 = (A - CGAL::sqrt(D))/B;
        x1 = (C - b*CGAL::sqrt(D))/(a*(a*a + b*b));
        x2 = (C + b*CGAL::sqrt(D))/(a*(a*a + b*b));
      }

      Hyperbolic_point_2 p1(x1, y1);
      Hyperbolic_point_2 p2(x2, y2);

      return std::make_pair(p1, p2);
    }

    std::pair<Hyperbolic_point_2, Hyperbolic_point_2> operator()(const Circle_2& c,
                                                                 const Euclidean_line_2& ell)
    {
      return operator()(ell, c);
    }

    std::pair<Hyperbolic_point_2, Hyperbolic_point_2> operator()(const Circle_2& c1,
                                                                 const Circle_2& c2)
    {
      double xa = CGAL::to_double(c1.center().x()),
             ya = CGAL::to_double(c1.center().y());
      double xb = CGAL::to_double(c2.center().x()),
             yb = CGAL::to_double(c2.center().y());
      double d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);

      double ra = CGAL::sqrt(CGAL::to_double(c1.squared_radius()));
      double rb = CGAL::sqrt(CGAL::to_double(c2.squared_radius()));
      double K  = CGAL::sqrt(((ra+rb)*(ra+rb)-d2)*(d2-(ra-rb)*(ra-rb)))/4.;

      double xbase = (xb + xa)/2. + (xb - xa)*(ra*ra - rb*rb)/d2/2.;
      double xdiff = 2.*(yb - ya)*K/d2;
      double x1 = xbase + xdiff;
      double x2 = xbase - xdiff;

      double ybase = (yb + ya)/2. + (yb - ya)*(ra*ra - rb*rb)/d2/2.;
      double ydiff = -2.*(xb - xa)*K/d2;
      double y1 = ybase + ydiff;
      double y2 = ybase - ydiff;

      Hyperbolic_point_2 res1(x1, y1);
      Hyperbolic_point_2 res2(x2, y2);
      return std::make_pair(res1, res2);
    }

    Hyperbolic_point_2 operator()(const Hyperbolic_segment_2& s1,
                                  const Hyperbolic_segment_2& s2)
    {
      if(Circular_arc_2* c1 = boost::get<Circular_arc_2>(&s1))
      {
        if(Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2))
        {
          std::pair<Hyperbolic_point_2, Hyperbolic_point_2> res = operator()(c1->circle(),
                                                                             c2->circle());
          Hyperbolic_point_2 p1 = res.first;
          if(p1.x()*p1.x() + p1.y()*p1.y() < FT(1))
            return p1;

          Hyperbolic_point_2 p2 = res.second;
          CGAL_assertion(p2.x()*p2.x() + p2.y()*p2.y() < FT(1));
          return p2;
        }
        else
        {
          Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
          std::pair<Hyperbolic_point_2, Hyperbolic_point_2> res = operator()(c1->circle(),
                                                                             ell2->supporting_line());
          Hyperbolic_point_2 p1 = res.first;
          if(p1.x()*p1.x() + p1.y()*p1.y() < FT(1))
            return p1;

          Hyperbolic_point_2 p2 = res.second;
          CGAL_assertion(p2.x()*p2.x() + p2.y()*p2.y() < FT(1));
          return p2;
        }
      }
      else
      {
        Euclidean_segment_2* ell1 = boost::get<Euclidean_segment_2>(&s1);
        if(Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2))
        {
          std::pair<Hyperbolic_point_2, Hyperbolic_point_2> res = operator()(ell1->supporting_line(),
                                                                             c2->circle());
          Hyperbolic_point_2 p1 = res.first;
          if(p1.x()*p1.x() + p1.y()*p1.y() < FT(1))
            return p1;

          Hyperbolic_point_2 p2 = res.second;
          CGAL_assertion(p2.x()*p2.x() + p2.y()*p2.y()) < FT(1);
          return p2;
        }
        else
        {
          Euclidean_segment_2* ell2 = boost::get<Euclidean_segment_2>(&s2);
          Hyperbolic_point_2 p1 = operator()(ell1->supporting_line(), ell2->supporting_line());
          CGAL_assertion(p1.x()*p1.x() + p1.y()*p1.y()) < FT(1);
          return p1;
        }
      }
    }

  private:
    const Traits& _gt;
  };



  template <typename Traits>
  class Construct_inexact_hyperbolic_circumcenter_2
  {
    typedef typename Traits::FT                             FT;
    typedef typename Traits::Hyperbolic_Voronoi_point_2     Hyperbolic_Voronoi_point_2;
    typedef typename Traits::Euclidean_circle_or_line_2     Euclidean_circle_or_line_2;
    typedef typename Traits::Hyperbolic_point_2             Hyperbolic_point_2;
    typedef typename Traits::Circle_2                       Circle_2;
    typedef typename Traits::Euclidean_line_2               Euclidean_line_2;

  public:
    Construct_inexact_hyperbolic_circumcenter_2(const Traits& gt = Traits()) : _gt(gt) {}

    typedef Hyperbolic_Voronoi_point_2                      result_type;

    Hyperbolic_Voronoi_point_2 operator()(const Hyperbolic_point_2& p,
                                          const Hyperbolic_point_2& q,
                                          const Hyperbolic_point_2& r)
    {
      Origin o;
      Hyperbolic_point_2 po = _gt.construct_hyperbolic_point_2_object()(o);
      Circle_2 l_inf = _gt.construct_circle_2_object()(po, FT(1));

      // Check if |p,O| = |q,O| = |r,O| -- then the circumcenter is the origin O
      if(_gt.compare_distance_2_object()(po,p,q) == EQUAL && _gt.compare_distance_2_object()(po,p,r) == EQUAL)
        return po;

      Euclidean_circle_or_line_2 bis_pq = _gt.construct_circle_or_line_supporting_bisector_2_object()(p,q);
      Euclidean_circle_or_line_2 bis_qr = _gt.construct_circle_or_line_supporting_bisector_2_object()(q,r);

      // now supporting objects cannot both be Euclidean lines

      Euclidean_line_2* l;
      Circle_2* c;

      if(Circle_2* c_pq = boost::get<Circle_2>(&bis_pq))
      {
        if(Circle_2* c_qr = boost::get<Circle_2>(&bis_qr))
        {
          std::pair<Hyperbolic_point_2, Hyperbolic_point_2> inters = _gt.construct_inexact_intersection_2_object()(*c_pq, *c_qr);

          if(_gt.has_on_bounded_side_2_object()(l_inf, inters.first))
            return inters.first;

          return inters.second;
        }
        // here bis_qr is a line
        l = boost::get<Euclidean_line_2>(&bis_qr);
        c = c_pq;
      }
      else
      {
        // here bis_pq is a line
        l = boost::get<Euclidean_line_2>(&bis_pq);
        c = boost::get<Circle_2>(&bis_qr);
      }

      std::pair<Hyperbolic_point_2, Hyperbolic_point_2> inters = _gt.construct_inexact_intersection_2_object()(*c, *l);

      if(_gt.has_on_bounded_side_2_object()(l_inf, inters.first))
        return inters.first;

      return inters.second;
    }

    template <typename Face_handle>
    Hyperbolic_Voronoi_point_2 operator()(const Face_handle fh)
    {
      return operator()(_gt.construct_hyperbolic_point_2_object()(fh.vertex(0)->point(), fh.translation(0)),
                        _gt.construct_hyperbolic_point_2_object()(fh.vertex(1)->point(), fh.translation(1)),
                        _gt.construct_hyperbolic_point_2_object()(fh.vertex(2)->point(), fh.translation(2)) );
    }

  private:
    const Traits& _gt;
  }; // end Construct_inexact_hyperbolic_circumcenter_2



template <typename Traits>
class Side_of_original_octagon
  {
    typedef typename Traits::FT                         FT;
    typedef typename Traits::Circle_2                   Circle_2;
    typedef typename Traits::Hyperbolic_point_2         Hyperbolic_point_2;

  public:
    Side_of_original_octagon(const Traits& gt = Traits()) : _gt(gt) {}

    template <class Point_2_template>
    CGAL::Bounded_side operator()(const Point_2_template& p)
    {
      FT n2(2);
      FT n4(4);
      FT sq2(CGAL::sqrt(n2));
      FT p14(CGAL::sqrt(CGAL::sqrt(n2*n2*n2))); // 2^{3/4}
      FT p34(CGAL::sqrt(CGAL::sqrt(n2*n2*n2))); // 2^{3/4}
      FT s1 (CGAL::sqrt(n2 + CGAL::sqrt(n2))); // CGAL::sqrt(2 + CGAL::sqrt(2))
      FT s2 (CGAL::sqrt(n2 - CGAL::sqrt(n2))); // CGAL::sqrt(2 - CGAL::sqrt(2))
      FT s3 (CGAL::sqrt(sq2 - FT(1))); // CGAL::sqrt(CGAL::sqrt(2) - 1)
      _gt.construct_point_2_object()(s1*p34/n4, -s2*p34/n4);
      Hyperbolic_point_2 V0 = _gt.construct_hyperbolic_point_2_object()(_gt.construct_point_2_object()(s1*p34/n4, -s2*p34/n4));
      Hyperbolic_point_2 V1 = _gt.construct_hyperbolic_point_2_object()(_gt.construct_point_2_object()(p14*(s1+s2)/n4, p14*(s1-s2)/n4));
      Hyperbolic_point_2 V2 = _gt.construct_hyperbolic_point_2_object()(_gt.construct_point_2_object()(s2*p34/n4, s1*p34/n4));

      Hyperbolic_point_2 M0 = _gt.construct_hyperbolic_point_2_object()(_gt.construct_point_2_object()(s3, FT(0)));
      Hyperbolic_point_2 M1 = _gt.construct_hyperbolic_point_2_object()(_gt.construct_point_2_object()(sq2*s3/n2, sq2*s3/n2));

      // Poincare disk (i.e., unit Euclidean disk)
      Circle_2 Poincare = _gt.construct_circle_2_object()(_gt.construct_hyperbolic_point_2_object()(_gt.construct_point_2_object()(0, 0)), 1*1);

      // Euclidean circle supporting the side [V0, V1]
      Circle_2 SuppCircle1 = _gt.construct_circle_2_object()(M0, V0, V1);

      // Euclidean circle supporting the side [V1, V2]
      Circle_2 SuppCircle2 = _gt.construct_circle_2_object()(M1, V1, V2);

      // This transformation brings the point in the first quadrant (positive x, positive y)
      FT x(FT(p.x()) > FT(0) ? p.x() : -p.x());
      FT y(FT(p.y()) > FT(0) ? p.y() : -p.y());

      // This brings the point in the first octant (positive x and y < x)
      if(y > x)
      {
        FT tmp = x;
        x = y;
        y = tmp;
      }

      // This tells us whether the point is on the side of the open boundary
      bool on_open_side = ((p.y() + std::tan(CGAL_PI / 8.) * p.x()) < 0.0);

      Hyperbolic_point_2 t = _gt.construct_hyperbolic_point_2_object()(x, y);

      CGAL::Bounded_side PoincareSide = _gt.bounded_side_2_object()(Poincare, t);
      CGAL::Bounded_side Circ1Side = _gt.bounded_side_2_object()(SuppCircle1, t);
      CGAL::Bounded_side Circ2Side = _gt.bounded_side_2_object()(SuppCircle2, t);

      // First off, the point needs to be inside the Poincare disk. if not, there's no hope.
      if(PoincareSide == CGAL::ON_BOUNDED_SIDE)
      {
        // Inside the Poincare disk, but still outside the original domain
        if(Circ1Side  == CGAL::ON_BOUNDED_SIDE || Circ2Side == CGAL::ON_BOUNDED_SIDE)
          return CGAL::ON_UNBOUNDED_SIDE;

        // Inside the Poincare disk and inside the original domain
        if(Circ1Side  == CGAL::ON_UNBOUNDED_SIDE && Circ2Side == CGAL::ON_UNBOUNDED_SIDE)
          return CGAL::ON_BOUNDED_SIDE;

        // This is boundary, but we only consider the upper half. The lower half means outside.
        if(on_open_side)
          return CGAL::ON_UNBOUNDED_SIDE;
        else
          return CGAL::ON_BOUNDED_SIDE;
      } else
      {
        return CGAL::ON_UNBOUNDED_SIDE;
      }
    }

  private:
    const Traits& _gt;
  };


} // end namespace internal


template <typename Kernel = Exact_predicates_exact_constructions_kernel_with_sqrt,
          template<typename> class Translation_type = Hyperbolic_octagon_translation>
class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2
  : public Hyperbolic_Delaunay_triangulation_traits_2<Kernel>
{
  typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>                              Base;
  typedef Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel, Translation_type> Self;

public:
  typedef typename Base::FT                                               FT;
  typedef Translation_type<FT>                                            Hyperbolic_translation;

  typedef typename Base::Hyperbolic_point_2                               Hyperbolic_point_2;
  typedef Hyperbolic_point_2                                              Hyperbolic_Voronoi_point_2;
  typedef typename Base::Line_2                                           Euclidean_line_2;
  typedef typename Base::Euclidean_circle_or_line_2                       Euclidean_circle_or_line_2;
  typedef typename Base::Circular_arc_2                                   Circular_arc_2;
  typedef typename Base::Euclidean_segment_2                              Euclidean_segment_2; //only used internally here
  typedef typename Base::Hyperbolic_segment_2                             Hyperbolic_segment_2;
  typedef Euclidean_segment_2                                             Line_segment_2; // kept for demo
  typedef Hyperbolic_segment_2                                            Segment_2; // kept for demo

  // the following types are only used internally in this traits class,
  // so they need not be documented, and they don't need _object()
  typedef typename Base::Construct_Euclidean_bisector_2                   Construct_Euclidean_bisector_2;
  typedef typename Base::Construct_intersection_2                         Construct_intersection_2;
  typedef typename Base::Construct_hyperbolic_bisector_2                  Construct_hyperbolic_bisector_2;
  typedef typename Base::Construct_circle_or_line_supporting_bisector     Construct_circle_or_line_supporting_bisector;
  typedef typename Base::Euclidean_collinear_2                            Euclidean_collinear_2;
  typedef typename Base::Compute_squared_Euclidean_distance_2             Compute_squared_Euclidean_distance_2;
  typedef typename Base::Has_on_bounded_side_2                            Has_on_bounded_side_2;

  // Wrappers for the translation adapter
  typedef internal::Hyperbolic_traits_with_translations_2_adaptor<
            Self, typename Base::Orientation_2>                           Orientation_2;

  typedef internal::Hyperbolic_traits_with_translations_2_adaptor<
            Self, typename Base::Side_of_oriented_circle_2>               Side_of_oriented_circle_2;

  typedef internal::Hyperbolic_traits_with_translations_2_adaptor<
            Self, typename Base::Construct_hyperbolic_circumcenter_2>     Construct_hyperbolic_circumcenter_2;

  typedef internal::Hyperbolic_traits_with_translations_2_adaptor<
            Self, typename Base::Construct_hyperbolic_segment_2>          Construct_hyperbolic_segment_2;

  typedef internal::Hyperbolic_traits_with_translations_2_adaptor<
            Self, typename Base::Construct_segment_2>                     Construct_segment_2;

  typedef internal::Hyperbolic_traits_with_translations_2_adaptor<
            Self, typename Base::Construct_triangle_2>                    Construct_triangle_2;

  typedef internal::Hyperbolic_traits_with_translations_2_adaptor<
            Self, typename Base::Compare_distance_2>                      Compare_distance_2;

  typedef internal::Periodic_4_construct_hyperbolic_point_2<
            Self, typename Kernel::Construct_point_2>                     Construct_hyperbolic_point_2;

  typedef internal::Hyperbolic_traits_with_translations_2_adaptor<
            Self, typename Base::Side_of_oriented_hyperbolic_segment_2>   Side_of_oriented_hyperbolic_segment_2;


  typedef internal::Compute_approximate_hyperbolic_diameter<Self>         Compute_approximate_hyperbolic_diameter;
  typedef internal::Construct_hyperbolic_line_2<Self>                     Construct_hyperbolic_line_2;
  typedef internal::Construct_inexact_intersection_2<Self>                Construct_inexact_intersection_2;
  typedef internal::Construct_inexact_hyperbolic_circumcenter_2<Self>     Construct_inexact_hyperbolic_circumcenter_2;
  typedef internal::Side_of_original_octagon<Self>                        Side_of_original_octagon;


public:

  Compute_approximate_hyperbolic_diameter
  compute_approximate_hyperbolic_diameter_object() const
  { return Compute_approximate_hyperbolic_diameter(*this); }


  Construct_hyperbolic_line_2
  construct_hyperbolic_line_2_object() const
  { return Construct_hyperbolic_line_2(*this); }


  Construct_inexact_intersection_2
  construct_inexact_intersection_2_object() const
  { return Construct_inexact_intersection_2(*this); }


  Construct_inexact_hyperbolic_circumcenter_2
  construct_inexact_hyperbolic_circumcenter_2_object() const
  { return Construct_inexact_hyperbolic_circumcenter_2(*this); }


  Side_of_original_octagon
  side_of_original_octagon_object() const
  { return Side_of_original_octagon(*this); }


  Side_of_oriented_hyperbolic_segment_2
  side_of_oriented_hyperbolic_segment_2_object() const
  { return Side_of_oriented_hyperbolic_segment_2(*this); }


  Construct_triangle_2
  construct_triangle_2_object() const
  { return Construct_triangle_2(); }

  Construct_hyperbolic_point_2
  construct_hyperbolic_point_2_object() const
  { return Construct_hyperbolic_point_2(); }

  Construct_hyperbolic_segment_2
  construct_hyperbolic_segment_2_object() const
  { return Construct_hyperbolic_segment_2(); }

  Construct_segment_2
  construct_segment_2_object() const
  { return Construct_segment_2(); }

  Orientation_2
  orientation_2_object() const
  { return Orientation_2(); }

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  { return Side_of_oriented_circle_2(); }


public:
  Periodic_4_hyperbolic_Delaunay_triangulation_traits_2(const Base& base = Base()) : Base(base) {}


}; // class Periodic_4_hyperbolic_Delaunay_triangulation_traits_2

} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
