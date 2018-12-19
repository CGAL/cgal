// Copyright (c) 2016  INRIA Nancy - Grand Est (France).
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
// $URL: $
// $Id:  $
//
// Author(s)     : Iordan Iordanov
//                 Monique Teillaud
//

#ifndef CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
#define CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/determinant.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/internal/Exact_complex.h>
#include <CGAL/Origin.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/utility.h>

#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"

using std::pair;
using std::make_pair;

namespace CGAL {

template<typename Kernel = CGAL::Cartesian<CORE::Expr> >
class Hyperbolic_Delaunay_triangulation_traits_2
    : public Kernel
{
  typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>  Self;

private:
  class Circular_arc_2_base
  {
    typedef typename Kernel::FT                               FT;
    typedef Exact_complex<FT>                                 Cplx;
    typedef typename Kernel::Point_2                          Point;
    typedef typename Kernel::Circle_2                         Circle;
    typedef typename Kernel::Orientation_2                    Orientation_2;

  private:
    Circle _c;
    Point _s, _t;

  public:
    Circular_arc_2_base()
      : _c(Point(FT(0),FT(0)), FT(0)), _s(FT(0),FT(0)), _t(FT(0),FT(0))
    {}

    Circular_arc_2_base(Circle c, Point source, Point target)
      : _c(c), _s(source), _t(target)
    {}

    Circular_arc_2_base(Point p1, Point p2)
    {
      Cplx p(p1), q(p2);
      Cplx O(0,0);
      Cplx inv;
      if(p == O)
        inv = q.invert_in_unit_circle();
      else
        inv = p.invert_in_unit_circle();

      Point ip(inv.real(), inv.imag());

      _c = Circle(p1, p2, ip);
      if(Orientation_2()(p1, p2, _c.center()) == LEFT_TURN)
      {
        _s = p1;
        _t = p2;
      }
      else
      {
        _s = p2;
        _t = p1;
      }

    }

    Circle supporting_circle() const { return _c; }
    Point source() const { return _s; }
    Point target() const { return _t; }
    FT squared_radius() const { return _c.squared_radius(); }
    Point center() const { return _c.center(); }
    Bbox_2 bbox(void) const { return typename Kernel::Construct_bbox_2()(*this); }
  };

public:
  typedef typename Kernel::FT                                 FT;
  typedef typename Kernel::RT                                 RT;
  typedef typename Kernel::Kernel_base                        Kernel_base;
  typedef typename Kernel::Point_2                            Hyperbolic_point_2;
  typedef Hyperbolic_point_2                                  Hyperbolic_Voronoi_point_2;
  typedef typename Kernel::Circle_2                           Circle_2;
  typedef typename Kernel::Line_2                             Euclidean_line_2;
  typedef boost::variant<Circle_2,Euclidean_line_2>           Euclidean_circle_or_line_2;
  typedef Self::Circular_arc_2_base                           Circular_arc_2;
  typedef typename Kernel::Segment_2                          Euclidean_segment_2; // only used internally here
  typedef boost::variant<Circular_arc_2, Euclidean_segment_2> Hyperbolic_segment_2;

  typedef typename Kernel::Triangle_2                         Hyperbolic_triangle_2;
  typedef typename Kernel::Orientation_2                      Orientation_2;
  typedef typename Kernel::Side_of_oriented_circle_2          Side_of_oriented_circle_2;

  // only kept for demo to please T2graphicsitems
  typedef Euclidean_segment_2                                 Line_segment_2;
  typedef Hyperbolic_segment_2                                Segment_2;

  typedef typename Kernel::Compare_x_2                        Compare_x_2;
  typedef typename Kernel::Compare_y_2                        Compare_y_2;

  typedef typename Kernel::Less_x_2                           Less_x_2;
  typedef typename Kernel::Less_y_2                           Less_y_2;

  // The objects Ray_2, Iso_rectangle_2 and Line_2 are needed by the CGAL::Qt::PainterOstream
  typedef typename Kernel::Direction_2                        Direction_2;
  typedef typename Kernel::Vector_2                           Vector_2;
  typedef typename Kernel::Ray_2                              Ray_2;
  typedef typename Kernel::Iso_rectangle_2                    Iso_rectangle_2;
  typedef Euclidean_line_2                                    Line_2;

  // the following types are only used internally in this traits class,
  // so they need not be documented, and they don't need _object()
  typedef typename Kernel::Collinear_2                        Euclidean_collinear_2;
  typedef typename Kernel::Construct_bisector_2               Construct_Euclidean_bisector_2;
  typedef typename Kernel::Construct_midpoint_2               Construct_Euclidean_midpoint_2;
  typedef typename Kernel::Construct_triangle_2               Construct_triangle_2;
  typedef typename Kernel::Compare_distance_2                 Compare_distance_2;
  typedef typename Kernel::Has_on_bounded_side_2              Has_on_bounded_side_2;
  typedef typename Kernel::Compute_squared_distance_2         Compute_squared_Euclidean_distance_2;

  // Can/should we keep those objects public?
public:
  class Construct_hyperbolic_segment_2
  {
    typedef typename Kernel::Construct_weighted_circumcenter_2 Construct_weighted_circumcenter_2;
    typedef typename Kernel::Weighted_point_2                  Weighted_point_2;
    typedef typename Kernel::Point_2                           Bare_point;

  public:
    typedef Hyperbolic_segment_2 result_type;

    Construct_hyperbolic_segment_2() {}

    Hyperbolic_segment_2 operator()(const Hyperbolic_point_2& p, const Hyperbolic_point_2& q) const
    {
      Origin o;
      if(Euclidean_collinear_2()(p, q, Hyperbolic_point_2(o)))
        return Euclidean_segment_2(p, q);

      Weighted_point_2 wp(p);
      Weighted_point_2 wq(q);
      Weighted_point_2 wo(Hyperbolic_point_2(o), FT(1)); // Poincaré circle

      Bare_point center = Construct_weighted_circumcenter_2()(wp, wo, wq);
      FT sq_radius = Compute_squared_Euclidean_distance_2()(p, center);

      Circle_2 circle(center, sq_radius);
      // uncomment!!!
      //assert(circle.has_on_boundary(p) && circle.has_on_boundary(q));

      if(Orientation_2()(p, q, center) == LEFT_TURN)
        return Circular_arc_2(circle, p, q);

      return Circular_arc_2(circle, q, p);
    }

  }; // end Construct_hyperbolic_segment_2

  Construct_hyperbolic_segment_2
  construct_hyperbolic_segment_2_object() const
  { return Construct_hyperbolic_segment_2(); }

  // wrong names kept for demo
  typedef Construct_hyperbolic_segment_2 Construct_segment_2;
  Construct_segment_2
  construct_segment_2_object() const
  { return Construct_hyperbolic_segment_2(); }

  class Construct_hyperbolic_circumcenter_2
  {
  public:

    typedef Hyperbolic_Voronoi_point_2                    result_type;

    Hyperbolic_Voronoi_point_2 operator()(Hyperbolic_point_2 p,
                                          Hyperbolic_point_2 q,
                                          Hyperbolic_point_2 r)
    {
      Origin o;
      Hyperbolic_point_2 po = Hyperbolic_point_2(o);
      Circle_2 l_inf(po, FT(1));

      Euclidean_circle_or_line_2 bis_pq = Construct_circle_or_line_supporting_bisector()(p, q);
      Euclidean_circle_or_line_2 bis_qr = Construct_circle_or_line_supporting_bisector()(q, r);

      if(Compare_distance_2()(po, p, q) == EQUAL && Compare_distance_2()(po, p,r) == EQUAL)
        return po;

      // now supporting objects cannot both be Euclidean lines
      Euclidean_line_2* l;
      Circle_2* c;

      if(Circle_2* c_pq = boost::get<Circle_2>(&bis_pq))
      {
        if(Circle_2* c_qr = boost::get<Circle_2>(&bis_qr))
        {
          std::pair<Hyperbolic_point_2, Hyperbolic_point_2> inters = Construct_intersection_2()(*c_pq, *c_qr);

          if(Has_on_bounded_side_2()(l_inf, inters.first))
            return inters.first;

          return inters.second;
        }

        l = boost::get<Euclidean_line_2>(&bis_qr);
        c = c_pq;
      }
      else
      {
        // here bis_pq is a line
        l = boost::get<Euclidean_line_2>(&bis_pq);
        c = boost::get<Circle_2>(&bis_qr);
      }

      std::pair<Hyperbolic_point_2, Hyperbolic_point_2> inters = Construct_intersection_2()(*c, *l);

      if(Has_on_bounded_side_2()(l_inf, inters.first))
        return inters.first;

      return inters.second;
    }
  }; // end Construct_hyperbolic_circumcenter_2

  Compare_x_2 compare_x_2_object() const { return Compare_x_2(); }
  Compare_y_2 compare_y_2_object() const { return Compare_y_2(); }
  Orientation_2 orientation_2_object() const { return Orientation_2(); }

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  { return Side_of_oriented_circle_2(); }

  Construct_hyperbolic_circumcenter_2
  construct_hyperbolic_circumcenter_2_object() const
  { return Construct_hyperbolic_circumcenter_2(); }

  class Construct_hyperbolic_bisector_2
  {
  public:
    Construct_hyperbolic_bisector_2() {}

    // constructs a hyperbolic line
    Hyperbolic_segment_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q) const
    {
      Origin o;
      Hyperbolic_point_2 po = Hyperbolic_point_2(o);
      Circle_2 l_inf = Circle_2(po,FT(1));

      if(Compare_distance_2()(po, p, q) == EQUAL)
      {
        Euclidean_line_2 l = Construct_Euclidean_bisector_2()(p, q);
        std::pair<Hyperbolic_point_2, Hyperbolic_point_2> inters = Construct_intersection_2()(l, l_inf);
        return Euclidean_segment_2(inters.first, inters.second);
      }

      Euclidean_circle_or_line_2 bis_pq = Construct_circle_or_line_supporting_bisector()(p, q);
      Circle_2* c = boost::get<Circle_2>(&bis_pq);
      std::pair<Hyperbolic_point_2, Hyperbolic_point_2> inters = Construct_intersection_2()(*c, l_inf);

      if(Orientation_2()(c->center(), inters.first, inters.second) == POSITIVE)
        return Circular_arc_2(*c, inters.first, inters.second);
      else
        return Circular_arc_2(*c, inters.second, inters.first);
    }

    // constructs the hyperbolic bisector of segment [p, q] limited by
    // circumcenter(p, q, r) on one side
    // and circumcenter(p, s, q) on the other side
    Hyperbolic_segment_2 operator()(Hyperbolic_point_2 p,
                                    Hyperbolic_point_2 q,
                                    Hyperbolic_point_2 r,
                                    Hyperbolic_point_2 s)
    {
      CGAL_triangulation_precondition((Orientation_2()(p, q, r) == ON_POSITIVE_SIDE) &&
                                      (Orientation_2()(p, s, q) == ON_POSITIVE_SIDE));
      CGAL_triangulation_precondition((Side_of_oriented_circle_2()(p, q, r,s) == ON_NEGATIVE_SIDE) &&
                                      (Side_of_oriented_circle_2()(p, s, q, r) == ON_NEGATIVE_SIDE));

      Origin o;
      Hyperbolic_point_2 po = Hyperbolic_point_2(o);

      // TODO MT this is non-optimal...
      // the bisector is already computed here
      // and it will be recomputed below
      Hyperbolic_point_2 a = Construct_hyperbolic_circumcenter_2()(p, q, r);
      Hyperbolic_point_2 b = Construct_hyperbolic_circumcenter_2()(p, s, q);

      if(Compare_distance_2()(po, p, q) == EQUAL)
      {
        Euclidean_line_2 l = Construct_Euclidean_bisector_2()(p, q);
        return Euclidean_segment_2(a, b);
      }

      Euclidean_circle_or_line_2 bis_pq = Construct_circle_or_line_supporting_bisector()(p, q);
      Circle_2* c_pq = boost::get<Circle_2>(&bis_pq);

      if(Compare_distance_2()(po, p, q) == POSITIVE) {
        // then p is inside the supporting circle
        return Circular_arc_2(*c_pq, b, a);
      }

      return Circular_arc_2(*c_pq, a, b);
    }

    // constructs the hyperbolic bisector of segment [p, q]
    // limited by hyperbolic circumcenter(p, q, r) on one side
    // and going to the infinite line on the other side
    Hyperbolic_segment_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 r)
    {
      CGAL_triangulation_precondition(Orientation_2()(p, q, r) == POSITIVE);

      Origin o;
      Hyperbolic_point_2 po = Hyperbolic_point_2(o);
      Circle_2 l_inf(po, FT(1));

      // TODO MT this is non-optimal...
      // the bisector is computed (at least) twice
      Hyperbolic_point_2 a = Construct_hyperbolic_circumcenter_2()(p, q, r);

      if(Compare_distance_2()(po, p, q) == EQUAL)
      {
        Euclidean_line_2 bis_pq = Construct_Euclidean_bisector_2()(p, q);

        std::pair<Hyperbolic_point_2,Hyperbolic_point_2> inters = Construct_intersection_2()(bis_pq, l_inf);
        if(Less_y_2()(p, q))
          return Euclidean_segment_2(a,inters.first);

        return Euclidean_segment_2(a,inters.second);
      }

      Euclidean_circle_or_line_2
          bis_pq = Construct_circle_or_line_supporting_bisector()(p, q);
      Circle_2* c_pq = boost::get<Circle_2>(&bis_pq);

      std::pair<Hyperbolic_point_2,Hyperbolic_point_2> inters = Construct_intersection_2()(*c_pq, l_inf);
      Hyperbolic_point_2 approx_pinf = inters.first;
      if(Orientation_2()(p, q,inters.first) == NEGATIVE)
      {
        if(Orientation_2()(c_pq->center(),a,inters.first) == POSITIVE)
          return Circular_arc_2(*c_pq, a, inters.first);

        return Circular_arc_2(*c_pq, inters.first, a);
      }

      if(Orientation_2()(c_pq->center(),a,inters.first) == POSITIVE)
        return Circular_arc_2(*c_pq, inters.second, a);

      return Circular_arc_2(*c_pq, a, inters.second);
    }
  }; // end Construct_hyperbolic_bisector_2

  Construct_hyperbolic_bisector_2
  construct_hyperbolic_bisector_2_object() const
  { return Construct_hyperbolic_bisector_2(); }

  Construct_Euclidean_bisector_2
  construct_Euclidean_bisector_2_object() const
  { return Construct_Euclidean_bisector_2(); }

  // For details see the JoCG paper (5:56-85, 2014)
  class Is_Delaunay_hyperbolic
  {
  public:
    typedef typename Kernel::Vector_3    Vector_3;
    typedef typename Kernel::Point_3     Point_3;

    bool operator()(const Hyperbolic_point_2& p0,
                    const Hyperbolic_point_2& p1,
                    const Hyperbolic_point_2& p2) const
    {
      Vector_3 v0 = Vector_3(p0.x()*p0.x() + p0.y()*p0.y(),
                             p1.x()*p1.x() + p1.y()*p1.y(),
                             p2.x()*p2.x() + p2.y()*p2.y());

      Vector_3 v1 = Vector_3(p0.x(), p1.x(), p2.x());
      Vector_3 v2 = Vector_3(p0.y(), p1.y(), p2.y());
      Vector_3 v3 = Vector_3(FT(1), FT(1), FT(1));

      FT dt0 = determinant(v0, v1, v3);
      FT dt1 = determinant(v0, v2, v3);
      FT dt2 = determinant(v0 - v3, v1, v2);

      return (dt0*dt0 + dt1*dt1 - dt2*dt2 < 0);
    }

    bool operator()(const Hyperbolic_point_2& p0,
                    const Hyperbolic_point_2& p1,
                    const Hyperbolic_point_2& p2,
                    int& ind) const
    {
      if(this->operator()(p0, p1, p2) == false)
      {
        ind = find_non_hyperbolic_edge(p0, p1, p2);
        return false;
      }

      return true;
    }

  private:

    // assume the face (p0, p1, p2) is non-hyperbolic
    int find_non_hyperbolic_edge(const Hyperbolic_point_2& p0,
                                 const Hyperbolic_point_2& p1,
                                 const Hyperbolic_point_2& p2) const
    {
      typedef typename Kernel::Direction_2 Direction_2;

      Vector_3 v0 = Vector_3(p0.x()*p0.x() + p0.y()*p0.y(),
                             p1.x()*p1.x() + p1.y()*p1.y(),
                             p2.x()*p2.x() + p2.y()*p2.y());

      Vector_3 v1 = Vector_3(p0.x(), p1.x(), p2.x());
      Vector_3 v2 = Vector_3(p0.y(), p1.y(), p2.y());
      Vector_3 v3 = Vector_3(FT(1), FT(1), FT(1));

      FT dt0 = determinant(v0, 2*v2, -v3);
      FT dt1 = determinant(2*v1, v0, -v3);
      FT dt2 = determinant(2*v1, 2*v2, -v3);

      Direction_2 d0(p0.x()*dt2 - dt0, p0.y()*dt2 - dt1);
      Direction_2 d1(p1.x()*dt2 - dt0, p1.y()*dt2 - dt1);
      Direction_2 d2(p2.x()*dt2 - dt0, p2.y()*dt2 - dt1);

      Direction_2 d(dt0, dt1);

      if(d.counterclockwise_in_between(d0, d1))
        return 2;

      if(d.counterclockwise_in_between(d1, d2))
        return 0;

      return 1;
    }
  }; // end Is_Delaunay_hyperbolic

  Is_Delaunay_hyperbolic
  is_Delaunay_hyperbolic_object() const
  { return Is_Delaunay_hyperbolic(); }

  // do not document
  // constructs the Euclidean circle or line supporting the hyperbolic
  // bisector of two points
  class Construct_circle_or_line_supporting_bisector
  {
  public:
    Construct_circle_or_line_supporting_bisector()
    {}

    Euclidean_circle_or_line_2 operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q) const
    {
      Origin o;
      Hyperbolic_point_2 po = Hyperbolic_point_2(o);
      typedef typename Kernel::Point_3     Point_3;

      if(Compare_distance_2()(po, p, q) == EQUAL)
        return Construct_Euclidean_bisector_2()(p, q);

      FT dop2 = p.x()*p.x() + p.y()*p.y();
      FT doq2 = q.x()*q.x() + q.y()*q.y();
      Point_3 p3(p.x(), p.y(), dop2);
      Point_3 q3(q.x(), q.y(), doq2);

      // TODO MT improve

      // The cirle belongs to the pencil with limit points p and q
      // p, q are zero-circles
      // (x, y, xˆ2 + yˆ2 - rˆ2) = alpha*(xp, yp, xpˆ2 + ypˆ2) + (1-alpha)*(xq, yq, xqˆ2 + yqˆ2)
      // xˆ2 + yˆ2 - rˆ2 = 1 (= radius of the Poincare disc)
      FT op = p.x()*p.x() + p.y()*p.y();
      FT oq = q.x()*q.x() + q.y()*q.y();
      FT alpha = (FT(1) - oq) / (op - oq);

      FT x = alpha*p.x() + (1-alpha)*q.x();
      FT y = alpha*p.y() + (1-alpha)*q.y();
      FT sq_radius = x*x + y*y - FT(1);

      // TODO MT improve
      // ?? orientation should depend on
      // Compare_distance(O,p, q)
      // so that p always on positive side
      // ???
      // CK does not care about orientation, circular arcs are
      // considered in CCW order in any case

      Euclidean_line_2 l = Construct_Euclidean_bisector_2()(p, q);
      Hyperbolic_point_2 middle = Construct_Euclidean_midpoint_2()(p, q);
      Hyperbolic_point_2 temp = middle + l.to_vector();

      if(Orientation_2()(middle, temp, Hyperbolic_point_2(x, y)) == ON_POSITIVE_SIDE)
        return Circle_2(Hyperbolic_point_2(x, y), sq_radius, CLOCKWISE);

      return Circle_2(Hyperbolic_point_2(x, y), sq_radius, COUNTERCLOCKWISE);
    }
  }; // end Construct_supporting_circle_of_bisector

  class Construct_intersection_2
  {
  public:
    Construct_intersection_2() {}

    Hyperbolic_point_2 operator()(Euclidean_line_2 ell1, Euclidean_line_2 ell2)
    {
      // The only point where two Euclidean lines can intersect in the Poincaré disk is the origin.
      return Hyperbolic_point_2(0,0);
    }

    std::pair<Hyperbolic_point_2, Hyperbolic_point_2> operator()(Euclidean_line_2 ell, Circle_2 c)
    {
      if(ell.b() == FT(0))
      {
        FT p = c.center().x();
        FT q = c.center().y();
        FT y1 = q + CGAL::sqrt(c.squared_radius() - p*p);
        FT y2 = q - CGAL::sqrt(c.squared_radius() - p*p);
        Hyperbolic_point_2 p1(FT(0), y1);
        Hyperbolic_point_2 p2(FT(0), y2);
        return make_pair(p1, p2);
      }

      FT lambda = -ell.a()/ell.b();
      FT mu = -ell.c()/ell.b();
      FT p = c.center().x();
      FT q = c.center().y();
      FT A = FT(1) + lambda*lambda;
      FT B = FT(2)*(lambda * mu - lambda*q - p);
      FT C = p*p + mu*mu + q*q - c.squared_radius() - FT(2)*q*mu;
      FT Delta = B*B - FT(4)*A*C;
      FT x1 = (-B + CGAL::sqrt(Delta))/(FT(2)*A);
      FT x2 = (-B - CGAL::sqrt(Delta))/(FT(2)*A);
      FT y1 = lambda*x1 + mu;
      FT y2 = lambda*x2 + mu;
      Hyperbolic_point_2 sol1(x1, y1);
      Hyperbolic_point_2 sol2(x2, y2);
      return make_pair(sol1, sol2);
    }

    std::pair<Hyperbolic_point_2, Hyperbolic_point_2> operator()(Circle_2 c, Euclidean_line_2 ell)
    {
      return operator()(ell, c);
    }

    std::pair<Hyperbolic_point_2, Hyperbolic_point_2> operator()(Circle_2 c1, Circle_2 c2)
    {
      FT xa = c1.center().x(), ya = c1.center().y();
      FT xb = c2.center().x(), yb = c2.center().y();
      FT d2 = squared_distance(c1.center(), c2.center());
      FT ra = CGAL::sqrt(c1.squared_radius());
      FT rb = CGAL::sqrt(c2.squared_radius());
      FT K  = CGAL::sqrt(((ra+rb)*(ra+rb)-d2)*(d2-(ra-rb)*(ra-rb)))/FT(4);

      FT xbase = (xb + xa)/FT(2) + (xb - xa)*(ra*ra - rb*rb)/d2/FT(2);
      FT xdiff = FT(2)*(yb - ya)*K/d2;
      FT x1 = xbase + xdiff;
      FT x2 = xbase - xdiff;

      FT ybase = (yb + ya)/FT(2) + (yb - ya)*(ra*ra - rb*rb)/d2/FT(2);
      FT ydiff = FT(-2)*(xb - xa)*K/d2;
      FT y1 = ybase + ydiff;
      FT y2 = ybase - ydiff;

      Hyperbolic_point_2 res1(x1, y1);
      Hyperbolic_point_2 res2(x2, y2);
      return make_pair(res1, res2);
    }

    Hyperbolic_point_2 operator()(Hyperbolic_segment_2 s1, Hyperbolic_segment_2 s2)
    {
      if(Circular_arc_2* c1 = boost::get<Circular_arc_2>(&s1))
      {
        if(Circular_arc_2* c2 = boost::get<Circular_arc_2>(&s2))
        {
          pair<Hyperbolic_point_2, Hyperbolic_point_2> res = operator()(c1->circle(), c2->circle());
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
          pair<Hyperbolic_point_2, Hyperbolic_point_2> res = operator()(c1->circle(), ell2->supporting_line());
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
          pair<Hyperbolic_point_2, Hyperbolic_point_2> res = operator()(ell1->supporting_line(), c2->circle());
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
          Hyperbolic_point_2 p1 = operator()(ell1->supporting_line(), ell2->supporting_line());
          CGAL_assertion(p1.x()*p1.x() + p1.y()*p1.y() < FT(1));
          return p1;
        }
      }
    }
  };

  Construct_intersection_2
  construct_intersection_2_object() const
  { return Construct_intersection_2(); }

  class Side_of_oriented_hyperbolic_segment_2
  {
    typedef typename Kernel::Construct_weighted_circumcenter_2  Construct_weighted_circumcenter_2;
    typedef typename Kernel::Weighted_point_2                   Weighted_point_2;
    typedef typename Kernel::Point_2                            Bare_point;

  public:
    Side_of_oriented_hyperbolic_segment_2() {}

    typedef Oriented_side                                       result_type;

    result_type operator()(Hyperbolic_point_2 p, Hyperbolic_point_2 q, Hyperbolic_point_2 query) const
    {
      // Check first if the points are collinear with the origin
      Circle_2 poincare(Hyperbolic_point_2(FT(0),FT(0)), FT(1));
      Hyperbolic_point_2 O(FT(0), FT(0));
      Orientation ori = orientation(p, q, O);
      if(ori == COLLINEAR)
      {
        Euclidean_line_2 seg(p, q);
        Orientation qori = orientation(p, q, query);
        if(qori == COLLINEAR)
        {
          return ON_ORIENTED_BOUNDARY;
        }
        else
        {
          // It is sufficient that these are consistent.
          if(qori == LEFT_TURN)
            return ON_POSITIVE_SIDE;
          else
            return ON_NEGATIVE_SIDE;
        }
      }

      Weighted_point_2 wp(p);
      Weighted_point_2 wq(q);
      Weighted_point_2 wo(O, FT(1)); // Poincaré circle

      Bare_point center = Construct_weighted_circumcenter_2()(wp, wo, wq);
      FT sq_radius = Compute_squared_Euclidean_distance_2()(p, center);

      Circle_2 circle(center, sq_radius);
      Bounded_side bs = circle.bounded_side(query);
      if(bs == ON_BOUNDARY)
      {
        return ON_ORIENTED_BOUNDARY;
      }
      else
      {
        if(bs == ON_BOUNDED_SIDE)
          return ON_POSITIVE_SIDE;
        else
          return ON_NEGATIVE_SIDE;
      }
    }
  };

  Side_of_oriented_hyperbolic_segment_2
  side_of_oriented_hyperbolic_segment_2_object() const
  { return Side_of_oriented_hyperbolic_segment_2(); }

public:
  Hyperbolic_Delaunay_triangulation_traits_2() {}
  Hyperbolic_Delaunay_triangulation_traits_2(const Hyperbolic_Delaunay_triangulation_traits_2 & other) {}

  Hyperbolic_Delaunay_triangulation_traits_2 &operator=(const Hyperbolic_Delaunay_triangulation_traits_2 &)
  {
    return *this;
  }

}; // class Hyperbolic_Delaunay_triangulation_traits_2

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
