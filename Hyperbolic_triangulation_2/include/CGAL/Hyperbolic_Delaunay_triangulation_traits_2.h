// Copyright (c) 2016-2018  INRIA Nancy - Grand Est (France).
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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
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

namespace internal {

template <typename Traits>
class HDT_2_Circular_arc_2
{
  typedef typename Traits::FT                               FT;
  typedef Exact_complex<FT>                                 Cplx;
  typedef typename Traits::Point_2                          Point;
  typedef typename Traits::Circle_2                         Circle;
  typedef typename Traits::Orientation_2                    Orientation_2;

private:
  Circle _c;
  Point _s, _t;
  const Traits& _gt;

public:
  HDT_2_Circular_arc_2(const Traits& gt = Traits())
    : _c(Point(FT(0),FT(0)), FT(0)),
      _s(FT(0),FT(0)),
      _t(FT(0),FT(0)),
      _gt(gt)
  {}

  HDT_2_Circular_arc_2(const Circle& c,
                       const Point& source, const Point& target,
                       const Traits& gt = Traits())
    : _c(c), _s(source), _t(target), _gt(gt)
  {}

  HDT_2_Circular_arc_2(const Point& p1, const Point& p2,
                       const Traits& gt = Traits())
    : _gt(gt)
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
    if(_gt.orientation_2_object()(p1, p2, _c.center()) == LEFT_TURN)
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
  Bbox_2 bbox(void) const { return _gt.construct_bbox_2_object()(*this); }
};

// only used internally to construct the Euclidean circle or line supporting the hyperbolic
// bisector of two points
template <typename Traits>
class Construct_circle_or_line_supporting_bisector
{
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Hyperbolic_point_2                 Hyperbolic_point_2;
  typedef typename Traits::Euclidean_line_2                   Euclidean_line_2;
  typedef typename Traits::Euclidean_circle_or_line_2         Euclidean_circle_or_line_2;
  typedef typename Traits::Circle_2                           Circle_2;
  typedef typename Traits::Point_3                            Point_3;

public:
  Construct_circle_or_line_supporting_bisector(const Traits& gt = Traits()) : _gt(gt) {}

  Euclidean_circle_or_line_2 operator()(const Hyperbolic_point_2& p,
                                        const Hyperbolic_point_2& q) const
  {
    Hyperbolic_point_2 po(CGAL::ORIGIN);

    if(_gt.compare_distance_2_object()(po, p, q) == EQUAL)
      return _gt.construct_bisector_2_object()(p, q);

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

    Euclidean_line_2 l = _gt.construct_bisector_2_object()(p, q);
    Hyperbolic_point_2 middle = _gt.construct_midpoint_2_object()(p, q);
    Hyperbolic_point_2 temp = middle + l.to_vector();

    if(_gt.orientation_2_object()(middle, temp, Hyperbolic_point_2(x, y)) == ON_POSITIVE_SIDE)
      return Circle_2(Hyperbolic_point_2(x, y), sq_radius, CLOCKWISE);

    return Circle_2(Hyperbolic_point_2(x, y), sq_radius, COUNTERCLOCKWISE);
  }

private:
  const Traits& _gt;
}; // end Construct_supporting_circle_of_bisector

// only used internally
template <typename Traits>
class Construct_intersection_2
{
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Hyperbolic_point_2                 Hyperbolic_point_2;
  typedef typename Traits::Hyperbolic_segment_2               Hyperbolic_segment_2;
  typedef typename Traits::Euclidean_segment_2                Euclidean_segment_2;
  typedef typename Traits::Euclidean_line_2                   Euclidean_line_2;
  typedef typename Traits::Circle_2                           Circle_2;
  typedef typename Traits::Circular_arc_2                     Circular_arc_2;

public:
  Construct_intersection_2(const Traits& gt = Traits()) : _gt(gt) { }

  Hyperbolic_point_2 operator()(const Euclidean_line_2& /*ell1*/,
                                const Euclidean_line_2& /*ell2*/)
  {
    // The only point where two Euclidean lines can intersect in the Poincaré disk is the origin.
    return Hyperbolic_point_2(0,0);
  }

  std::pair<Hyperbolic_point_2, Hyperbolic_point_2> operator()(const Euclidean_line_2& ell,
                                                               const Circle_2& c)
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

  std::pair<Hyperbolic_point_2, Hyperbolic_point_2> operator()(const Circle_2& c,
                                                               const Euclidean_line_2& ell)
  {
    return operator()(ell, c);
  }

  std::pair<Hyperbolic_point_2, Hyperbolic_point_2> operator()(const Circle_2& c1, const Circle_2& c2)
  {
    FT xa = c1.center().x(), ya = c1.center().y();
    FT xb = c2.center().x(), yb = c2.center().y();
    FT d2 = _gt.compute_squared_distance_2_object()(c1.center(), c2.center());
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
        pair<Hyperbolic_point_2, Hyperbolic_point_2> res = operator()(c1->circle(),
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
        pair<Hyperbolic_point_2, Hyperbolic_point_2> res = operator()(ell1->supporting_line(),
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
        Hyperbolic_point_2 p1 = operator()(ell1->supporting_line(), ell2->supporting_line());
        CGAL_assertion(p1.x()*p1.x() + p1.y()*p1.y() < FT(1));
        return p1;
      }
    }
  }
private:
  const Traits& _gt;
};

template <typename Traits>
class Construct_hyperbolic_segment_2
{
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Point_2                            Bare_point;
  typedef typename Traits::Weighted_point_2                   Weighted_point_2;
  typedef typename Traits::Euclidean_segment_2                Euclidean_segment_2;
  typedef typename Traits::Hyperbolic_point_2                 Hyperbolic_point_2;
  typedef typename Traits::Hyperbolic_segment_2               Hyperbolic_segment_2;
  typedef typename Traits::Circle_2                           Circle_2;
  typedef typename Traits::Circular_arc_2                     Circular_arc_2;

  typedef typename Traits::Construct_weighted_circumcenter_2  Construct_weighted_circumcenter_2;

public:
  Construct_hyperbolic_segment_2(const Traits& gt = Traits()) : _gt(gt) {}

  typedef Hyperbolic_segment_2                                result_type;

  Hyperbolic_segment_2 operator()(const Hyperbolic_point_2& p,
                                  const Hyperbolic_point_2& q) const
  {
    Origin o;
    if(_gt.collinear_2_object()(p, q, Hyperbolic_point_2(o)))
      return Euclidean_segment_2(p, q);

    Weighted_point_2 wp(p);
    Weighted_point_2 wq(q);
    Weighted_point_2 wo(Hyperbolic_point_2(o), FT(1)); // Poincaré circle

    Bare_point center = _gt.construct_weighted_circumcenter_2_object()(wp, wo, wq);
    FT sq_radius = _gt.compute_squared_distance_2_object()(p, center);

    Circle_2 circle = _gt.construct_circle_2_object()(center, sq_radius);
    // uncomment!!!
    //assert(circle.has_on_boundary(p) && circle.has_on_boundary(q));

    if(_gt.orientation_2_object()(p, q, center) == LEFT_TURN)
      return Circular_arc_2(circle, p, q);

    return Circular_arc_2(circle, q, p);
  }

private:
  const Traits& _gt;
}; // end Construct_hyperbolic_segment_2

template <typename Traits>
class Construct_hyperbolic_circumcenter_2
{
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Circle_2                           Circle_2;
  typedef typename Traits::Hyperbolic_point_2                 Hyperbolic_point_2;
  typedef typename Traits::Hyperbolic_Voronoi_point_2         Hyperbolic_Voronoi_point_2;
  typedef typename Traits::Euclidean_line_2                   Euclidean_line_2;
  typedef typename Traits::Euclidean_circle_or_line_2         Euclidean_circle_or_line_2;

public:
  typedef Hyperbolic_Voronoi_point_2                          result_type;

  Construct_hyperbolic_circumcenter_2(const Traits& gt = Traits()) : _gt(gt) {}

  Hyperbolic_Voronoi_point_2 operator()(const Hyperbolic_point_2& p,
                                        const Hyperbolic_point_2& q,
                                        const Hyperbolic_point_2& r) const
  {
    Hyperbolic_point_2 po(CGAL::ORIGIN);
    Circle_2 l_inf(po, FT(1));

    Construct_circle_or_line_supporting_bisector<Traits> cclsb(_gt);
    Construct_intersection_2<Traits> ci(_gt);

    Euclidean_circle_or_line_2 bis_pq = cclsb(p, q);
    Euclidean_circle_or_line_2 bis_qr = cclsb(q, r);

    if(_gt.compare_distance_2_object()(po, p, q) == EQUAL &&
       _gt.compare_distance_2_object()(po, p, r) == EQUAL)
      return po;

    // now supporting objects cannot both be Euclidean lines
    Euclidean_line_2* l;
    Circle_2* c;

    if(Circle_2* c_pq = boost::get<Circle_2>(&bis_pq))
    {
      if(Circle_2* c_qr = boost::get<Circle_2>(&bis_qr))
      {
        std::pair<Hyperbolic_point_2, Hyperbolic_point_2> inters = ci(*c_pq, *c_qr);

        if(_gt.has_on_bounded_side_2_object()(l_inf, inters.first))
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

    std::pair<Hyperbolic_point_2, Hyperbolic_point_2> inters = ci(*c, *l);

    if(_gt.has_on_bounded_side_2_object()(l_inf, inters.first))
      return inters.first;

    return inters.second;
  }

private:
  const Traits& _gt;
}; // end Construct_hyperbolic_circumcenter_2

template <typename Traits>
class Construct_hyperbolic_bisector_2
{
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Hyperbolic_point_2                 Hyperbolic_point_2;
  typedef typename Traits::Hyperbolic_segment_2               Hyperbolic_segment_2;
  typedef typename Traits::Euclidean_segment_2                Euclidean_segment_2;
  typedef typename Traits::Euclidean_line_2                   Euclidean_line_2;
  typedef typename Traits::Euclidean_circle_or_line_2         Euclidean_circle_or_line_2;
  typedef typename Traits::Circle_2                           Circle_2;
  typedef typename Traits::Circular_arc_2                     Circular_arc_2;

public:
  Construct_hyperbolic_bisector_2(const Traits& gt = Traits())
    : _gt(gt)
  {}

  // constructs a hyperbolic line
  Hyperbolic_segment_2 operator()(const Hyperbolic_point_2& p,
                                  const Hyperbolic_point_2& q) const
  {
    Construct_circle_or_line_supporting_bisector<Traits> cclsb(_gt);
    Construct_intersection_2<Traits> ci(_gt);

    Hyperbolic_point_2 po(CGAL::ORIGIN);
    Circle_2 l_inf(po, FT(1));

    if(_gt.compare_distance_2_object()(po, p, q) == EQUAL)
    {
      Euclidean_line_2 l = _gt.construct_bisector_2_object()(p, q);
      std::pair<Hyperbolic_point_2, Hyperbolic_point_2> inters = ci(l, l_inf);
      return Euclidean_segment_2(inters.first, inters.second);
    }

    Euclidean_circle_or_line_2 bis_pq = cclsb(p, q);
    Circle_2* c = boost::get<Circle_2>(&bis_pq);
    std::pair<Hyperbolic_point_2, Hyperbolic_point_2> inters = ci(*c, l_inf);

    if(_gt.orientation_2_object()(c->center(), inters.first, inters.second) == POSITIVE)
      return Circular_arc_2(*c, inters.first, inters.second);
    else
      return Circular_arc_2(*c, inters.second, inters.first);
  }

  // constructs the hyperbolic bisector of segment [p, q] limited by
  // circumcenter(p, q, r) on one side
  // and circumcenter(p, s, q) on the other side
  Hyperbolic_segment_2 operator()(const Hyperbolic_point_2& p,
                                  const Hyperbolic_point_2& q,
                                  const Hyperbolic_point_2& r,
                                  const Hyperbolic_point_2& s) const
  {
    CGAL_triangulation_precondition((_gt.orientation_2_object()(p, q, r) == ON_POSITIVE_SIDE) &&
                                    (_gt.orientation_2_object()(p, s, q) == ON_POSITIVE_SIDE));
    CGAL_triangulation_precondition((_gt.side_of_oriented_circle_2_object()(p, q, r,s) == ON_NEGATIVE_SIDE) &&
                                    (_gt.side_of_oriented_circle_2_object()(p, s, q, r) == ON_NEGATIVE_SIDE));

    Construct_hyperbolic_circumcenter_2<Traits> chc(_gt);
    Construct_circle_or_line_supporting_bisector<Traits> cclsb(_gt);

    Hyperbolic_point_2 po(CGAL::ORIGIN);

    // TODO MT this is non-optimal...
    // the bisector is already computed here
    // and it will be recomputed below
    Hyperbolic_point_2 a = chc(p, q, r);
    Hyperbolic_point_2 b = chc(p, s, q);

    if(_gt.compare_distance_2_object()(po, p, q) == EQUAL)
    {
      Euclidean_line_2 l = _gt.construct_bisector_2_object()(p, q);
      return Euclidean_segment_2(a, b);
    }

    Euclidean_circle_or_line_2 bis_pq = cclsb(p, q);
    Circle_2* c_pq = boost::get<Circle_2>(&bis_pq);

    if(_gt.compare_distance_2_object()(po, p, q) == POSITIVE) {
      // then p is inside the supporting circle
      return Circular_arc_2(*c_pq, b, a);
    }

    return Circular_arc_2(*c_pq, a, b);
  }

  // constructs the hyperbolic bisector of segment [p, q]
  // limited by hyperbolic circumcenter(p, q, r) on one side
  // and going to the infinite line on the other side
  Hyperbolic_segment_2 operator()(const Hyperbolic_point_2& p,
                                  const Hyperbolic_point_2& q,
                                  const Hyperbolic_point_2& r) const
  {
    CGAL_triangulation_precondition(_gt.orientation_2_object()(p, q, r) == POSITIVE);

    Construct_circle_or_line_supporting_bisector<Traits> cclsb(_gt);
    Construct_hyperbolic_circumcenter_2<Traits> chc(_gt);
    Construct_intersection_2<Traits> ci(_gt);

    Hyperbolic_point_2 po(CGAL::ORIGIN);
    Circle_2 l_inf(po, FT(1));

    // TODO MT this is non-optimal...
    // the bisector is computed (at least) twice
    Hyperbolic_point_2 a = chc(p, q, r);

    if(_gt.compare_distance_2_object()(po, p, q) == EQUAL)
    {
      Euclidean_line_2 bis_pq = _gt.construct_bisector_2_object()(p, q);

      std::pair<Hyperbolic_point_2,Hyperbolic_point_2> inters = ci(bis_pq, l_inf);
      if(_gt.less_y_2_object()(p, q))
        return Euclidean_segment_2(a,inters.first);

      return Euclidean_segment_2(a,inters.second);
    }

    Euclidean_circle_or_line_2 bis_pq = cclsb(p, q);
    Circle_2* c_pq = boost::get<Circle_2>(&bis_pq);

    std::pair<Hyperbolic_point_2,Hyperbolic_point_2> inters = ci(*c_pq, l_inf);
    Hyperbolic_point_2 approx_pinf = inters.first;
    if(_gt.orientation_2_object()(p, q,inters.first) == NEGATIVE)
    {
      if(_gt.orientation_2_object()(c_pq->center(),a,inters.first) == POSITIVE)
        return Circular_arc_2(*c_pq, a, inters.first);

      return Circular_arc_2(*c_pq, inters.first, a);
    }

    if(_gt.orientation_2_object()(c_pq->center(),a,inters.first) == POSITIVE)
      return Circular_arc_2(*c_pq, inters.second, a);

    return Circular_arc_2(*c_pq, a, inters.second);
  }

private:
  const Traits& _gt;
}; // end Construct_hyperbolic_bisector_2

// For details see the JoCG paper (5:56-85, 2014)
template <typename Traits>
class Is_Delaunay_hyperbolic
{
public:
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Hyperbolic_point_2                 Hyperbolic_point_2;
  typedef typename Traits::Direction_2                        Direction_2;
  typedef typename Traits::Point_3                            Point_3;
  typedef typename Traits::Vector_3                           Vector_3;

  Is_Delaunay_hyperbolic(const Traits& gt = Traits())
    : _gt(gt)
  {}

  bool operator()(const Hyperbolic_point_2& p0,
                  const Hyperbolic_point_2& p1,
                  const Hyperbolic_point_2& p2) const
  {
    // @todo use _gt
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

private:
  const Traits& _gt;
}; // end Is_Delaunay_hyperbolic

template <typename Traits>
class Side_of_oriented_hyperbolic_segment_2
{
  typedef typename Traits::FT                                 FT;
  typedef typename Traits::Point_2                            Bare_point;
  typedef typename Traits::Weighted_point_2                   Weighted_point_2;
  typedef typename Traits::Hyperbolic_point_2                 Hyperbolic_point_2;
  typedef typename Traits::Circle_2                           Circle_2;
  typedef typename Traits::Euclidean_line_2                   Euclidean_line_2;

  typedef typename Traits::Construct_weighted_circumcenter_2  Construct_weighted_circumcenter_2;

public:
  typedef Oriented_side                                       result_type;

  Side_of_oriented_hyperbolic_segment_2(const Traits& gt = Traits()) : _gt(gt) {}

  result_type operator()(const Hyperbolic_point_2& p,
                         const Hyperbolic_point_2& q,
                         const Hyperbolic_point_2& query) const
  {
    // Check first if the points are collinear with the origin
    Circle_2 poincare(Hyperbolic_point_2(FT(0),FT(0)), FT(1));
    Hyperbolic_point_2 O(FT(0), FT(0));
    Orientation ori = _gt.orientation_2_object()(p, q, O);
    if(ori == COLLINEAR)
    {
      Euclidean_line_2 seg(p, q);
      Orientation qori = _gt.orientation_2_object()(p, q, query);
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

    Bare_point center = _gt.construct_weighted_circumcenter_2_object()(wp, wo, wq);
    FT sq_radius = _gt.compute_squared_distance_2_object()(p, center);

    Circle_2 circle = _gt.construct_circle_2_object()(center, sq_radius);
    Bounded_side bs = _gt.bounded_side_2_object()(circle, query);
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

private:
  const Traits& _gt;
};

} // end namespace internal

template<typename Kernel = CGAL::Cartesian<CORE::Expr> >
class Hyperbolic_Delaunay_triangulation_traits_2
  : public Kernel
{
  typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>    Self;
  typedef Kernel                                                Base;

public:
  typedef typename Kernel::FT                                   FT;
  typedef typename Kernel::Point_2                              Hyperbolic_point_2;
  typedef Hyperbolic_point_2                                    Hyperbolic_Voronoi_point_2;
  typedef typename Kernel::Circle_2                             Circle_2;
  typedef typename Kernel::Line_2                               Euclidean_line_2;
  typedef boost::variant<Circle_2,Euclidean_line_2>             Euclidean_circle_or_line_2;
  typedef internal::HDT_2_Circular_arc_2<Self>                  Circular_arc_2;
  typedef typename Kernel::Segment_2                            Euclidean_segment_2; // only used internally here
  typedef boost::variant<Circular_arc_2, Euclidean_segment_2>   Hyperbolic_segment_2;

  typedef typename Kernel::Triangle_2                           Hyperbolic_triangle_2;

  // only kept for demo to please T2graphicsitems
  typedef Euclidean_segment_2                                   Line_segment_2;
  typedef Hyperbolic_segment_2                                  Segment_2;

  // The objects Ray_2, Iso_rectangle_2 and Line_2 are needed by the CGAL::Qt::PainterOstream
  typedef typename Kernel::Direction_2                          Direction_2;
  typedef typename Kernel::Ray_2                                Ray_2;
  typedef Euclidean_line_2                                      Line_2;
  typedef typename Kernel::Iso_rectangle_2                      Iso_rectangle_2;

  // wrong names kept for demo
  typedef internal::Construct_hyperbolic_segment_2<Self>        Construct_hyperbolic_segment_2;
  typedef Construct_hyperbolic_segment_2                        Construct_segment_2;

  typedef internal::Construct_hyperbolic_circumcenter_2<Self>   Construct_hyperbolic_circumcenter_2;
  typedef internal::Construct_hyperbolic_bisector_2<Self>       Construct_hyperbolic_bisector_2;
  typedef internal::Is_Delaunay_hyperbolic<Self>                Is_Delaunay_hyperbolic;
  typedef internal::Side_of_oriented_hyperbolic_segment_2<Self> Side_of_oriented_hyperbolic_segment_2;

public:
  Hyperbolic_Delaunay_triangulation_traits_2(const Base& kernel = Base())
    : Base(kernel)
  {}

  Construct_hyperbolic_segment_2
  construct_hyperbolic_segment_2_object() const
  { return Construct_hyperbolic_segment_2(*this); }

  Construct_segment_2
  construct_segment_2_object() const
  { return Construct_hyperbolic_segment_2(*this); }

  Construct_hyperbolic_circumcenter_2
  construct_hyperbolic_circumcenter_2_object() const
  { return Construct_hyperbolic_circumcenter_2(*this); }

  Construct_hyperbolic_bisector_2
  construct_hyperbolic_bisector_2_object() const
  { return Construct_hyperbolic_bisector_2(*this); }

  Is_Delaunay_hyperbolic
  is_Delaunay_hyperbolic_2_object() const
  { return Is_Delaunay_hyperbolic(*this); }

  Side_of_oriented_hyperbolic_segment_2
  side_of_oriented_hyperbolic_segment_2_object() const
  { return Side_of_oriented_hyperbolic_segment_2(*this); }

}; // class Hyperbolic_Delaunay_triangulation_traits_2

} // namespace CGAL

#endif // CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
