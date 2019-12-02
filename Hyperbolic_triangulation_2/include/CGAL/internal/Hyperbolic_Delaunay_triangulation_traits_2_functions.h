// Copyright (c) 2010-2018  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mikhail Bogdanov
//                 Monique Teillaud <Monique.Teillaud@inria.fr>

#ifndef CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_FUNCTIONS
#define CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_FUNCTIONS

#include <CGAL/Origin.h>
#include <CGAL/enum.h>

namespace CGAL {

namespace internal {

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
    Hyperbolic_point_2 po = CGAL::ORIGIN;

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

} // end namespace CGAL


#endif // CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_FUNCTIONS
