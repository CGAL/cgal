// Copyright (c) 2010-2018   INRIA Sophia Antipolis, INRIA Nancy (France).
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
// Author(s)     : Mikhail Bogdanov
//                 Monique Teillaud <Monique.Teillaud@inria.fr>

#ifndef CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_OLD_TRAITS_2_H
#define CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_OLD_TRAITS_2_H

//#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/determinant.h>
#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"

namespace CGAL {

template < class R >
class Hyperbolic_Delaunay_triangulation_old_traits_2
  : public R
  {
public:
  typedef typename R::FT          FT;

  typedef typename R::Point_2     Point_2;
  typedef Point_2                 Hyperbolic_point_2;
  typedef Hyperbolic_point_2      Hyperbolic_Voronoi_point_2;
  typedef typename R::Circle_2    Circle_2;
  typedef typename R::Triangle_2  Hyperbolic_triangle_2;

  typedef boost::tuple<Circle_2, Point_2, Point_2>    Arc_2;
  typedef typename R::Segment_2                       Euclidean_segment_2; //only used internally here
  typedef boost::variant<Arc_2, Euclidean_segment_2>  Hyperbolic_segment_2;

  typedef typename R::Compare_x_2                Compare_x_2;
  typedef typename R::Compare_y_2                Compare_y_2;
  typedef typename R::Orientation_2              Orientation_2;
  typedef typename R::Side_of_oriented_circle_2  Side_of_oriented_circle_2;

  // only kept for demo to please T2graphicsitems
  typedef Euclidean_segment_2  Line_segment_2;
  typedef Hyperbolic_segment_2 Segment_2;

  // the following types are only used internally in this traits class,
  // so they need not be documented, and they don't need _object()
  typedef typename R::Collinear_2                Euclidean_collinear_2;
  typedef typename R::Construct_bisector_2       Construct_Euclidean_bisector_2;
  typedef typename R::Construct_midpoint_2       Construct_Euclidean_midpoint_2;
  typedef typename R::Compute_squared_distance_2 Compute_squared_Euclidean_distance_2;
  typedef typename R::Line_2      Euclidean_line_2;
  typedef typename R::Vector_2    Vector_2;
  // used by Is_hyperbolic
  typedef typename R::Vector_3    Vector_3;
  typedef typename R::Point_3     Point_3;

  // MT useless?
  //  typedef Hyperbolic_Delaunay_triangulation_traits_2<R> Self;
  //  typedef typename R::RT          RT;
  //  typedef R Kernel;
  //  typedef R Rep;
  //  typedef typename R::Triangle_2  Triangle_2;
  //  typedef typename R::Line_2      Line_2;
  //  typedef typename R::Ray_2       Ray_2; // why would we need Eucldean rays??
  //  typedef typename R::Iso_rectangle_2 Iso_rectangle_2;
  //  typedef typename R::Angle_2                       Angle_2;

  //  typedef typename R::Less_x_2                   Less_x_2;
  //  typedef typename R::Less_y_2                   Less_y_2;
  //  typedef typename R::Compare_distance_2         Compare_distance_2;
  //  typedef typename R::Construct_triangle_2       Construct_triangle_2;
  //  typedef typename R::Construct_direction_2      Construct_direction_2;

public:

  class Side_of_oriented_hyperbolic_segment_2 {};


  class Construct_hyperbolic_segment_2
  {
    //typedef typename CGAL::Regular_triangulation_filtered_traits_2<R> Regular_geometric_traits_2;
    typedef R  Regular_geometric_traits_2;
    typedef typename Regular_geometric_traits_2::Construct_weighted_circumcenter_2 Construct_weighted_circumcenter_2;
    typedef typename Regular_geometric_traits_2::Weighted_point_2 Weighted_point_2;
    typedef typename Regular_geometric_traits_2::Bare_point Bare_point;

  public:
        Construct_hyperbolic_segment_2()
        {
        }

    Hyperbolic_segment_2 operator()(const Point_2& p, const Point_2& q) const
    {
      Origin o;
      if(Euclidean_collinear_2()(p, q, Point_2(o))){
        return Euclidean_segment_2(p, q);
      }

      Weighted_point_2 wp(p);
      Weighted_point_2 wq(q);
      Weighted_point_2 wo(Point_2(o), FT(1)); // Poincaré circle

      Bare_point center = Construct_weighted_circumcenter_2()(wp, wo, wq);
      FT radius = Compute_squared_Euclidean_distance_2()(p, center);

      Circle_2 circle( center, radius);
      // uncomment!!!
      //assert(circle.has_on_boundary(p) && circle.has_on_boundary(q));

      if(Orientation_2()(p, q, center) == LEFT_TURN) {
        return Arc_2(circle, p, q);
      }
      return Arc_2(circle, q, p);
    }

  };

  Construct_hyperbolic_segment_2
    construct_hyperbolic_segment_2_object() const
  { return Construct_hyperbolic_segment_2(); }

  // wrong names kept for demo
  typedef Construct_hyperbolic_segment_2 Construct_segment_2;
  Construct_segment_2
    construct_segment_2_object() const
  { return Construct_hyperbolic_segment_2(); }

  class Construct_circumcenter_2
  {
  public:

    // TODO: improve this function
    Point_2	operator()(Point_2 p, Point_2 q, Point_2 r)
    {
      CGAL_triangulation_assertion_code(Origin oo; Point_2 poo(oo); Circle_2 co(poo,FT(1)));
      CGAL_triangulation_assertion(co.bounded_side(p) == ON_BOUNDED_SIDE);
      CGAL_triangulation_assertion(co.bounded_side(q) == ON_BOUNDED_SIDE);
      CGAL_triangulation_assertion(co.bounded_side(r) == ON_BOUNDED_SIDE);

      Circle_2 circle(p, q, r);
      // circle must be inside the unit one
      CGAL_triangulation_assertion(do_intersect(co, circle) == false);

      Origin o;
      Point_2 po = Point_2(o);
      if(circle.center() == po)
  { return po; }

      FT x0 = circle.center().x(), y0 = circle.center().y();
      // a*alphaˆ2 + b*alpha + c = 0;
      FT a = x0*x0 + y0*y0;
      FT b = a - circle.squared_radius() + FT(1); // Poincare disc has radius 1
      FT D = b*b - 4*a;

      FT alpha = (b - CGAL::sqrt(to_double(D)))/(2*a);

      Point_2 center(x0*alpha, y0*alpha);
      if(!circle.has_on_bounded_side(center))
      { std::cout << "Center does not belong to the pencil of spheres!!!" << std::endl;} ;
      return center;
    }

  };

  Construct_circumcenter_2
    construct_circumcenter_2_object()
  { return Construct_circumcenter_2(); }

  Hyperbolic_Delaunay_triangulation_old_traits_2()
  {}

  Hyperbolic_Delaunay_triangulation_old_traits_2(const Hyperbolic_Delaunay_triangulation_old_traits_2 & other)
  {}

  Hyperbolic_Delaunay_triangulation_old_traits_2 &operator=
  (const Hyperbolic_Delaunay_triangulation_old_traits_2 &)
  {
    return *this;
  }

  Compare_x_2
    compare_x_2_object() const
  { return Compare_x_2();}

  Compare_y_2
    compare_y_2_object() const
  { return Compare_y_2();}

  Orientation_2
    orientation_2_object() const
  { return Orientation_2();}

  Side_of_oriented_circle_2
    side_of_oriented_circle_2_object() const
  { return Side_of_oriented_circle_2(); }

  Construct_circumcenter_2
    construct_circumcenter_2_object() const
  { return Construct_circumcenter_2(); }

  class Construct_hyperbolic_bisector_2
  {
  public:
  Construct_hyperbolic_bisector_2()
    {}

    Hyperbolic_segment_2 operator()(Point_2 p, Point_2 q) const
    {
      // If two points are almost of the same distance to the origin, then
      // the bisector is supported by the circle of huge radius etc.
      // This circle is computed inexactly.
      // At present time, in this case the bisector is supported by the line.

      Compute_squared_Euclidean_distance_2 dist = Compute_squared_Euclidean_distance_2();
      Origin o;
      Point_2 po = Point_2(o);
      FT dif = dist(po, p) - dist(po, q);
      FT eps = 0.0000000001;

      // Bisector is straight in euclidean sense
      if(dif > -eps && dif < eps){

        // ideally
        //if(Compare_distance_2()(origin, p, q) == EQUAL){

        // TODO: calling R::Construct_bisector
        Euclidean_line_2 l = Construct_Euclidean_bisector_2()(p, q);
        // compute the ending points
        std::pair<Point_2, Point_2> points = find_intersection(l);
        // TODO: improve
        Vector_2 v(points.first, points.second);
        if(v*l.to_vector() > 0){
          return Euclidean_segment_2(points.first, points.second);
        }
        return Euclidean_segment_2(points.second, points.first);
      }

      Circle_2 c =  construct_supporting_circle_of_bisector(p, q);
      // compute the ending points
      std::pair<Point_2, Point_2> points = find_intersection(c);

      if(Orientation_2()(points.first, points.second, c.center()) == LEFT_TURN) {
        return Arc_2(c, points.first, points.second);
      }
      return Arc_2(c, points.second, points.first);
    }

  private:
    // The cirle belongs to the pencil with limit points p and q
    Circle_2 construct_supporting_circle_of_bisector(Point_2 p, Point_2 q) const
    {
      // p, q are zero-circles
      // (x, y, xˆ2 + yˆ2 - rˆ2) = alpha*(xp, yp, xpˆ2 + ypˆ2) + (1-alpha)*(xq, yq, xqˆ2 + yqˆ2)
      // xˆ2 + yˆ2 - rˆ2 = Rˆ2, where R - is a radius of the given unit circle
      FT op = p.x()*p.x() + p.y()*p.y();
      FT oq = q.x()*q.x() + q.y()*q.y();
      FT alpha = (FT(1) - oq) / (op - oq); // Poincare disc has radius 1

      FT x = alpha*p.x() + (1-alpha)*q.x();
      FT y = alpha*p.y() + (1-alpha)*q.y();
      FT radius = x*x + y*y - FT(1);

      //improve
      Euclidean_line_2 l = Construct_Euclidean_bisector_2()(p, q);
      Point_2 middle = Construct_Euclidean_midpoint_2()(p, q);
      Point_2 temp = middle + l.to_vector();
      if(Orientation_2()(middle, temp, Point_2(x, y)) == ON_POSITIVE_SIDE){
        return Circle_2(Point_2(x, y), radius, CLOCKWISE);
      }

      return Circle_2(Point_2(x, y), radius, COUNTERCLOCKWISE);
    }

    // Find intersection of an input circle orthogonal to the Poincare disk
    // and the circle representing this disk

    // TODO: sqrt(to_double()?)
    std::pair<Point_2, Point_2> find_intersection(Circle_2& circle) const
    {
      FT x = circle.center().x(), y = circle.center().y();

      // axˆ2 + 2bˆx + c = 0;
      FT a = x*x + y*y;
      /* FT b = -_unit_circle.squared_radius() * x; */
      /* FT c = _unit_circle.squared_radius()*_unit_circle.squared_radius() - _unit_circle.squared_radius()*y*y; */
      FT b = -x;
      FT c = 1-y*y;
      assert(b*b - a*c > 0);
      FT D = CGAL::sqrt(to_double(b*b - a*c));

      FT x1 = (-b - D)/a;
      FT x2 = (-b + D)/a;
      FT y1 = (FT(1) - x1*x)/y;
      FT y2 = (FT(1) - x2*x)/y;

      return std::make_pair(Point_2(x1, y1), Point_2(x2, y2));
    }

    // Find intersection of an input line orthogonal to the Poincare disk
    // and the circle representing this disk

    // TODO: sqrt(to_double()?)
    std::pair<Point_2, Point_2> find_intersection(Euclidean_line_2& l) const
    {
      typedef typename R::Vector_2 Vector_2;
      Vector_2 v = l.to_vector();

      // normalize the vector
      FT squared_coeff = FT(1)/v.squared_length();
      FT coeff = CGAL::sqrt(to_double(squared_coeff));

      Point_2 p1(coeff*v.x(), coeff*v.y());
      Point_2 p2(-p1.x(), -p1.y());
      return std::make_pair(p1, p2);
    }

  };

  Construct_hyperbolic_bisector_2
  construct_hyperbolic_bisector_2_object() const
  { return Construct_hyperbolic_bisector_2(); }

  Construct_Euclidean_bisector_2
  construct_Euclidean_bisector_2_object() const
  { return Construct_Euclidean_bisector_2(); }

  class Construct_ray_2
  {
  public:
    Construct_ray_2()
      {}

    Hyperbolic_segment_2 operator()(Point_2 p, Hyperbolic_segment_2 l) const
    {
      if(Euclidean_segment_2* s = boost::get<Euclidean_segment_2>(&l)){
        return operator()(p, *s);
      }
      if(Arc_2* arc = boost::get<Arc_2>(&l)){
        if(get<0>(*arc).orientation() == CLOCKWISE){
          get<1>(*arc) = p;
          return *arc;
        }
        get<2>(*arc) = p;
        return *arc;
      }
      assert(false);
      return Hyperbolic_segment_2();
    }

    Hyperbolic_segment_2 operator()(Point_2 p, Euclidean_segment_2 s) const
    {
      return Euclidean_segment_2(p, s.target());
    }

  };

  Construct_ray_2
    construct_ray_2_object() const
  { return Construct_ray_2(); }

  // For details see the JoCG paper (5:56-85, 2014)
  class Is_Delaunay_hyperbolic
  {
  public:

    bool operator() (const Point_2& p0, const Point_2& p1, const Point_2& p2) const
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

      return dt0*dt0 + dt1*dt1 - dt2*dt2 < 0;
    }

    bool operator() (const Point_2& p0, const Point_2& p1, const Point_2& p2, int& ind) const
    {
      if(this->operator()(p0, p1, p2) == false) {
        ind = find_non_hyperbolic_edge(p0, p1, p2);
        return false;
      }
      return true;
    }

  private:

    // assume the face (p0, p1, p2) is non-hyperbolic
    int find_non_hyperbolic_edge(const Point_2& p0, const Point_2& p1, const Point_2& p2) const
    {
      typedef typename R::Direction_2 Direction_2;

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

      if(d.counterclockwise_in_between(d0, d1)) {
        return 2;
      }

      if(d.counterclockwise_in_between(d1, d2)) {
        return 0;
      }

      return 1;
    }
  };

  Is_Delaunay_hyperbolic
    Is_Delaunay_hyperbolic_object() const
  { return Is_Delaunay_hyperbolic(); }
};

// Take out the code below to some separate file

#ifdef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
template  <>
struct Triangulation_structural_filtering_traits< Hyperbolic_Delaunay_triangulation_old_traits_2<Epeck> > {
  typedef Tag_true Use_structural_filtering_tag;
};
#endif // CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H

#ifdef CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
template <>
struct Triangulation_structural_filtering_traits< Hyperbolic_Delaunay_triangulation_old_traits_2<Epick> > {
  typedef Tag_true Use_structural_filtering_tag;
};
#endif // CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H

} //namespace CGAL

#endif // CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_TRAITS_2_H
