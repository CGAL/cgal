// Copyright (c) 2010   INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://mbogdanov@scm.gforge.inria.fr/svn/cgal/trunk/Triangulation_2/include/CGAL/Triangulation_hyperbolic_traits_2.h $
// $Id: Triangulation_hyperbolic_traits_2.h 57323 2010-07-05 10:07:39Z sloriot $
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_TRIANGULATION_HYPERBOLIC_TRAITS_2_H
#define CGAL_TRIANGULATION_HYPERBOLIC_TRAITS_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>

#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include "boost/tuple/tuple.hpp"
#include "boost/variant.hpp"

namespace CGAL {



  template < class R >
  class Triangulation_hyperbolic_traits_2 {
  public:
    typedef Triangulation_hyperbolic_traits_2<R> Self;
    
    typedef R Triangulation_euclidean_traits_2;
    
    typedef R Rep;
    typedef typename R::RT          RT;
    typedef typename R::Point_2     Point_2;
    typedef typename R::Vector_2    Vector_2;
    typedef typename R::Triangle_2  Triangle_2;
    typedef typename R::Line_2      Line_2;
    typedef typename R::Ray_2       Ray_2;
    
    typedef typename R::Vector_3    Vector_3;
    typedef typename R::Point_3     Point_3;
    
    typedef typename R::Less_x_2                   Less_x_2;
    typedef typename R::Less_y_2                   Less_y_2;
    typedef typename R::Compare_x_2                Compare_x_2;
    typedef typename R::Compare_y_2                Compare_y_2;
    typedef typename R::Orientation_2              Orientation_2;
    typedef typename R::Side_of_oriented_circle_2  Side_of_oriented_circle_2;
    typedef typename R::Construct_bisector_2       Construct_bisector_2;
    typedef typename R::Compare_distance_2         Compare_distance_2;
    typedef typename R::Construct_triangle_2       Construct_triangle_2;
    typedef typename R::Construct_direction_2      Construct_direction_2;
    
    typedef typename R::Angle_2                       Angle_2;
    typedef typename R::Construct_midpoint_2          Construct_midpoint_2;
    typedef typename R::Compute_squared_distance_2    Compute_squared_distance_2;
    
    typedef typename R::Iso_rectangle_2 Iso_rectangle_2;
    typedef typename R::Circle_2 Circle_2;
    
    typedef boost::tuple<Circle_2, Point_2, Point_2> Arc_2;
    typedef typename R::Segment_2 Line_segment_2;
    typedef boost::variant<Arc_2, Line_segment_2> Segment_2;
    
    typedef typename R::Line_2 Euclidean_line_2;
    
  private:
    // Poincaré disk
    const Circle_2 _unit_circle;
  
  public:
    const Circle_2& unit_circle() const
    {
      return _unit_circle;
    }
    
    Angle_2
    angle_2_object() const
    { return Angle_2(); }
    
    Compute_squared_distance_2
    compute_squared_distance_2_object() const
    { return Compute_squared_distance_2(); }
    
    class Construct_segment_2
    {
      typedef typename CGAL::Regular_triangulation_filtered_traits_2<R> Regular_geometric_traits_2;
      typedef typename Regular_geometric_traits_2::Construct_weighted_circumcenter_2 Construct_weighted_circumcenter_2;
      typedef typename Regular_geometric_traits_2::Weighted_point_2 Weighted_point_2;
      typedef typename Regular_geometric_traits_2::Bare_point Bare_point;
      
    public:
      Construct_segment_2(const Circle_2& c) : _unit_circle(c)
      {
      }
      
      Segment_2 operator()(const Point_2& p, const Point_2& q) const
      {
        typedef typename R::Collinear_2 Collinear_2;
        if(Collinear_2()(p, q, _unit_circle.center())){
          return Line_segment_2(p, q);
        }
        
        Weighted_point_2 wp(p);
        Weighted_point_2 wq(q);
        Weighted_point_2 wo(_unit_circle.center(), _unit_circle.squared_radius());
        
        Bare_point center = Construct_weighted_circumcenter_2()(wp, wo, wq);
        FT radius = Compute_squared_distance_2()(p, center);
        
        Circle_2 circle( center, radius);
        // uncomment!!!
        //assert(circle.has_on_boundary(p) && circle.has_on_boundary(q));
                
        if(Orientation_2()(p, q, center) == LEFT_TURN) {
          return Arc_2(circle, p, q);
        }
        return Arc_2(circle, q, p);
      }
      
    private:
      const Circle_2& _unit_circle;
    };
    
    Construct_segment_2
    construct_segment_2_object() const
    {
      return Construct_segment_2(_unit_circle);
    }
    
    class Construct_circumcenter_2
    {
    public:
      Construct_circumcenter_2(const Circle_2& c) : _unit_circle(c)
      {}
      
      // TODO: improve this function
      Point_2	operator()(Point_2 p, Point_2 q, Point_2 r)
      {        
        assert(_unit_circle.bounded_side(p) == ON_BOUNDED_SIDE);
        assert(_unit_circle.bounded_side(q) == ON_BOUNDED_SIDE);
        assert(_unit_circle.bounded_side(r) == ON_BOUNDED_SIDE);
        
        Circle_2 circle(p, q, r);
        // circle must be inside the unit one
        assert(CGAL::do_intersect(_unit_circle, circle) == false);
               
        if(circle.center() <= _unit_circle.center() && circle.center() >= _unit_circle.center()){
          return _unit_circle.center();
        }
        
        FT x0 = circle.center().x(), y0 = circle.center().y();
        // a*alphaˆ2 + b*alpha + c = 0;
        FT a = x0*x0 + y0*y0;
        FT b = a - circle.squared_radius() + _unit_circle.squared_radius();
        FT c = _unit_circle.squared_radius();
        FT D = b*b - 4*a*c;
        
        FT alpha = (b - CGAL::sqrt(to_double(D)))/(2*a);
        
        Point_2 center(x0*alpha, y0*alpha);
        if(!circle.has_on_bounded_side(center)) 
        { std::cout << "Center does not belong to the pencil of spheres!!!" << std::endl;} ;
        return center;
      }
      
    private:
      const Circle_2 _unit_circle;
    };
    
    Construct_circumcenter_2
    construct_circumcenter_2_object()
    {
      Construct_circumcenter_2(_unit_circle);
    }
    
    Construct_midpoint_2
    construct_midpoint_2_object() const
    { return Construct_midpoint_2(); }
    
    //for natural_neighbor_coordinates_2
    typedef typename R::FT                         FT;
    typedef typename R::Equal_x_2                  Equal_x_2;
    typedef typename R::Compute_area_2             Compute_area_2;
    Compute_area_2 compute_area_2_object () const 
    {
      return Compute_area_2();
    }
    
    // for compatibility with previous versions
    typedef Point_2      Point;
    typedef Segment_2    Segment;
    typedef Triangle_2   Triangle;
    typedef Ray_2        Ray;
    //typedef Line_2       Line;
    
    Triangulation_hyperbolic_traits_2() : 
      _unit_circle(Point_2(0, 0), 1*1)
    {}
    
    Triangulation_hyperbolic_traits_2(FT r) : 
      _unit_circle(Point_2(0, 0), r*r)
    {}    
    
    Triangulation_hyperbolic_traits_2(const Triangulation_hyperbolic_traits_2 & other) : 
      _unit_circle(other._unit_circle)
    {}
    
    Triangulation_hyperbolic_traits_2 &operator=
    (const Triangulation_hyperbolic_traits_2 &)
    {
      return *this;
    }
    
    Less_x_2
    less_x_2_object() const
    { return Less_x_2();}
    
    Less_y_2
    less_y_2_object() const
    { return Less_y_2();}
    
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
    {return Side_of_oriented_circle_2();}
    
    Construct_circumcenter_2
    construct_circumcenter_2_object() const
    {
      return Construct_circumcenter_2(_unit_circle);
    }
    
    class Construct_hyperbolic_bisector_2
    {    
    public:      
      Construct_hyperbolic_bisector_2(const Circle_2& unit_circle) :
        _unit_circle(unit_circle) {}
      
      Segment_2 operator()(Point_2 p, Point_2 q) const
      {
        // If two points are almost of the same distance to the origin, then
        // the bisector is supported by the circle of huge radius etc.
        // This circle is computed inexactly.
        // At present time, in this case the bisector is supported by the line.
        
        Compute_squared_distance_2 dist = Compute_squared_distance_2();
        Point origin = _unit_circle.center();
        FT dif = dist(origin, p) - dist(origin, q);
        FT eps = 0.0000000001;
        
        // Bisector is straight in euclidean sense
        if(dif > -eps && dif < eps){
          
        // ideally
        //if(Compare_distance_2()(_unit_circle.center(), p, q) == EQUAL){
        
          // TODO: calling R::Construct_bisector
          Euclidean_line_2 l = Construct_bisector_2()(p, q);
          // compute the ending points
          std::pair<Point_2, Point_2> points = find_intersection(l);
          // TODO: improve
          Vector_2 v(points.first, points.second);
          if(v*l.to_vector() > 0){
            return Line_segment_2(points.first, points.second);
          }
          return Line_segment_2(points.second, points.first);
        }
        
        Circle_2 c =  construct_supporting_circle(p, q);
        // compute the ending points
        std::pair<Point_2, Point_2> points = find_intersection(c);
        
        if(Orientation_2()(points.first, points.second, c.center()) == LEFT_TURN) {
          return Arc_2(c, points.first, points.second);
        }
        return Arc_2(c, points.second, points.first);
      }
      
    private:
      // The cirle belongs to the pencil with limit points p and q
      Circle_2 construct_supporting_circle(Point_2 p, Point_2 q) const
      {
        // p, q are zero-circles
        // (x, y, xˆ2 + yˆ2 - rˆ2) = alpha*(xp, yp, xpˆ2 + ypˆ2) + (1-alpha)*(xq, yq, xqˆ2 + yqˆ2)
        // xˆ2 + yˆ2 - rˆ2 = Rˆ2, where R - is a radius of the given unit circle
        FT op = p.x()*p.x() + p.y()*p.y();
        FT oq = q.x()*q.x() + q.y()*q.y();
        FT alpha = (_unit_circle.squared_radius() - oq) / (op - oq);
        
        FT x = alpha*p.x() + (1-alpha)*q.x();
        FT y = alpha*p.y() + (1-alpha)*q.y();
        FT radius = x*x + y*y - _unit_circle.squared_radius();
        
        //improve
        typename R::Line_2 l = typename R::Construct_bisector_2()(p, q);
        Point_2 middle = Construct_midpoint_2()(p, q);
        Point_2 temp = middle + l.to_vector();
        if(Orientation_2()(middle, temp, Point_2(x, y)) == ON_POSITIVE_SIDE){
           return Circle_2(Point_2(x, y), radius, CLOCKWISE);
        }
        
        return Circle_2(Point_2(x, y), radius, COUNTERCLOCKWISE);
      }
      
      // Find intersection of an input circle orthogonal to the Poincaré disk
      // and the circle representing this disk
      
      // TODO: sqrt(to_double()?)
      std::pair<Point_2, Point_2> find_intersection(Circle_2& circle) const
      {
        FT x = circle.center().x(), y = circle.center().y();
        
        // axˆ2 + 2bˆx + c = 0;
        FT a = x*x + y*y;
        FT b = -_unit_circle.squared_radius() * x;
        FT c = _unit_circle.squared_radius()*_unit_circle.squared_radius() - _unit_circle.squared_radius()*y*y;
        assert(b*b - a*c > 0);
        FT D = CGAL::sqrt(to_double(b*b - a*c));
        
        FT x1 = (-b - D)/a;
        FT x2 = (-b + D)/a;
        FT y1 = (_unit_circle.squared_radius() - x1*x)/y;
        FT y2 = (_unit_circle.squared_radius() - x2*x)/y;
        
        return std::make_pair(Point_2(x1, y1), Point_2(x2, y2));
      }
      
      // Find intersection of an input line orthogonal to the Poincaré disk
      // and the circle representing this disk
      
      // TODO: sqrt(to_double()?)
      std::pair<Point_2, Point_2> find_intersection(Euclidean_line_2& l) const
      {
        typedef typename R::Vector_2 Vector_2;
        Vector_2 v = l.to_vector();
        
        // normalize the vector
        FT squared_coeff = _unit_circle.squared_radius()/v.squared_length();
        FT coeff = CGAL::sqrt(to_double(squared_coeff));
        
        Point_2 p1(coeff*v.x(), coeff*v.y());
        Point_2 p2(-p1.x(), -p1.y());
        return std::make_pair(p1, p2);
      }
      
    private:
      const Circle_2 _unit_circle;
    };
    
    Construct_hyperbolic_bisector_2
    construct_hyperbolic_bisector_2_object() const
    { return Construct_hyperbolic_bisector_2(_unit_circle);}
    
    Construct_bisector_2
    construct_bisector_2_object() const
    {return Construct_bisector_2();}
    
    Compare_distance_2
    compare_distance_2_object() const
    {return Compare_distance_2();}
    
    Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}
    
    Construct_direction_2  construct_direction_2_object() const
    {return Construct_direction_2();}
    
    class Construct_ray_2
    {
    public:
      Construct_ray_2(Circle_2 c) : 
        _unit_circle(c) {}
      
      Segment_2 operator()(Point_2 p, Segment_2 l) const
      { 
        if(typename R::Segment_2* s = boost::get<typename R::Segment_2>(&l)){
          return operator()(p, *s);
        }
        if(Arc_2* arc = boost::get<Arc_2>(&l)){
          if(arc->get<0>().orientation() == CLOCKWISE){
            arc->get<1>() = p;
            return *arc;
          }
          arc->get<2>() = p;
          return *arc;
        }
        assert(false);
        return Segment_2();
      }
      
      Segment_2 operator()(Point_2 p, typename R::Segment_2 s) const
      {
        return typename R::Segment_2(p, s.target());
      }
      
    private:
      
      const Circle_2 _unit_circle;
     };
    
    Construct_ray_2  construct_ray_2_object() const
    {return Construct_ray_2(_unit_circle);}
  };

  
// Take out the code below to some separate file
  
#ifdef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
template  <>
struct Triangulation_structural_filtering_traits< Triangulation_hyperbolic_traits_2<Epeck> > {
   typedef Tag_true Use_structural_filtering_tag;
};
#endif // CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H

#ifdef CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
template <>
struct Triangulation_structural_filtering_traits< Triangulation_hyperbolic_traits_2<Epick> > {
  typedef Tag_true Use_structural_filtering_tag;
};
#endif // CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
  
} //namespace CGAL 

#endif // CGAL_TRIANGULATION_HYPERBOLIC_TRAITS_2_H
