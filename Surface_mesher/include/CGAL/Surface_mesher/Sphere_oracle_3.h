// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
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
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_SURFACE_MESHER_SPHERE_ORACLE_3_H
#define CGAL_SURFACE_MESHER_SPHERE_ORACLE_3_H

#include <CGAL/Surface_mesher/Null_oracle_visitor.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/Origin.h>

#include <boost/tuple/tuple.hpp>

namespace CGAL {

  namespace Surface_mesher {

  template <
    class GT,
    class Point_creator = Creator_uniform_3<typename GT::FT,
                                            typename GT::Point_3>,
    class Visitor = Null_oracle_visitor
    >
  class Sphere_oracle_3
  {
    // private types
    typedef Sphere_oracle_3<GT, Point_creator, Visitor> Self;
    
    typedef typename GT::Point_3 Point;
    typedef typename GT::FT FT;
    typedef typename GT::Sphere_3 Sphere_3;

  public:

    // Public types
    typedef GT Geom_traits;
    typedef typename GT::Point_3 Point_3;
    typedef typename GT::Segment_3 Segment_3;
    typedef typename GT::Ray_3 Ray_3;
    typedef typename GT::Line_3 Line_3;

    typedef Sphere_3 Surface_3;

    typedef Point Intersection_point;
  private:
    // Private members
    Visitor visitor; // a visitor that can modify a point, before returning it.

  public:

    // Constructors
    Sphere_oracle_3 (Visitor visitor_ = Visitor() ) :
      visitor(visitor_)
    {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CONSTRUCTORS
#  ifndef CGAL_SURFACE_MESHER_IMPLICIT_SURFACE_ORACLE_3_H
      std::cerr << "CONS: Sphere_oracle_3\n";
#  endif
#endif
    }

    const Visitor& get_visitor() const
    {
      return visitor;
    }

    // Predicates and Constructions

    bool is_in_volume(const Surface_3& sphere, const Point& p) const
    {
      typename GT::Has_on_bounded_side_3 on_bounded_side_of_sphere =
        GT().has_on_bounded_side_3_object();

      return on_bounded_side_of_sphere(sphere, p);
    }

    class Intersect_3 
    {
      const Self& oracle;

      boost::tuple<int, FT, FT> 
      intersection_line_sphere_lambda(const Surface_3& sphere,
                                      const Point& a,
                                      const Point& b) const
      {
        /*
          Let the vectorial line equation:
            m = a + lambda * ( b - a )
          (a, b, and m are points, and lambda if a real.)
          
          Let c be the center of the sphere, of radius r.
          The intersection of the line and the sphere is given by:
            (c-m)^2 = r^2
          That is:
                       ((c-a)^2 - r^2)
            - 2 lambda (c-a)*(b-a)
            + lambda^2 (b-a)^2 == 0

            (second degre equation)

          deltaprime = delta/4 = ((c-a)(b-a))^2 - (b-a)^2 * ( (c-a)^2 -r^2 )

          if delta > 0, root_1 = ((c-a)(b-a) - \sqrt(delta/4)) / (b-a)^2
                        root_2 = ((c-a)(b-a) + \sqrt(delta/4)) / (b-a)^2
                 (root_1 < root_2)
        */

        typedef typename GT::Vector_3 Vector_3;

        typename GT::Construct_vector_3 vector =
          GT().construct_vector_3_object();
        typename GT::Compute_scalar_product_3 scalar_product = 
          GT().compute_scalar_product_3_object();
        typename GT::Compute_squared_distance_3 squared_distance = 
          GT().compute_squared_distance_3_object();
        typename GT::Construct_center_3 center =
          GT().construct_center_3_object();
        typename GT::Compute_squared_radius_3 squared_radius =
          GT().compute_squared_radius_3_object();

        const Point c = center(sphere);
        const Vector_3 ab = vector(a, b);
        const Vector_3 ac = vector(a, c);
        const FT ab_ac = scalar_product(ab, ac);
        const FT ab2 = squared_distance(a, b);
        const FT ac2 = squared_distance(a, c);
        const FT r2 = squared_radius(sphere);
        const FT deltaprime = ab_ac * ab_ac - ab2 * ( ac2 - r2 );

        switch( CGAL::sign(deltaprime) )
        {
        case ZERO:
          return boost::make_tuple(1, ab_ac / ab2, 0);
        case POSITIVE:
	  {
	    const FT sqrt_deltaprime = CGAL::sqrt(deltaprime);
	    return boost::make_tuple(2,
				     (ab_ac - sqrt_deltaprime) / ab2,
				     (ab_ac + sqrt_deltaprime) / ab2);
	  }
        case NEGATIVE:
          break;
        }
        return boost::make_tuple(0, 0, 0);
      } //end intersection_line_sphere_lambda

      template <class Assert_on_lambda>
      Object private_intersection(const Surface_3& sphere,
                                  const Point& a,
                                  const Point& b,
                                  Assert_on_lambda test) const
      {
        typedef typename GT::Vector_3 Vector;
        
        typename GT::Construct_vector_3 vector =
          GT().construct_vector_3_object();
        typename GT::Construct_scaled_vector_3 scaled_vector = 
          GT().construct_scaled_vector_3_object();
        typename GT::Construct_translated_point_3 translated_point = 
          GT().construct_translated_point_3_object();

        int number_of_roots;
        FT root_1, root_2;
        boost::tie(number_of_roots, root_1, root_2) = 
          intersection_line_sphere_lambda(sphere, a, b);

        const Vector ab = vector(a, b);
        if(number_of_roots > 0 && test(root_1))
        {
          Point p = translated_point(a, scaled_vector(ab, root_1));
          oracle.get_visitor().new_point(p);
          return make_object(p);
        }
        else if (number_of_roots > 1 && test(root_2))
        {
          Point p = translated_point(a, scaled_vector(ab, root_2));
          oracle.get_visitor().new_point(p);
          return make_object(p);
        }
        // else
        return Object();
      } // end private_intersection

      struct Lambda_between_0_and_1 : public std::unary_function<FT, bool> 
      {
        bool operator()(const FT x) const
        {
          return FT(0) <= x && x <= FT(1);
        }
      };
        
      struct Lambda_positive : public std::unary_function<FT, bool> 
      {
        bool operator()(const FT x) const
        {
          return FT(0) <= x;
        }
      };
        
      struct Always_true : public std::unary_function<FT, bool> 
      {
        bool operator()(const FT) const
        {
          return true;
        }
      };
        
    public:
      Intersect_3(const Self& oracle) : oracle(oracle)
      {
      }

      Object operator()(const Surface_3& sphere, Segment_3 s) const
      {
        typename GT::Construct_point_on_3 point_on =
          GT().construct_point_on_3_object();

        const Point& a = point_on(s, 0);
        const Point& b = point_on(s, 1);
        
        return private_intersection(sphere, a, b, Lambda_between_0_and_1());
      } // end operator()(Surface_3, Segment_3)

      Object operator()(const Surface_3& sphere, const Ray_3& r) const {
        typename GT::Construct_point_on_3 point_on =
          GT().construct_point_on_3_object();

        const Point& a = point_on(r, 0);
        const Point& b = point_on(r, 1);
        
        return private_intersection(sphere, a, b, Lambda_positive());
      } // end operator()(Surface_3, Ray_3)

      Object operator()(const Surface_3& sphere, const Line_3& l) const {
        typename GT::Construct_point_on_3 point_on =
          GT().construct_point_on_3_object();

        const Point& a = point_on(l, 0);
        const Point& b = point_on(l, 1);
        
        return private_intersection(sphere, a, b, Always_true());
      } // end operator()(Surface_3, Line_3)

      /** Modifies s = [a, b] by clipping it to sphere.
          Return false iff s is outside sphere. */
      bool clip_segment(const Surface_3& sphere,
                        Point_3& a,
                        Point_3& b) const
      {
        typedef typename GT::Vector_3 Vector;
        
        typename GT::Has_on_bounded_side_3 on_bounded_side_of_sphere =
          GT().has_on_bounded_side_3_object();
        typename GT::Construct_vector_3 vector =
          GT().construct_vector_3_object();
        typename GT::Construct_scaled_vector_3 scaled_vector = 
          GT().construct_scaled_vector_3_object();
        typename GT::Construct_translated_point_3 translated_point = 
          GT().construct_translated_point_3_object();

        const bool a_in_sphere = on_bounded_side_of_sphere(sphere, a);
        const bool b_in_sphere = on_bounded_side_of_sphere(sphere, b);        

        if( a_in_sphere && b_in_sphere )
          return true;

        int number_of_roots;
        FT root_1, root_2;
        
        boost::tie(number_of_roots, root_1, root_2) = 
          intersection_line_sphere_lambda(sphere, a, b);

#ifdef CGAL_SURFACE_MESHER_DEBUG_IMPLICIT_ORACLE
        std::cerr << "Clip segment. Roots=("
                  << root_1 << ", " << root_2 << ")\n";
#endif
        if( number_of_roots < 2 )
          return false;

        if( root_1 > FT(1) ) // root_x \in ]1,\infinity[
          return false;      // no intersection

        if( root_1 >= FT(0) ) // root_1 \in [0,1[
        {                     // move point a
          const Point original_a = a;
          const Vector ab = vector(a, b);
          a = translated_point(original_a, scaled_vector(ab, root_1));
          if( root_2 <= FT(1) ) /// move b iif root_2 <=1
          {
            b = translated_point(original_a, scaled_vector(ab, root_2));
          }
          return true;
        }
        else // root_1 in ]-\infinity, 0[
        {    // do not move point a
          if( root_2 < FT(0) ) // root_x in ]-\infinity, 0[
            return false;      // no intersection
          else 
          {
            const Vector ab = vector(a, b);
            if( root_2 <= FT(1) )
              b = translated_point(a, scaled_vector(ab, root_2));
            return true;
          }
        }
      }

      /** The return value s is r clipped to sphere.
          Return false iff r does not intersect sphere. */
      bool clip_ray(const Surface_3& sphere,
                    const Ray_3& r,
                    Point_3& a,
                    Point_3& b) const
      {
        typedef typename GT::Vector_3 Vector;
        
        typename GT::Construct_point_on_3 point_on =
          GT().construct_point_on_3_object();
        typename GT::Construct_vector_3 vector =
          GT().construct_vector_3_object();
        typename GT::Construct_scaled_vector_3 scaled_vector = 
          GT().construct_scaled_vector_3_object();
        typename GT::Construct_translated_point_3 translated_point = 
          GT().construct_translated_point_3_object();

        a = point_on(r, 0);
        b = point_on(r, 1);

        int number_of_roots;
        FT root_1, root_2;
        
        boost::tie(number_of_roots, root_1, root_2) = 
          intersection_line_sphere_lambda(sphere, a, b);

        if( number_of_roots == 2 && root_2 > FT(0) )
        {
          const Vector ab = vector(a, b);
          b = translated_point(a, scaled_vector(ab, root_2));
          if(root_1 > FT(0))
            a = translated_point(a, scaled_vector(ab, root_1));
          // if root_1 <= 0, a is in the ball
          return true;
        }
        // else r does not intersect the sphere
        return false;
      } // end clip_ray

      /** The return value s=(ab) is l clipped to sphere.
          Return false iff l does not intersect sphere. */
      bool clip_line(const Surface_3& sphere, const Line_3& l,
                     Point& a,
                     Point& b) const
      {
        typedef typename GT::Vector_3 Vector;
        
        typename GT::Construct_point_on_3 point_on =
          GT().construct_point_on_3_object();
        typename GT::Construct_vector_3 vector =
          GT().construct_vector_3_object();
        typename GT::Construct_scaled_vector_3 scaled_vector = 
          GT().construct_scaled_vector_3_object();
        typename GT::Construct_translated_point_3 translated_point = 
          GT().construct_translated_point_3_object();

        a = point_on(l, 0);
        b = point_on(l, 1);

        int number_of_roots;
        FT root_1, root_2;
        
        boost::tie(number_of_roots, root_1, root_2) = 
          intersection_line_sphere_lambda(sphere, a, b);

        if( number_of_roots == 2 )
        {
          const Point original_a = a;
          const Vector ab = vector(a, b);
          a = translated_point(original_a, scaled_vector(ab, root_1));
          b = translated_point(original_a, scaled_vector(ab, root_2));
          return true;
        }
        // else l does not intersect the sphere
        return false;
      } // end clip_line

    }; // end nested class Intersect_3

    class Construct_initial_points
    {
      const Self& oracle;
    public:
      Construct_initial_points(const Self& oracle) : oracle(oracle)
      {
      }
      
      // Random points
      template <typename OutputIteratorPoints>
      OutputIteratorPoints operator() (const Surface_3& sphere, 
                                       OutputIteratorPoints out, 
                                       int n = 20) const // WARNING: why 20?
      {
        const Point center = 
          GT().construct_center_3_object()(sphere);
        const FT squared_radius = 
          GT().compute_squared_radius_3_object()(sphere);
        const double radius_in_double = 
          CGAL::sqrt(CGAL::to_double(squared_radius));

        typename CGAL::Random_points_on_sphere_3<Point,
          Point_creator> random_point_on_sphere(radius_in_double);
        typename GT::Construct_vector_3 vector_3 =
          GT().construct_vector_3_object();
        typename GT::Construct_translated_point_3 translate =
          GT().construct_translated_point_3_object();

        while (n-->0)
        {
          Point p = translate(*random_point_on_sphere++,
                             vector_3(CGAL::ORIGIN, center));
          oracle.get_visitor().new_point(p);
          *out++ = p;
        }
        return out;
      }
    }; // end nested class Construct_initial_points

    Construct_initial_points construct_initial_points_object() const
    {
      return Construct_initial_points(*this);
    }

    Intersect_3 intersect_3_object() const
    {
      return Intersect_3(*this);
    }
  };  // end Sphere_oracle_3

  }  // namespace Surface_mesher

} // namespace CGAL


#endif  // CGAL_SURFACE_MESHER_SPHERE_ORACLE_3_H
