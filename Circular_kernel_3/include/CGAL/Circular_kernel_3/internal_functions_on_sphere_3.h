// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_SPHERE_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_SPHERE_3_H

#include <CGAL/Circular_kernel_3/Intersection_traits.h>

namespace CGAL {
  namespace SphericalFunctors {

    template < class SK >
    typename SK::Algebraic_kernel::Polynomials_for_line_3
    get_equation( const typename SK::Line_3 & l)
    {
      typedef typename SK::Algebraic_kernel Algebraic_kernel;
      return Algebraic_kernel().construct_polynomials_for_line_3_object()
	  (l.to_vector().x(), l.point().x(), 
           l.to_vector().y(), l.point().y(),
           l.to_vector().z(), l.point().z());
    }

    template < class SK >
    typename SK::Polynomial_1_3
    get_equation( const typename SK::Plane_3 & s )
    {
      typedef typename SK::Algebraic_kernel Algebraic_kernel;
      return Algebraic_kernel().construct_polynomial_1_3_object()
	  ( s.a(), s.b(), s.c(), s.d() );
    }

    template < class SK >
    typename SK::Polynomial_for_spheres_2_3
    get_equation( const typename SK::Sphere_3 & s )
    {
      typedef typename SK::Point_3    Point_3;
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      Point_3 center = s.center();
      return Algebraic_kernel().construct_polynomial_for_spheres_2_3_object()
	  ( center.x(), center.y(), center.z(), s.squared_radius() );
    }

    template < class SK >
    typename SK::Polynomials_for_circle_3
    get_equation( const typename SK::Circle_3 & c )
    {
      typedef typename SK::Algebraic_kernel Algebraic_kernel;
      return std::make_pair ( Algebraic_kernel().construct_polynomial_for_spheres_2_3_object()
                               (c.center().x(), c.center().y(), 
                                c.center().z(), c.squared_radius()),
                              Algebraic_kernel().construct_polynomial_1_3_object()
                               (c.supporting_plane().a(),
                                c.supporting_plane().b(),
                                c.supporting_plane().c(),
                                c.supporting_plane().d()));
    }

    template < class SK >
    typename SK::Sphere_3
    construct_sphere_3(const typename SK::Polynomial_for_spheres_2_3 &eq)
    {
      typedef typename SK::Sphere_3 Sphere_3;
      typedef typename SK::Point_3  Point_3;
      return Sphere_3(Point_3(eq.a(),eq.b(),eq.c()),eq.r_sq());
    }

   
    template < class SK >
    inline
    typename SK::Linear_kernel::Bounded_side_3::result_type
    bounded_side(const typename SK::Sphere_3 &s,
                 const typename SK::Circular_arc_point_3 &p) {
      typedef typename SK::Algebraic_kernel Algebraic_kernel;
      typedef typename SK::Polynomial_for_spheres_2_3 Equation;
      Equation equation = get_equation<SK>(s);
      Sign sign = Algebraic_kernel().sign_at_object()(equation,p.rep().coordinates());
      if(sign == NEGATIVE) return ON_BOUNDED_SIDE;
      else if(sign == POSITIVE) return ON_UNBOUNDED_SIDE;
      else return ON_BOUNDARY;
    }

    template < class SK >
    inline
    typename SK::Linear_kernel::Bounded_side
    bounded_side(const typename SK::Circle_3 &c,
                 const typename SK::Circular_arc_point_3 &p) {
      typedef typename SK::Algebraic_kernel Algebraic_kernel;
      typedef typename SK::Polynomial_for_spheres_2_3 Equation;
      CGAL_kernel_assertion(SK().has_on_3_object()(c.supporting_plane(),p));
      Equation equation = get_equation<SK>(c.diametral_sphere());
      Sign sign = Algebraic_kernel().sign_at_object()(equation,p.rep().coordinates());
      if(sign == NEGATIVE) return ON_BOUNDED_SIDE;
      else if(sign == POSITIVE) return ON_UNBOUNDED_SIDE;
      else return ON_BOUNDARY;
    }

    template < class SK >
    inline
    bool
    non_oriented_equal(const typename SK::Sphere_3 & s1,
                       const typename SK::Sphere_3 & s2)
    {
      // Should we compare anyway even if they are degenerated?
      CGAL_kernel_assertion(!(s1.is_degenerate() || s2.is_degenerate()));
      return s1.center() == s2.center() &&
      s1.squared_radius() == s2.squared_radius();
    }

    template < class SK >
    inline
    bool
    non_oriented_equal(const typename SK::Plane_3 & p1,
                       const typename SK::Plane_3 & p2)
    {
      // Should we compare anyway even if they are degenerated?
      CGAL_kernel_assertion(!(p1.is_degenerate() || p2.is_degenerate()));
      if(is_zero(p1.a())) {
        if(!is_zero(p2.a())) return false;
        if(is_zero(p1.b())) {
          if(!is_zero(p2.b())) return false;
          return p1.c() * p2.d() == p1.d() * p2.c();
        }
        return (p2.c() * p1.b() == p1.c() * p2.b()) &&
               (p2.d() * p1.b() == p1.d() * p2.b());
      }
      return (p2.b() * p1.a() == p1.b() * p2.a()) &&
             (p2.c() * p1.a() == p1.c() * p2.a()) &&
             (p2.d() * p1.a() == p1.d() * p2.a());
    }

    template < class SK >
    inline
    bool
    non_oriented_equal(const typename SK::Circle_3 & c1,
                       const typename SK::Circle_3 & c2) {
      // We see degeneracies on the other non_oriented_equal functions
      if(!non_oriented_equal<SK>(c1.diametral_sphere(), c2.diametral_sphere())) return false;
      if(!non_oriented_equal<SK>(c1.supporting_plane(), c2.supporting_plane())) return false;
      return true;
    }

    template< class SK>
    bool
    non_oriented_equal(const typename SK::Line_3 &l1,
                       const typename SK::Line_3 &l2)
    {
      typedef typename SK::Vector_3 Vector_3;
      if(!SK().has_on_3_object()(l1, l2.point())) return false;

      const Vector_3& v1 = l1.to_vector();
      const Vector_3& v2 = l2.to_vector();

      if(v1.x() * v2.y() != v1.y() * v2.x()) return false;
      if(v1.x() * v2.z() != v1.z() * v2.x()) return false;
      if(v1.y() * v2.z() != v1.z() * v2.y()) return false;

      return true;
    }

    template< class SK>
    bool
    non_oriented_equal( const typename SK::Line_arc_3 &l1,
                        const typename SK::Line_arc_3 &l2)
    {
      if(!non_oriented_equal<SK>(l1.supporting_line(),
                                 l2.supporting_line())) return false;
      return (l1.lower_xyz_extremity() == l2.lower_xyz_extremity()) &&
             (l1.higher_xyz_extremity() == l2.higher_xyz_extremity());
    }

    template< class SK>
    bool
    non_oriented_equal( const typename SK::Circular_arc_3 &c1,
                        const typename SK::Circular_arc_3 &c2)
    {
      if(!non_oriented_equal<SK>(c1.supporting_circle(),
                                 c2.supporting_circle())) return false;
      if(c1.rep().is_full() && c2.rep().is_full())
        return true;
      return (c1.source() == c2.source()) &&
             (c1.target() == c2.target());
    }

   namespace internal {
     // we need to pass the result_type, to hack around the fact that
     // object doesn't support operator=(const T&) and so we keep backwards compatibility
     template<typename SK, typename RT>
     struct pair_transform {
       RT operator()(const std::pair< typename SK::Root_for_spheres_2_3, unsigned >& p) {
         return RT(std::make_pair(typename SK::Circular_arc_point_3(p.first), p.second));
       }
     };

     template<typename SK>
     struct pair_transform<SK, CGAL::Object> {
       CGAL::Object operator()(const std::pair< typename SK::Root_for_spheres_2_3, unsigned >& p) {
         return CGAL::make_object(std::make_pair(typename SK::Circular_arc_point_3(p.first), p.second));
       }
     };

   // obscure trick: calculate the regular intersection between two
   // objects, throw it away if empty, else check if the result was a
   // point_3, if so turn in into a (circular_arc_point, unsigned)
   // pair, otherwise just dump the value with conversion to RT
   // (again: converting to RT before assigning to the Iterator is
   // just to keep object working)
   template<typename SK, typename RT, typename OutputIterator>
   struct Point_conversion_visitor : public boost::static_visitor<OutputIterator> {
     Point_conversion_visitor(const OutputIterator& it) : it(it) {}
     template<typename T>
     OutputIterator operator()(const T& t) { *it++ = RT(t); return it; }

     OutputIterator operator()(const typename SK::Point_3& p) { 
       // 2 multiplicities
       *it++ = RT(std::make_pair(typename SK::Circular_arc_point_3(p), 2u));
       return it;
     }
     OutputIterator it;
   };

   } // namespace internal

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Sphere_3 & s, 
                const typename SK::Line_3 & l, 
	        OutputIterator res)
    {
      typedef typename SK3_Intersection_traits<SK, typename SK::Sphere_3, typename SK::Line_3>
        ::type result_type;
      typedef typename SK::Algebraic_kernel                          Algebraic_kernel;
      typedef typename SK::Polynomial_for_spheres_2_3  Equation_sphere; 
      typedef typename SK::Polynomials_for_line_3      Equation_line; 
      typedef typename SK::Root_for_spheres_2_3        Root_for_spheres_2_3;
      CGAL_kernel_precondition(!s.is_degenerate());
      CGAL_kernel_precondition(!l.is_degenerate());      
      Equation_sphere e1 = get_equation<SK>(s);
      Equation_line e2 = get_equation<SK>(l);
      typedef std::vector< std::pair < Root_for_spheres_2_3, unsigned > > 
        solutions_container;
      solutions_container solutions;
      Algebraic_kernel().solve_object()(e1, e2, std::back_inserter(solutions)); 

      return std::transform(solutions.begin(), solutions.end(), res, internal::pair_transform<SK, result_type>());
    }

    // The special 3 object functions
    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Sphere_3 & s1, 
                const typename SK::Sphere_3 & s2, 
	        const typename SK::Sphere_3 & s3, 
                OutputIterator res)
    {
       // circle_3, sphere_3, (circular_arc_point_3, uint)
       typedef typename SK::Polynomial_for_spheres_2_3  Equation_sphere;
       typedef typename SK::Root_for_spheres_2_3  Root_for_spheres_2_3;
       typedef typename SK::Circular_arc_point_3  Circular_arc_point_3;
       typedef typename SK::Circle_3  Circle_3;
       typedef typename SK::Point_3  Point_3;
       typedef typename SK::Sphere_3  Sphere_3;
       typedef typename SK::Algebraic_kernel  Algebraic_kernel;
       typedef typename SK3_Intersection_traits<SK, Sphere_3, Sphere_3, Sphere_3>::type result_type;

       CGAL_kernel_precondition(!s1.is_degenerate());
       CGAL_kernel_precondition(!s2.is_degenerate());
       CGAL_kernel_precondition(!s3.is_degenerate());
       if(non_oriented_equal<SK>(s1,s2) && non_oriented_equal<SK>(s2,s3)) {
         #if CGAL_INTERSECTION_VERSION < 2
         *res++ = make_object(s1);
         #else
         *res++ = result_type(s1);
         #endif
         return res;
       }
       if(non_oriented_equal<SK>(s1,s2)) {
         if(typename Intersection_traits<SK, Sphere_3, Sphere_3>::result_type v = 
            SK().intersect_3_object()(s1, s3)) {
           #if CGAL_INTERSECTION_VERSION < 2
           if( const Point_3* p = object_cast<Point_3>(&v) )
             *res++ = make_object(std::make_pair(Circular_arc_point_3(*p), 2u));
           else
             *res++ = v;
           #else
           internal::Point_conversion_visitor<SK, result_type, OutputIterator> visitor(res);
           return boost::apply_visitor(visitor,
             *v);
           #endif
         }
         return res;
       }
       if(non_oriented_equal<SK>(s1,s3) || non_oriented_equal<SK>(s2,s3)) {
         if(typename Intersection_traits<SK, Sphere_3, Sphere_3>::result_type v = 
            SK().intersect_3_object()(s1, s2)) {
           #if CGAL_INTERSECTION_VERSION < 2
           if( const Point_3* p = object_cast<Point_3>(&v) )
             *res++ = make_object(std::make_pair(Circular_arc_point_3(*p), 2u));
           else
             *res++ = v;
           #else
           internal::Point_conversion_visitor<SK, result_type, OutputIterator> visitor(res);
           return boost::apply_visitor(
             visitor,
             *v);
           #endif
         }
         return res;
       }
       if(SK().collinear_3_object()(s1.center(),s2.center(),s3.center())) {
         typename Intersection_traits<SK, Sphere_3, Sphere_3>::result_type v = 
           SK().intersect_3_object()(s1, s2);
         if(!v) return res;
         if(const Point_3* p = CGAL::internal::intersect_get<Point_3>(v)) {
            if(SK().has_on_3_object()(s3, *p)) {
              #if CGAL_INTERSECTION_VERSION < 2
              *res++ = make_object(std::make_pair(Circular_arc_point_3(*p),2u));
              #else
              *res++ = result_type(std::make_pair(Circular_arc_point_3(*p),2u));
              #endif
            }
             return res;
         }
         if(const Circle_3* c = CGAL::internal::intersect_get<Circle_3>(v)) {
            if(SK().has_on_3_object()(s3, *c)) {
              #if CGAL_INTERSECTION_VERSION < 2
              *res++ = make_object(*c);
              #else
              *res++ = result_type(*c);
              #endif
            }
           return res;
         }
         return res;
       }
       Equation_sphere e1 = get_equation<SK>(s1);
       Equation_sphere e2 = get_equation<SK>(s2);
       Equation_sphere e3 = get_equation<SK>(s3);
       typedef std::vector< std::pair < Root_for_spheres_2_3, unsigned > > 
         algebraic_solutions_container;
       algebraic_solutions_container solutions;
       Algebraic_kernel().solve_object()(e1, e2, e3, std::back_inserter(solutions)); 

       return std::transform(solutions.begin(), solutions.end(), res, internal::pair_transform<SK, result_type>());
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Plane_3 & p, 
                const typename SK::Sphere_3 & s1, 
	        const typename SK::Sphere_3 & s2, 
                OutputIterator res)
    {
      typedef typename boost::variant< std::pair< typename SK::Circular_arc_point_3, unsigned int>,
                                       typename SK::Circle_3 > result_type;
      typedef typename SK::Root_for_spheres_2_3  Root_for_spheres_2_3;
      typedef typename SK::Polynomial_for_spheres_2_3  Equation_sphere;
      typedef typename SK::Polynomial_1_3  Equation_plane;
      typedef typename SK::Plane_3  Plane_3;
      typedef typename SK::Sphere_3 Sphere_3;
      typedef typename SK::Algebraic_kernel  Algebraic_kernel;
      #if CGAL_INTERSECTION_VERSION < 2
      typedef typename SK::Circular_arc_point_3  Circular_arc_point_3;
      #endif
      CGAL_kernel_precondition(!p.is_degenerate());
      CGAL_kernel_precondition(!s1.is_degenerate());
      CGAL_kernel_precondition(!s2.is_degenerate());
      if(non_oriented_equal<SK>(s1,s2)) {
        if(typename Intersection_traits<SK, Plane_3, Sphere_3>::result_type v = 
            SK().intersect_3_object()(p, s1)) {
           #if CGAL_INTERSECTION_VERSION < 2
           if( const typename SK::Point_3* p = CGAL::object_cast<typename SK::Point_3>(&v) )
             *res++ = make_object(std::make_pair(Circular_arc_point_3(*p), 2u));
           else
             *res++ = v;
           #else
           internal::Point_conversion_visitor<SK, result_type, OutputIterator> visitor(res);
           return boost::apply_visitor(
             visitor,
             *v);
           #endif
         }
         return res;
      }
      Plane_3 radical_p = SK().construct_radical_plane_3_object()(s1,s2);
      if(non_oriented_equal<SK>(p,radical_p)) {
        if(typename Intersection_traits<SK, Plane_3, Sphere_3>::result_type v = 
            SK().intersect_3_object()(p, s1)) {
           #if CGAL_INTERSECTION_VERSION < 2
           if( const typename SK::Point_3* p = CGAL::object_cast<typename SK::Point_3>(&v) )
             *res++ = make_object(std::make_pair(Circular_arc_point_3(*p), 2u));
           else
             *res++ = v;
           #else
           internal::Point_conversion_visitor<SK, result_type, OutputIterator> visitor(res);
           return boost::apply_visitor(
             visitor,
             *v);
           #endif
         }
         return res;
      }
      Equation_sphere e1 = get_equation<SK>(s1);
      Equation_sphere e2 = get_equation<SK>(s2);
      Equation_plane e3 = get_equation<SK>(p);
      typedef std::vector< std::pair < Root_for_spheres_2_3, unsigned > > 
        algebraic_solutions_container;
      algebraic_solutions_container solutions;
      Algebraic_kernel().solve_object()(e1, e2, e3, std::back_inserter(solutions)); 
      return std::transform(solutions.begin(), solutions.end(), res, internal::pair_transform<SK, result_type>());
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Plane_3 & p1, 
                const typename SK::Plane_3 & p2, 
	        const typename SK::Sphere_3 & s, 
                OutputIterator res)
    { 
      typedef typename boost::variant< std::pair< typename SK::Circular_arc_point_3, unsigned int>, 
                                       typename SK::Circle_3 > result_type;
      typedef typename SK::Root_for_spheres_2_3  Root_for_spheres_2_3;
      typedef typename SK::Polynomial_for_spheres_2_3  Equation_sphere;
      typedef typename SK::Polynomial_1_3  Equation_plane;
      typedef typename SK::Algebraic_kernel  Algebraic_kernel;
      typedef typename SK::Plane_3 Plane_3;
      typedef typename SK::Sphere_3 Sphere_3;
      CGAL_kernel_precondition(!p1.is_degenerate());
      CGAL_kernel_precondition(!p2.is_degenerate());
      CGAL_kernel_precondition(!s.is_degenerate());      
      if(non_oriented_equal<SK>(p1,p2)) {
        if(typename Intersection_traits<SK, Plane_3, Sphere_3>::result_type v = 
            SK().intersect_3_object()(p1, s)) {
           #if CGAL_INTERSECTION_VERSION < 2
           typedef typename SK::Circular_arc_point_3  Circular_arc_point_3;
           if( const typename SK::Point_3* p = CGAL::object_cast<typename SK::Point_3>(&v) )
             *res++ = make_object(std::make_pair(Circular_arc_point_3(*p), 2u));
           else
             *res++ = v;
           #else
           internal::Point_conversion_visitor<SK, result_type, OutputIterator> visitor(res);
           return boost::apply_visitor(
             visitor,
             *v);
           #endif
         }
         return res;
      }
      Equation_plane e1 = get_equation<SK>(p1);
      Equation_plane e2 = get_equation<SK>(p2);
      Equation_sphere e3 = get_equation<SK>(s);
      typedef std::vector< std::pair < Root_for_spheres_2_3, unsigned > > 
        algebraic_solutions_container;
      algebraic_solutions_container solutions;
      Algebraic_kernel().solve_object()(e1, e2, e3, std::back_inserter(solutions)); 
      return std::transform(solutions.begin(), solutions.end(), res, internal::pair_transform<SK, result_type>());
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Circle_3 & c, 
                const typename SK::Plane_3 & p,
                OutputIterator res)
    { 
      return intersect_3<SK>(p,c.supporting_plane(),c.diametral_sphere(),res);
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Circle_3 & c, 
                const typename SK::Sphere_3 & s,
                OutputIterator res)
    { 
      return intersect_3<SK>(c.supporting_plane(),s,c.diametral_sphere(),res);
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Circle_3 & c1, 
                const typename SK::Circle_3 & c2,
                OutputIterator res)
    { 
      typedef typename SK::Root_for_spheres_2_3     Root_for_spheres_2_3;
      typedef typename SK::Polynomials_for_circle_3 Equation_circle;
      typedef typename SK::Algebraic_kernel         Algebraic_kernel;
      typedef typename SK::Circle_3                 Circle_3;

      typedef typename SK3_Intersection_traits<SK, Circle_3, Circle_3>
        ::type result_type;

      if(non_oriented_equal<SK>(c1,c2)) {
         *res++ = CGAL::internal::sk3_intersection_return<result_type>(c1);
         return res;
      }
      Equation_circle e1 = get_equation<SK>(c1);
      Equation_circle e2 = get_equation<SK>(c2);
      typedef std::vector< std::pair < Root_for_spheres_2_3, unsigned > > 
        algebraic_solutions_container;
      algebraic_solutions_container solutions;
      Algebraic_kernel().solve_object()(e1, e2, std::back_inserter(solutions)); 
      return std::transform(solutions.begin(), solutions.end(), res, internal::pair_transform<SK, result_type>());
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Circle_3 & c, 
                const typename SK::Line_3 & l,
                OutputIterator res)
    { 
      typedef typename SK::Root_for_spheres_2_3     Root_for_spheres_2_3;
      typedef typename SK::Polynomials_for_circle_3 Equation_circle;
      typedef typename SK::Polynomials_for_line_3   Equation_line;
      typedef typename SK::Circle_3                 Circle_3;

      typedef typename SK3_Intersection_traits<SK, Circle_3, typename SK::Line_3>
        ::type result_type;

      typedef typename SK::Algebraic_kernel  Algebraic_kernel;
      CGAL_kernel_precondition(!l.is_degenerate());
      Equation_circle e1 = get_equation<SK>(c);
      Equation_line e2 = get_equation<SK>(l);
      typedef std::vector< std::pair < Root_for_spheres_2_3, unsigned > > 
        algebraic_solutions_container;
      algebraic_solutions_container solutions;
      Algebraic_kernel().solve_object()(e1, e2, std::back_inserter(solutions)); 
      return std::transform(solutions.begin(), solutions.end(), res, internal::pair_transform<SK, result_type>());
    }

    // At the moment we dont need those functions
    // But in the future maybe (some make_x_monotone? etc..)
    template <class SK>
    typename SK::Circular_arc_point_3
    x_extremal_point(const typename SK::Sphere_3 & c, bool i)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().x_critical_points_object()(typename SK::Get_equation()(c),i);
    }

    template <class SK,class OutputIterator>
    OutputIterator
    x_extremal_points(const typename SK::Sphere_3 & c, OutputIterator res)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().x_critical_points_object()(typename SK::Get_equation()(c),res);
    }

    template <class SK>
    typename SK::Circular_arc_point_3
    y_extremal_point(const typename SK::Sphere_3 & c, bool i)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().y_critical_points_object()(typename SK::Get_equation()(c),i);
    }

    template <class SK,class OutputIterator>
    OutputIterator
    y_extremal_points(const typename SK::Sphere_3 & c, OutputIterator res)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().y_critical_points_object()(typename SK::Get_equation()(c),res);
    }

    template <class SK>
    typename SK::Circular_arc_point_3
    z_extremal_point(const typename SK::Sphere_3 & c, bool i)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().z_critical_points_object()(typename SK::Get_equation()(c),i);
    }

    template <class SK,class OutputIterator>
    OutputIterator
    z_extremal_points(const typename SK::Sphere_3 & c, OutputIterator res)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().z_critical_points_object()(typename SK::Get_equation()(c),res);
    }

    template <class SK>
    typename SK::Circular_arc_point_3
    x_extremal_point(const typename SK::Circle_3 & c, bool i)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().x_critical_points_object()(typename SK::Get_equation()(c),i);
    }

    template <class SK,class OutputIterator>
    OutputIterator
    x_extremal_points(const typename SK::Circle_3 & c, OutputIterator res)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().x_critical_points_object()(typename SK::Get_equation()(c),res);
    }

    template <class SK>
    typename SK::Circular_arc_point_3
    y_extremal_point(const typename SK::Circle_3 & c, bool i)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().y_critical_points_object()(typename SK::Get_equation()(c),i);
    }

    template <class SK,class OutputIterator>
    OutputIterator
    y_extremal_points(const typename SK::Circle_3 & c, OutputIterator res)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().y_critical_points_object()(typename SK::Get_equation()(c),res);
    }

    template <class SK>
    typename SK::Circular_arc_point_3
    z_extremal_point(const typename SK::Circle_3 & c, bool i)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().z_critical_points_object()(typename SK::Get_equation()(c),i);
    }

    template <class SK,class OutputIterator>
    OutputIterator
    z_extremal_points(const typename SK::Circle_3 & c, OutputIterator res)
    {
      typedef typename SK::Algebraic_kernel   Algebraic_kernel;
      return Algebraic_kernel().z_critical_points_object()(typename SK::Get_equation()(c),res);
    }

  }//SphericalFunctors
}//CGAL

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_SPHERE_3_H
