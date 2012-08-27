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
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado, 
//             Sebastien Loriot, Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_3_H
#define CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_3_H

#include <CGAL/Circular_kernel_3/internal_functions_on_circular_arc_3.h>
#include <boost/tuple/tuple.hpp>


namespace CGAL {
  namespace internal{
    template < class SK >
    class Circular_arc_3 {
      typedef typename SK::Plane_3              Plane_3;
      typedef typename SK::Circle_3             Circle_3;
      typedef typename SK::Sphere_3             Sphere_3;
      typedef typename SK::Point_3              Point_3;
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef typename SK::Line_3               Line_3;
      typedef typename SK::FT                   FT;

    private:

      typedef boost::tuple<Circle_3, Circular_arc_point_3, 
                               Circular_arc_point_3> Rep;
      typedef typename SK::template Handle<Rep>::type Base;

      Base base;
      mutable bool _full;
      // It is the sign of the cross product 
      // of the vector (Center -> S) x (Center -> T)
      // it saves execution time for the has_on functor
      Sign _sign_cross_product;

    public:

      const Sphere_3& reference_sphere(){
        return get_ref_sphere(get(base).template get<0>());
      };

        
      Circular_arc_3()
      {}

      // The arc of circle goes from s to t in counterclockwise orientation
      // in relation to the normal vector N of the supporting plane
      // such that 
      // if N.x != 0 then N.x > 0
      // else if N.y != 0 -> N.y > 0
      // else N.z > 0
      // Interesting thing is that if _sign_cross_product is negative
      // the arc is the bigger one (angle > pi)
      Circular_arc_3(const Circle_3 &c, 
                     const Circular_arc_point_3 &s,
                     const Circular_arc_point_3 &t) 
      : _full(false)
      {
        // l must pass through s and t, and s != t
        CGAL_kernel_precondition(SK().has_on_3_object()(c,s));
        CGAL_kernel_precondition(SK().has_on_3_object()(c,t));
        CGAL_kernel_precondition(s != t);
        base = Rep(c,s,t);
        // we can optimize the computations of the sign (for the has_on functor), 
        // by computing the vector s-c and t-s, in order to use them directly on 
        // another compute_sign_of_cross_product function
        // we can save time computing the substractions
        // the problem is: more memory space is needed
        _sign_cross_product =
          CGAL::SphericalFunctors::compute_sign_of_cross_product<SK>(s,t,c.center());
      }

      Circular_arc_3(const Circle_3 &c, 
                     const Point_3 &s,
                     const Circular_arc_point_3 &t) 
      : _full(false)
      {
        // l must pass through s and t, and s != t
        CGAL_kernel_precondition(SK().has_on_3_object()(c,s));
        CGAL_kernel_precondition(SK().has_on_3_object()(c,t));
        CGAL_kernel_precondition(Circular_arc_point_3(s) != t);
        base = Rep(c,s,t);
        _sign_cross_product = 
          CGAL::SphericalFunctors::compute_sign_of_cross_product<SK>(s,t,c.center());
      }

      Circular_arc_3(const Circle_3 &c, 
                     const Circular_arc_point_3 &s,
                     const Point_3 &t) 
      : _full(false)
      {
        // l must pass through s and t, and s != t
        CGAL_kernel_precondition(SK().has_on_3_object()(c,s));
        CGAL_kernel_precondition(SK().has_on_3_object()(c,t));
        CGAL_kernel_precondition(s != Circular_arc_point_3(t));
        base = Rep(c,s,t);
        _sign_cross_product = 
          CGAL::SphericalFunctors::compute_sign_of_cross_product<SK>(s,t,c.center());
      }

      Circular_arc_3(const Circle_3 &c, 
                     const Point_3 &s,
                     const Point_3 &t) 
      : _full(false)
      {
        // l must pass through s and t, and s != t
        CGAL_kernel_precondition(SK().has_on_3_object()(c,s));
        CGAL_kernel_precondition(SK().has_on_3_object()(c,t));
        CGAL_kernel_precondition(Circular_arc_point_3(s) != 
                                 Circular_arc_point_3(t));
        base = Rep(c,s,t);
        _sign_cross_product = 
          CGAL::SphericalFunctors::compute_sign_of_cross_product<SK>(s,t,c.center());
      }

      // This is the one of the two cases we want that s == t
      // that makes the is_full() correct and complete
      Circular_arc_3(const Circle_3 &c)
      : _full(true)
      {
        const Plane_3 &p = c.supporting_plane();
        if(is_zero(p.b()) && is_zero(p.c())) {
          const Circular_arc_point_3 v = 
	    SphericalFunctors::y_extremal_point<SK>(c,true);
          base = Rep(c,v,v);
        } else {
          const Circular_arc_point_3 v = 
	    SphericalFunctors::x_extremal_point<SK>(c,true);
          base = Rep(c,v,v);
        }
        /* don't matter
        _sign_cross_product = 0;
        */
      }

      // This is the second case where we want that s == t
      // that makes the is_full() correct and complete
      Circular_arc_3(const Circle_3 &c,const Circular_arc_point_3& point)
      : base(Rep(c,point,point)),_full(true)
      {CGAL_kernel_precondition(SK().has_on_3_object()(c,point));}
      
      Circular_arc_3(const Circle_3 &c, 
                     const Sphere_3 &s1, bool less_xyz_s1,
                     const Sphere_3 &s2, bool less_xyz_s2) 
      {
         std::vector<Object> sols1, sols2;
         // The spheres must not include the circle
         CGAL_kernel_precondition(!SK().has_on_3_object()(s1,c));
         CGAL_kernel_precondition(!SK().has_on_3_object()(s2,c));
         SK().intersect_3_object()(c, s1, std::back_inserter(sols1));
         SK().intersect_3_object()(c, s2, std::back_inserter(sols2));
         // l must intersect s1 and s2
         CGAL_kernel_precondition(sols1.size() > 0);
         CGAL_kernel_precondition(sols2.size() > 0);
         const std::pair<typename SK::Circular_arc_point_3, unsigned>& pair1=
            *object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(
              &sols1[(sols1.size()==1)?(0):(less_xyz_s1?0:1)]
            );
         const std::pair<typename SK::Circular_arc_point_3, unsigned>& pair2=
            *object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(
              &sols2[(sols2.size()==1)?(0):(less_xyz_s2?0:1)]
            );        
         // the source and target must be different
         CGAL_kernel_precondition(pair1.first != pair2.first);
         *this = Circular_arc_3(c, pair1.first, pair2.first);
      }

      Circular_arc_3(const Circle_3 &c, 
                     const Plane_3 &p1, bool less_xyz_p1,
                     const Plane_3 &p2, bool less_xyz_p2) 
      {
         std::vector<Object> sols1, sols2;
         // The planes must not include the circle
         CGAL_kernel_precondition(!SK().has_on_3_object()(p1,c));
         CGAL_kernel_precondition(!SK().has_on_3_object()(p2,c));
         SK().intersect_3_object()(c, p1, std::back_inserter(sols1));
         SK().intersect_3_object()(c, p2, std::back_inserter(sols2));
         // l must intersect s1 and s2
         CGAL_kernel_precondition(sols1.size() > 0);
         CGAL_kernel_precondition(sols2.size() > 0);
         const std::pair<typename SK::Circular_arc_point_3, unsigned>& pair1=
            *object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(
              &sols1[(sols1.size()==1)?(0):(less_xyz_p1?0:1)]
            );
         const std::pair<typename SK::Circular_arc_point_3, unsigned>& pair2=
            *object_cast<std::pair<typename SK::Circular_arc_point_3, unsigned> >(
              &sols2[(sols2.size()==1)?(0):(less_xyz_p2?0:1)]
            );                
         // the source and target must be different
         CGAL_kernel_precondition(pair1.first != pair2.first);
         *this = Circular_arc_3(c, pair1.first, pair2.first);
      }

      Circular_arc_3(const Point_3 &begin,
		     const Point_3 &middle,
		     const Point_3 &end)
      {
	CGAL_kernel_precondition(!CGAL::collinear(begin, middle, end));
	const Circle_3 c = Circle_3(begin, middle, end);
	base = Rep(c,begin,end);
        _sign_cross_product =
          CGAL::SphericalFunctors::compute_sign_of_cross_product<SK>
	         (begin,end,c.center());
      }

      const Circle_3& supporting_circle() const 
      {
        return get(base).template get<0>();
      }

      const Circular_arc_point_3& source() const 
      {
        return get(base).template get<1>();
      }

      const Circular_arc_point_3& target() const 
      {
        return get(base).template get<2>();
      }

      Plane_3 supporting_plane() const {
        return supporting_circle().supporting_plane();
      }

      Point_3 center() const {
        return supporting_circle().center();
      }

      FT squared_radius() const {
        return supporting_circle().squared_radius();
      }

      Sphere_3 diametral_sphere() const {
        return supporting_circle().diametral_sphere();
      }

      bool is_full() const {
        return _full;
      }

      Sign sign_cross_product() const {
        return _sign_cross_product;
      }

      double approximate_angle() const {
        if(is_full()) return 2.0*CGAL_PI;
        const double x1 = to_double(source().x());
        const double y1 = to_double(source().y());
        const double z1 = to_double(source().z());
        const double x2 = to_double(target().x());
        const double y2 = to_double(target().y());
        const double z2 = to_double(target().z());
        const double dx = x2-x1;
        const double dy = y2-y1;
        const double dz = z2-z1;
        const double d_sq = dx*dx + dy*dy + dz*dz;
        const double r_sq = to_double(squared_radius());
        const double ap_ang = 2.0 * std::asin(0.5 * std::sqrt(d_sq / r_sq));
        if(sign_cross_product() == NEGATIVE) return 2.0 * CGAL_PI - ap_ang;
        else return ap_ang;
      }

      double approximate_squared_length() const {
        const double ang = approximate_angle();
        return ang * ang * to_double(squared_radius());
      }

      // It is of course possible to increase the precision
      // maybe it will be done after
      CGAL::Bbox_3 bbox() const {
        return supporting_circle().bbox();
      }

      bool operator==(const Circular_arc_3 &) const;
      bool operator!=(const Circular_arc_3 &) const;

    };

    template < class SK >
    CGAL_KERNEL_INLINE
    bool
    Circular_arc_3<SK>::operator==(const Circular_arc_3<SK> &t) const
    {
      if (CGAL::identical(base, t.base))
        return true;		
      return CGAL::SphericalFunctors::non_oriented_equal<SK>(*this, t);
    }

    template < class SK >
    CGAL_KERNEL_INLINE
    bool
    Circular_arc_3<SK>::operator!=(const Circular_arc_3<SK> &t) const
    {
      return !(*this == t);
    }

  }
}

#endif
