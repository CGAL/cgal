// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//             Pedro Machado    <tashimir@gmail.com>
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_SPHERICAL_KERNEL_CIRCLE_3_H
#define CGAL_SPHERICAL_KERNEL_CIRCLE_3_H

#include <CGAL/Circular_kernel_3/internal_functions_on_sphere_3.h>
#include <CGAL/Sphere_with_radius_3.h>
#include <boost/utility/enable_if.hpp>

namespace CGAL {
  namespace CGALi{
    
    //default implementation suppose we have two arbitrary spheres
    template <class T1,class T2,class SK>
    class Circle_representation_3{};
      
    template <class SK>
    class Circle_representation_3<CGAL::Sphere_with_radius_3<SK>,CGAL::Sphere_with_radius_3<SK>,SK>{
      public:
      typedef typename CGAL::Sphere_with_radius_3<SK> Sphere_3;
      private:
      typedef typename SK::Plane_3      Plane_3;
      typedef typename SK::Point_3      Point_3;
      typedef typename SK::Vector_3     Vector_3;
      typedef typename SK::Direction_3  Direction_3;
      typedef typename SK::FT           FT;  
      
      typedef std::pair<Sphere_3,Sphere_3> Rep;
      typedef typename SK::template Handle<Rep>::type  Base;

      Base base;  
      public:
      Circle_representation_3() {}
      Circle_representation_3(const Sphere_3& S1,const Sphere_3& S2):base(S1,S2){}
      
      Circle_representation_3(const Point_3& center, const FT& squared_r, const Direction_3& d){CGAL_precondition(false);}//use enable_if or trick like this to invalidate them
      Circle_representation_3(const Point_3& center, const FT& squared_r, const Vector_3& normal){CGAL_precondition(false);}
      Circle_representation_3(const Point_3& center, const FT& squared_r, const Plane_3& p){CGAL_precondition(false);}
      Circle_representation_3(const Plane_3 &p,const Sphere_3 &s){CGAL_precondition(false);}        
        
      Plane_3 supporting_plane() const {
        return SK().construct_radical_plane_3_object()(get(base).first,get(base).second);
      }

      const Sphere_3& supporting_sphere() const {
        return get(base).first;
      }
      
      const Sphere_3& reference_sphere() const {
        return get(base).second;
      }      
      
      //~ Point_3 center() const {
        //~ Plane_3 p=supporting_plane();
        //~ return p.projection(get(base).first.center());
      //~ }
      
      FT get_circle_center_coeff(const Point_3& scenter,const FT& r2) const{
        //~ return (((FT) 0.5) +  (reference_sphere().squared_radius() - r2)/(FT)(2* (my_pow(scenter.x(),2) 
        //~ + my_pow(scenter.y(),2) + my_pow(scenter.z(),2) )) );    
        return ((FT) 0.5) +  (reference_sphere().squared_radius() - r2)/(2* (scenter-CGAL::ORIGIN).squared_length());    
      }
 
      Point_3 center() const {
        FT coeff=get_circle_center_coeff(supporting_sphere().center(),supporting_sphere().squared_radius());
        return CGAL::ORIGIN+(supporting_sphere().center()-CGAL::ORIGIN)*coeff;    
      };
      
      FT squared_radius() const {
        const Point_3& c=center();
        //~ Point_3 p=get(base).first.center();
        //~ return get(base).first.squared_radius()-CGAL::squared_distance(p,c);
        return reference_sphere().squared_radius()-(c-CGAL::ORIGIN).squared_length();
      }

      Sphere_3 diametral_sphere() const {
        return Sphere_3(center(),squared_radius());
      }  
    };
      
    //specialization using one plane and one diametral sphere
    template <class SK>
    class Circle_representation_3<typename CGAL::Sphere_3<SK>,typename CGAL::Plane_3<SK>,SK >{
      public:
      typedef typename CGAL::Sphere_3<SK> Sphere_3;
      private:
      typedef typename SK::Plane_3      Plane_3;
      typedef typename SK::Point_3      Point_3;
      typedef typename SK::Vector_3     Vector_3;
      typedef typename SK::Direction_3  Direction_3;
      typedef typename SK::FT           FT;  

      
      typedef std::pair<Sphere_3,typename SK::Plane_3> Rep;
      typedef typename SK::template Handle<Rep>::type  Base;

      Base base;  
      public:

      Circle_representation_3() {}
        
      Circle_representation_3(const Point_3& center, const FT& squared_r, const Direction_3& d) 
      {
        // It is not allowed non-positive radius 
        // Should we keep this pre-condition?
        CGAL_kernel_assertion(squared_r > 0);
        base = Rep(Sphere_3(center,squared_r), 
                   plane_from_point_direction(center, d));
      }

      Circle_representation_3(const Point_3& center, const FT& squared_r, const Vector_3& normal) 
      {
        // It is not allowed non-positive radius 
        // Should we keep this pre-condition?
        CGAL_kernel_assertion(squared_r > 0);
        base = Rep(Sphere_3(center,squared_r),plane_from_point_direction(center, normal.direction()));
      }

      Circle_representation_3(const Point_3& center, const FT& squared_r, const Plane_3& p)
      {
        // the plane contains the center
        CGAL_kernel_assertion((p.a() * center.x() +
                               p.b() * center.y() +
                               p.c() * center.z() +
                               p.d()) == CGAL::ZERO);
        // It is not allowed non-positive radius 
        // Should we keep this pre-condition?
        CGAL_kernel_assertion(squared_r > 0);
        base = Rep(Sphere_3(center,squared_r), p);
      }

      Circle_representation_3(const Sphere_3 &s1, 
               const Sphere_3 &s2) {
         std::vector<Object> sols;
         SK().intersect_3_object()(s1, s2, std::back_inserter(sols));
         // s1,s2 must intersect
         CGAL_kernel_precondition(sols.size() != 0);
         typename SK::Circle_3 circle;
         // the intersection must be a circle (no point allowed)
         CGAL_kernel_precondition(assign(circle,sols[0]));
         assign(circle,sols[0]);
         base=Rep(circle.diametral_sphere(),circle.supporting_plane());
         //~ *this = circle.rep();
      }

      Circle_representation_3(const Plane_3 &p, 
               const Sphere_3 &s) {
         std::vector<Object> sols;
         SK().intersect_3_object()(p, s, std::back_inserter(sols));
         // s1,s2 must intersect
         CGAL_kernel_precondition(sols.size() != 0);
         typename SK::Circle_3 circle;
         // the intersection must be a circle (no point allowed)
         CGAL_kernel_precondition(assign(circle,sols[0]));
         assign(circle,sols[0]);
         base=Rep(circle.diametral_sphere(),circle.supporting_plane());
      }
        
      const Plane_3& supporting_plane() const {
        return get(base).second;
      }
      
      const Sphere_3& supporting_sphere() const {
        return diametral_sphere();
      }
      
      Point_3 center() const {
        return diametral_sphere().center();
      }

      FT squared_radius() const {
        return diametral_sphere().squared_radius();
      }

      const Sphere_3& diametral_sphere() const {
        return get(base).first;
      }
    };

    //~ template <class SK,class Sphere= typename CGAL::Sphere_3<SK> > class Circle_3 {
    template <class SK,class Container=Circle_representation_3<typename CGAL::Sphere_3<SK>,typename CGAL::Plane_3<SK>,SK > >
    //~ template <class SK,class Container=Circle_representation_3<typename SK::Sphere_3,typename CGAL::Plane_3<SK>,SK > >
    class Circle_3 {    

      typedef typename SK::Plane_3      Plane_3;
      typedef typename Container::Sphere_3         Sphere_3;
      typedef typename SK::Point_3      Point_3;
      typedef typename SK::Vector_3     Vector_3;
      typedef typename SK::Direction_3  Direction_3;
      typedef typename SK::FT           FT;

    private:
      //~ typedef std::pair<Sphere_3, Plane_3>  Rep;
      //~ typedef typename SK::template Handle<Rep>::type  Base;

      //~ Base base;
      Container base;
    
    public:
      
 
    
      Circle_3() {}
      Circle_3(const Point_3& center, const FT& squared_r, const Direction_3& d):base(center,squared_r,d){}
      Circle_3(const Point_3& center, const FT& squared_r, const Vector_3& normal):base(center,squared_r,normal){}
      Circle_3(const Point_3& center, const FT& squared_r, const Plane_3& p):base(center,squared_r,p){}
      Circle_3(const Plane_3 &p,const Sphere_3 &s):base(p,s){}
      Circle_3(const Sphere_3 &s1, const Sphere_3 &s2):base(s1,s2){}


      const Plane_3& supporting_plane() const { return base.supporting_plane(); }
      const Sphere_3& supporting_sphere() const { return base.supporting_sphere(); }
      Point_3 center() const { return base.center(); }
      FT squared_radius() const { return base.squared_radius(); }
      const Sphere_3& diametral_sphere() const { return base.diametral_sphere(); }

      FT area_divided_by_pi() const {
        return squared_radius();
      }

      FT squared_length_divided_by_pi_square() const {
        return 4 * squared_radius();
      }

      static double pi;

      double approximate_area() const {
        return pi * to_double(squared_radius());
      }

      double approximate_squared_length() const {
        return pi * pi * 4.0 * to_double(squared_radius());
      }

      // this bbox function
      // can be optimize by doing different cases
      // for each variable = 0 (cases with is_zero)
      CGAL::Bbox_3 bbox() const {
        typedef CGAL::Interval_nt<false> Interval;
        CGAL::Interval_nt<false>::Protector ip;
        const Plane_3 &plane = supporting_plane();
        const Sphere_3 &s = diametral_sphere();
        const Point_3 &p = s.center();
        const FT &sq_r = s.squared_radius();
        const Interval a = CGAL::to_interval(plane.a());
        const Interval b = CGAL::to_interval(plane.b());
        const Interval c = CGAL::to_interval(plane.c());
        const Interval x = CGAL::to_interval(p.x());
        const Interval y = CGAL::to_interval(p.y());
        const Interval z = CGAL::to_interval(p.z());
        const Interval r2 = CGAL::to_interval(sq_r);
        const Interval r = CGAL::sqrt(r2); // maybe we can work with r2
                                           // in order to save this operation
                                           // but if the coefficients are to high
                                           // the multiplication would lead to inf results
        const Interval a2 = CGAL::square(a);
        const Interval b2 = CGAL::square(b);
        const Interval c2 = CGAL::square(c);
        const Interval sqr_sum = a2 + b2 + c2;
        const Interval mx = r * CGAL::sqrt((sqr_sum - a2)/sqr_sum);
        const Interval my = r * CGAL::sqrt((sqr_sum - b2)/sqr_sum);
        const Interval mz = r * CGAL::sqrt((sqr_sum - c2)/sqr_sum);
        return CGAL::Bbox_3((x-mx).inf(),(y-my).inf(), (z-mz).inf(),
                            (x+mx).sup(),(y+my).sup(), (z+mz).sup());
      }

      bool operator==(const Circle_3 &) const;
      bool operator!=(const Circle_3 &) const;

    };

    template < class SK, class Sphere >
    double Circle_3<SK,Sphere>::pi = CGAL_PI;
  //#PB : operateur for Circle_3 and Circle_on_sphere_3
    template < class SK, class Sphere >
    CGAL_KERNEL_INLINE
    bool
    Circle_3<SK,Sphere>::operator==(const Circle_3<SK,Sphere> &t) const
    {
      if (CGAL::identical(base, t.base))
        return true;
      return CGAL::SphericalFunctors::non_oriented_equal<SK>(*this, t);
    } 

    template < class SK, class Sphere >
    CGAL_KERNEL_INLINE
    bool
    Circle_3<SK,Sphere>::operator!=(const Circle_3<SK,Sphere> &t) const
    {
      return !(*this == t);
    }

  }
}

#endif
