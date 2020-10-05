// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado,
//             Sebastien Loriot, Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a
// STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCLE_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCLE_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/Circle_type.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_plane_3.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_circular_arc_point_3.h>
#include <CGAL/Circular_kernel_3/Intersection_traits.h>
#include <CGAL/Root_of_traits.h>

namespace CGAL {
  namespace SphericalFunctors {

    template < class SK >
    typename SK::Circle_3
    construct_circle_3(const typename SK::Polynomials_for_circle_3 &eq)
    {
      typedef typename SK::Point_3 Point_3;
      typedef typename SK::Plane_3 Plane_3;
      typedef typename SK::Circle_3 Circle_3;
      typedef typename SK::Sphere_3 Sphere_3;
      typedef typename SK::FT FT;
      Sphere_3 s = construct_sphere_3<SK>(eq.first);
      Plane_3 p = construct_plane_3<SK>(eq.second);
      const FT d2 = CGAL::square(p.a()*s.center().x() +
                                 p.b()*s.center().y() +
                                 p.c()*s.center().z() + p.d()) /
       (CGAL::square(p.a()) + CGAL::square(p.b()) + CGAL::square(p.c()));
      // We do not accept circles with radius 0 (should we?)
      CGAL_kernel_precondition(d2 <  s.squared_radius());
      // d2 < s.squared_radius()
      Point_3 center = p.projection(s.center());
      return Circle_3(center,s.squared_radius() - d2,p);
    }


    template < class SK >
    CGAL::Circle_type
    classify_circle_3(const typename SK::Circle_3& circle,const typename SK::Sphere_3& sphere)
    {
      typedef typename SK::Algebraic_kernel::Root_for_spheres_2_3   Root_for_spheres_2_3;
      typedef typename SK::Circular_arc_point_3                     Circular_arc_point_3;

      CGAL_kernel_precondition(SK().has_on_3_object()(sphere,circle));

      //if circle is a great circle, it can only be a bipolar or a threaded.
      if (circle.center()==sphere.center()){
        if (circle.supporting_plane().orthogonal_vector().z()==0)
          return CGAL::BIPOLAR;
        return CGAL::THREADED;
      }

      typename SK::Root_of_2 radius=CGAL::make_sqrt(sphere.squared_radius());
      Circular_arc_point_3 north_pole( Root_for_spheres_2_3(sphere.center().x(),sphere.center().y(),sphere.center().z()+radius) );
      Circular_arc_point_3 south_pole( Root_for_spheres_2_3(sphere.center().x(),sphere.center().y(),sphere.center().z()-radius) );


      const typename SK::Sphere_3& supp_sphere=circle.diametral_sphere();
      typename SK::Bounded_side_3 bounded_side=SK().bounded_side_3_object();


      CGAL::Bounded_side side_of_north=bounded_side(supp_sphere,north_pole);
      CGAL::Bounded_side side_of_south=bounded_side(supp_sphere,south_pole);

      if (side_of_north==ON_BOUNDARY || side_of_south==ON_BOUNDARY)
        return CGAL::POLAR;

      if (side_of_north==ON_BOUNDED_SIDE || side_of_south==ON_BOUNDED_SIDE)
        return CGAL::THREADED;

      CGAL_kernel_precondition(side_of_north==ON_UNBOUNDED_SIDE && side_of_south==ON_UNBOUNDED_SIDE);

      return CGAL::NORMAL;
    }


    template < class SK >
    inline typename SK::FT
    extremal_points_z_coordinate(const typename SK::Circle_3& circle,const typename SK::Sphere_3& sphere)
    {
      CGAL_kernel_precondition(SK().has_on_3_object()(sphere,circle));
      CGAL_kernel_precondition(classify_circle_3<SK>(circle,sphere)==CGAL::NORMAL);

      const typename SK::Point_3& circle_center=circle.center();
      const typename SK::Point_3& sphere_center=sphere.center();

      return
      typename SK::FT(2) * (circle_center-sphere_center).z() * sphere.squared_radius()
               / ( SK().compute_squared_distance_3_object()(circle_center,sphere_center) + sphere.squared_radius()-circle.squared_radius() )
               + sphere_center.z();
    }

    template < class SK, class Output_iterator >
    Output_iterator theta_extremal_points(const typename SK::Circle_3& circle,const typename SK::Sphere_3& sphere,Output_iterator out_it){
      CGAL_kernel_precondition(classify_circle_3<SK>(circle,sphere)==NORMAL);
      CGAL_kernel_precondition(SK().has_on_3_object()(sphere,circle));

      typename SK::FT z_coord=extremal_points_z_coordinate<SK>(circle,sphere);

      typename SK::Plane_3 plane(0,0,1,-z_coord);
      std::vector<typename SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Plane_3>::type > inters;

      intersect_3<SK>(circle,plane,std::back_inserter(inters));
      CGAL_kernel_precondition(inters.size()==2);
      const std::pair<typename SK::Circular_arc_point_3,unsigned>* pt[2]={nullptr,nullptr};
      pt[0]=CGAL::Intersections::internal::intersect_get<std::pair<typename SK::Circular_arc_point_3,unsigned> >(inters[0]);
      pt[1]=CGAL::Intersections::internal::intersect_get<std::pair<typename SK::Circular_arc_point_3,unsigned> >(inters[1]);
      CGAL_kernel_precondition(pt[0]!=nullptr);
      CGAL_kernel_precondition(pt[1]!=nullptr);

      if ( compare_theta_of_pts<SK>(pt[0]->first,pt[1]->first,sphere) == SMALLER){
        *out_it++=pt[0]->first;
        *out_it++=pt[1]->first;
      }
      else{
        *out_it++=pt[1]->first;
        *out_it++=pt[0]->first;
      }

      return out_it;
    }

    template < class SK >
    typename SK::Circular_arc_point_3 theta_extremal_point(const typename SK::Circle_3& circle,const typename SK::Sphere_3& sphere,bool is_smallest){
      typename SK::Circular_arc_point_3 pts[2];
      theta_extremal_points(circle,sphere,pts);
      if (is_smallest)
        return pts[0];
      return pts[1];
    }

    template < class SK,class Output_iterator>
    Output_iterator make_circle_theta_monotone(const typename SK::Circle_3& circle,const typename SK::Sphere_3& sphere,Output_iterator out_it){
      CGAL::Circle_type type=classify_circle_3<SK>(circle,sphere);
      switch (type){
        case THREADED:
        {
          *out_it++=typename SK::Circular_arc_3(circle);
          break;
        }
        case POLAR:{
          typename SK::Vector_3 ortho=circle.center()-sphere.center();
          CGAL_kernel_precondition(ortho.z()!=0);
          bool is_north_pole=ortho.z()>0;
          typename SK::Root_of_2 radius = (is_north_pole?1:-1)* make_sqrt(sphere.squared_radius());
          typename SK::Circular_arc_point_3 source_target(
            typename SK::Algebraic_kernel::Root_for_spheres_2_3(
              sphere.center().x(),
              sphere.center().y(),
              sphere.center().z()+radius
            )
          );
          *out_it++=typename SK::Circular_arc_3(circle,source_target);
          break;
        }
        case NORMAL:{
          typename SK::Circular_arc_point_3 ints[2];
          theta_extremal_points(circle,sphere,ints);
          *out_it++=typename SK::Circular_arc_3(circle,ints[0],ints[1]);
          *out_it++=typename SK::Circular_arc_3(circle,ints[1],ints[0]);
        }
        break;
        case BIPOLAR:
          CGAL_kernel_precondition_msg(false, "This function does not accept bipolar circle as input.");
      }
      return out_it;
    }



  }//SphericalFunctors
}//CGAL

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCLE_3_H
