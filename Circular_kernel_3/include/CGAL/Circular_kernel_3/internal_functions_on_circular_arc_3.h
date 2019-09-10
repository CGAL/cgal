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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado, 
//             Sebastien Loriot, Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/Circular_kernel_3/internal_function_has_on_spherical_kernel.h>
#include <CGAL/Circular_kernel_3/internal_functions_on_circle_3.h>

namespace CGAL {
  namespace SphericalFunctors {

    template< class SK>
    bool
    equal( const typename SK::Circular_arc_3 &c1,
           const typename SK::Circular_arc_3 &c2)
    {
      return c1.rep() == c2.rep();
    }

    template <class SK>
    inline
    bool
    do_overlap(const typename SK::Circular_arc_3 &c1,
               const typename SK::Circular_arc_3 &c2,
               const bool known_equal_supporting_circle = false)
    { 
      if(!known_equal_supporting_circle) {
        if(!non_oriented_equal<SK>(c1.supporting_circle(), 
                                   c2.supporting_circle()))
          return false;
      }
      if(c1.rep().is_full()) return true;
      if(c2.rep().is_full()) return true;
      if((SK().has_on_3_object()(c1,c2.target(),true)) || 
         (SK().has_on_3_object()(c1,c2.source(),true))) return true;
      return SK().has_on_3_object()(c2,c1.source(),true);
    }

    template < class SK >
    void
    split(const typename SK::Circular_arc_3 &c,
	  const typename SK::Circular_arc_point_3 &p,
	  typename SK::Circular_arc_3 &c1,
	  typename SK::Circular_arc_3 &c2)
    {
      // The point must be on the circular arc 
      CGAL_kernel_precondition(SK().has_on_3_object()(c, p));
      typedef typename SK::Circular_arc_3  Circular_arc_3;
      // It doesn't make sense to split an arc on an extremity
      CGAL_kernel_precondition(c.source() != p);
      CGAL_kernel_precondition(c.target() != p);
      const Circular_arc_3 &rc1 = 
        Circular_arc_3(c.supporting_circle(), c.source(), p);
      const Circular_arc_3 &rc2 = 
        Circular_arc_3(c.supporting_circle(), p, c.target());
      if ( SK().compare_xyz_3_object()(rc1.source(), rc2.source()) != 
           SMALLER) {
        c1 = rc2; c2 = rc1;
      } else { c1 = rc1; c2 = rc2; }
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Line_3 & l, 
                const typename SK::Circular_arc_3 & ca, 
	        OutputIterator res)
    { 
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<
        typename SK3_Intersection_traits<SK, typename SK::Line_3, typename SK::Circle_3>::type 
        > solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      solutions_container solutions;

      SK().intersect_3_object()(l, ca.supporting_circle(),
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
        const Solution& sol=*CGAL::internal::intersect_get<Solution>(solutions[0]);
         if(SK().has_on_3_object()(ca,sol.first,true))
           *res++ = solutions[0];
      } else {
         const Solution& sol1=*CGAL::internal::intersect_get<Solution>(solutions[0]);
         const Solution& sol2=*CGAL::internal::intersect_get<Solution>(solutions[1]);        
         if(SK().has_on_3_object()(ca,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(ca,sol2.first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Circle_3 & c, 
                const typename SK::Circular_arc_3 & ca, 
	        OutputIterator res)
    { 
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef typename SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Circular_arc_3>::type result_type;

      typedef std::vector<typename SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Circle_3
        >::type > solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      if(non_oriented_equal<SK>(c, ca.supporting_circle())) {
        *res++ = result_type(ca);
      }

      solutions_container solutions;

      SK().intersect_3_object()(ca.supporting_circle(), c, 
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         const Solution& sol=*CGAL::internal::intersect_get<Solution>(solutions[0]);
         if(SK().has_on_3_object()(ca,sol.first,true))
           *res++ = solutions[0];
      } else {
         const Solution& sol1=*CGAL::internal::intersect_get<Solution>(solutions[0]);
         const Solution& sol2=*CGAL::internal::intersect_get<Solution>(solutions[1]);        
         if(SK().has_on_3_object()(ca,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(ca,sol2.first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Sphere_3 & s,
                const typename SK::Circular_arc_3 & c,
	       OutputIterator res)
    {
      typedef typename SK3_Intersection_traits<SK, typename SK::Sphere_3, typename SK::Circular_arc_3
        >::type result_type;

      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<typename SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Sphere_3>::type
      > solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      if(SK().has_on_3_object()(s, c.supporting_circle())) {
        *res++ = result_type(c);
      }

      solutions_container solutions;

      SK().intersect_3_object()(c.supporting_circle(), s,
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         const Solution& sol=*CGAL::internal::intersect_get<Solution>(solutions[0]);
         if(SK().has_on_3_object()(c,sol.first,true))
           *res++ = solutions[0];
      } else {
         const Solution& sol1=*CGAL::internal::intersect_get<Solution>(solutions[0]);        
         const Solution& sol2=*CGAL::internal::intersect_get<Solution>(solutions[1]);        
         if(SK().has_on_3_object()(c,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(c,sol2.first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Plane_3 & p, 
                const typename SK::Circular_arc_3 & ca, 
	        OutputIterator res)
    {
      typedef typename SK3_Intersection_traits<SK, typename SK::Plane_3, typename SK::Circular_arc_3
        >::type result_type;
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<
        typename SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Plane_3>::type
      > solutions_container;

      typedef std::pair<Circular_arc_point_3, unsigned> Solution;
      if(SK().has_on_3_object()(p,ca.supporting_circle())) {
        *res++ = CGAL::internal::sk3_intersection_return<result_type>(ca);
      }
      solutions_container solutions;

      SK().intersect_3_object()(ca.supporting_circle(), p,
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         const Solution& sol=*CGAL::internal::intersect_get<Solution>(solutions[0]);
         if(SK().has_on_3_object()(ca,sol.first,true))
           *res++ = CGAL::internal::sk3_intersection_return<result_type>(sol);
      } else {
         const Solution& sol1=*CGAL::internal::intersect_get<Solution>(solutions[0]);
         const Solution& sol2=*CGAL::internal::intersect_get<Solution>(solutions[1]);        
         if(SK().has_on_3_object()(ca,sol1.first,true))
           *res++ = CGAL::internal::sk3_intersection_return<result_type>(sol1);
         if(SK().has_on_3_object()(ca,sol2.first,true))
           *res++ = CGAL::internal::sk3_intersection_return<result_type>(sol2);
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Line_arc_3 & la, 
                const typename SK::Circular_arc_3 & ca, 
	        OutputIterator res)
    { 
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef std::vector<
        typename SK3_Intersection_traits<SK, typename SK::Line_3, typename SK::Line_3>::type> 
        solutions_container;
      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      solutions_container solutions;

      SK().intersect_3_object()(la.supporting_line(), ca.supporting_circle(),
                                std::back_inserter(solutions) );
      if(solutions.size() == 0) return res;
      if(solutions.size() == 1) {
         const Solution& sol=*CGAL::internal::intersect_get<Solution>(solutions[0]);
         if(SK().has_on_3_object()(ca,sol.first,true) &&
            SK().has_on_3_object()(la,sol.first,true))
           *res++ = solutions[0];
      } else {
         const Solution& sol1=*CGAL::internal::intersect_get<Solution>(solutions[0]);
         const Solution& sol2=*CGAL::internal::intersect_get<Solution>(solutions[1]);        
         if(SK().has_on_3_object()(ca,sol1.first,true) &&
            SK().has_on_3_object()(la,sol1.first,true))
           *res++ = solutions[0];
         if(SK().has_on_3_object()(ca,sol2.first,true) &&
            SK().has_on_3_object()(la,sol2.first,true))
           *res++ = solutions[1];
      }
      return res;
    }

    template < class SK, class OutputIterator >
    OutputIterator
    intersect_3(const typename SK::Circular_arc_3 & a1, 
                const typename SK::Circular_arc_3 & a2, 
	        OutputIterator res)
    { 
      typedef typename SK3_Intersection_traits<SK, typename SK::Circular_arc_3, typename SK::Circular_arc_3
        >::type result_type;
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef typename SK::Circular_arc_3 Circular_arc_3;
      typedef std::vector< typename SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Circle_3>
        ::type> solutions_container;

      typedef std::pair<Circular_arc_point_3, unsigned> Solution;

      if(non_oriented_equal<SK>(a1.supporting_circle(), a2.supporting_circle())) {
        if(a1.rep().is_full()) {
          *res++ = CGAL::internal::sk3_intersection_return<result_type>(a2); 
          //return res;
        }
        else if(a2.rep().is_full()) {
          *res++ = CGAL::internal::sk3_intersection_return<result_type>(a1); 
          //return res;
        } else {
          bool t2_in_a1 = SK().has_on_3_object()(a1,a2.target(),true);
          bool s2_in_a1 = SK().has_on_3_object()(a1,a2.source(),true);
          if(t2_in_a1 && s2_in_a1) {
            bool t1_in_a2 = SK().has_on_3_object()(a2,a1.target(),true);
            bool s1_in_a2 = SK().has_on_3_object()(a2,a1.source(),true);
            if(t1_in_a2 && s1_in_a2) {
              const Comparison_result comp = 
                SK().compare_xyz_3_object()(a1.source(), a2.source());
              if(comp < 0) {
                if(a1.source() == a2.target()) {
                  *res++ = CGAL::internal::sk3_intersection_return<result_type>(std::make_pair(a1.source(),1u));
                } else {
                  const Circular_arc_3 & arc =
	          Circular_arc_3(a1.supporting_circle(),a1.source(),a2.target());
	          *res++ = CGAL::internal::sk3_intersection_return<result_type>(arc);
                }
                if(a2.source() == a1.target()) {
                  *res++ = CGAL::internal::sk3_intersection_return<result_type>(std::make_pair(a2.source(),1u));
                } else {
                  const Circular_arc_3 & arc =
	          Circular_arc_3(a1.supporting_circle(),a2.source(),a1.target());
	          *res++ = CGAL::internal::sk3_intersection_return<result_type>(arc);
                }
              } else if(comp > 0) {
                if(a2.source() == a1.target()) {
                  *res++ = CGAL::internal::sk3_intersection_return<result_type>(std::make_pair(a2.source(),1u));
                } else {
                  const Circular_arc_3 & arc =
	          Circular_arc_3(a1.supporting_circle(),a2.source(),a1.target());
	          *res++ = CGAL::internal::sk3_intersection_return<result_type>(arc);
                }
                if(a1.source() == a2.target()) {
                  *res++ = CGAL::internal::sk3_intersection_return<result_type>(std::make_pair(a1.source(),1u));
                } else {
                  const Circular_arc_3 & arc =
	          Circular_arc_3(a1.supporting_circle(),a1.source(),a2.target());
	          *res++ = CGAL::internal::sk3_intersection_return<result_type>(arc);
                } 
              } else { 
                *res++ = CGAL::internal::sk3_intersection_return<result_type>(a1);
              }
            } else {
              *res++ = CGAL::internal::sk3_intersection_return<result_type>(a2);
            }
          } else if(t2_in_a1) {
            if(a1.source() == a2.target()) 
              *res++ = CGAL::internal::sk3_intersection_return<result_type>(std::make_pair(a1.source(),1u));
            else {
              const Circular_arc_3 & arc =
	        Circular_arc_3(a1.supporting_circle(),a1.source(),a2.target());
	      *res++ = CGAL::internal::sk3_intersection_return<result_type>(arc);
            } //return res;
          } else if(s2_in_a1) {
            if(a2.source() == a1.target()) {
              *res++ = CGAL::internal::sk3_intersection_return<result_type>(std::make_pair(a2.source(),1u));
            } else {
              const Circular_arc_3 & arc =
	        Circular_arc_3(a1.supporting_circle(),a2.source(),a1.target());
	      *res++ = CGAL::internal::sk3_intersection_return<result_type>(arc);
            }
          } else if(SK().has_on_3_object()(a2,a1.source(),true)) {
              *res++ = CGAL::internal::sk3_intersection_return<result_type>(a1);
          } 
        }
      } else {
        solutions_container solutions;

        SK().intersect_3_object()(a1.supporting_circle(), a2.supporting_circle(), 
                                  std::back_inserter(solutions) );
        if(solutions.size() == 0) return res;
        if(solutions.size() == 1) {
          const Solution& sol=*CGAL::internal::intersect_get<Solution>(solutions[0]);
          if(SK().has_on_3_object()(a1,sol.first,true) &&
             SK().has_on_3_object()(a2,sol.first,true))
            *res++ = solutions[0];
        } else {
          const Solution& sol1=*CGAL::internal::intersect_get<Solution>(solutions[0]);
          const Solution& sol2=*CGAL::internal::intersect_get<Solution>(solutions[1]);          
          if(SK().has_on_3_object()(a1,sol1.first,true) &&
             SK().has_on_3_object()(a2,sol1.first,true))
            *res++ = solutions[0];
          if(SK().has_on_3_object()(a1,sol2.first,true) &&
             SK().has_on_3_object()(a2,sol2.first,true))
            *res++ = solutions[1];
        }
      }
      return res;
    }

    template < class SK>
    bool
    is_theta_monotone_3(const typename SK::Circular_arc_3 & arc,const typename SK::Sphere_3& sphere)
    {
      CGAL_kernel_precondition(SK().has_on_3_object()(sphere,arc));
      CGAL::Circle_type type=classify_circle_3<SK>(arc.supporting_circle(),sphere);
      CGAL_kernel_precondition(type!=CGAL::BIPOLAR);
      if (type==THREADED)
        return true;
      if (type==POLAR){
        bool circle_contains_north = arc.supporting_circle().center().z() > sphere.center().z();
        typename SK::Root_of_2 radius=make_sqrt(sphere.squared_radius());
        typename SK::Circular_arc_point_3 pole (
          typename SK::Algebraic_kernel::Root_for_spheres_2_3( sphere.center().x(),
                                                              sphere.center().y(),
                                                              sphere.center().z()+(circle_contains_north?1:-1)*radius
          )
        );
        
        if (arc.source().z()==pole.z() || arc.target().z()==pole.z())
          return true;
        
        return !has_on<SK>(arc,pole);
      }
      
      
      typename SK::FT z_coord=extremal_points_z_coordinate<SK>(arc.supporting_circle(),sphere);
      typename SK::Plane_3 plane(0,0,1,-z_coord);
      std::vector<typename SK3_Intersection_traits<SK, typename SK::Plane_3, typename SK::Circular_arc_3>
        ::type> inters;
      
      intersect_3<SK>(plane,arc,std::back_inserter(inters));
      
      if (inters.empty())
        return true;
      
      //checks whether circular arc endpoints are theta extremal
      unsigned nb_extrem = (arc.source().z()-z_coord == 0)? 1:0;
      if (arc.target().z()-z_coord == 0) ++nb_extrem;
      
      if (inters.size()==nb_extrem)
        return true;
      
      return false;
    }
    
    template < class SK,class Output_iterator>
    Output_iterator 
    make_circular_arc_theta_monotone( const typename SK::Circular_arc_3& arc,
                                      const typename SK::Sphere_3& sphere,
                                      Output_iterator out_it)
    {
      CGAL::Circle_type type=classify_circle_3<SK>(arc.supporting_circle(),sphere);
      CGAL_kernel_precondition(type!=BIPOLAR);
      switch (type){
        case THREADED:
        case POLAR:{
          bool circle_contains_north = arc.supporting_circle().center().z() > sphere.center().z();
          typename SK::Root_of_2 radius=make_sqrt(sphere.squared_radius());
          typename SK::Circular_arc_point_3 pole (
            typename SK::Algebraic_kernel::Root_for_spheres_2_3( sphere.center().x(),
                                                                sphere.center().y(),
                                                                sphere.center().z()+(circle_contains_north?1:-1)*radius
            )
          );
          
          if (arc.source().z()==pole.z() || arc.target().z()==pole.z())
            *out_it++=arc;
          else{
            if ( has_on<SK>(arc,pole) ){
              *out_it++=typename SK::Circular_arc_3(arc.supporting_circle(),arc.source(),pole);
              *out_it++=typename SK::Circular_arc_3(arc.supporting_circle(),pole,arc.target());
            }
            else
              *out_it++=arc;
          }
        }
        break;
        case NORMAL:{
          typename SK::FT z_coord=extremal_points_z_coordinate<SK>(arc.supporting_circle(),sphere);
          typename SK::Plane_3 plane(0,0,1,-z_coord);
          std::vector<typename SK3_Intersection_traits< SK, typename SK::Plane_3, typename SK::Circular_arc_3 >::type
          > inters;
          
          intersect_3<SK>(plane,arc,std::back_inserter(inters));
          
          //No intersection with horizontal plane: theta-monotone
          if (inters.empty()){
            *out_it++=arc;
            break;
          }
          
          //check if endpoints of circular arc are theta extremal points
          unsigned nb_extrem = (arc.source().z()-z_coord == 0)? 1:0;
          if (arc.target().z()-z_coord == 0) ++nb_extrem;
          
          if (inters.size()==nb_extrem){
            *out_it++=arc;
            break;
          }

          //one endpoint is extremal: just split the arc
          if (nb_extrem==1){
            const std::pair<typename SK::Circular_arc_point_3,unsigned>* pt[2]={NULL,NULL};
            pt[0]=CGAL::internal::intersect_get<std::pair<typename SK::Circular_arc_point_3,unsigned> >(inters[0]);
            pt[1]=CGAL::internal::intersect_get<std::pair<typename SK::Circular_arc_point_3,unsigned> >(inters[1]);
            CGAL_kernel_precondition(pt[0]!=NULL);
            CGAL_kernel_precondition(pt[1]!=NULL);
            const typename SK::Circular_arc_point_3& midpt=(arc.source()==pt[0]->first || arc.target()==pt[0]->first)?pt[1]->first:pt[0]->first;
            *out_it++=typename SK::Circular_arc_3(arc.supporting_circle(),arc.source(),midpt);
            *out_it++=typename SK::Circular_arc_3(arc.supporting_circle(),midpt,arc.target());
            break;
          }
          
          CGAL_kernel_precondition(nb_extrem==0);
          
          //only one intersection points
          if (inters.size()==1){
            const std::pair<typename SK::Circular_arc_point_3,unsigned>* midpt=NULL;
            midpt=CGAL::internal::intersect_get<std::pair<typename SK::Circular_arc_point_3,unsigned> >(inters[0]);
            CGAL_kernel_precondition(midpt!=NULL);
            *out_it++=typename SK::Circular_arc_3(arc.supporting_circle(),arc.source(),midpt->first);
            *out_it++=typename SK::Circular_arc_3(arc.supporting_circle(),midpt->first,arc.target());
            break;
          }
          
          //three arcs are defined by two intersection points
          const std::pair<typename SK::Circular_arc_point_3,unsigned>* pt[2]={NULL,NULL};
          pt[0]=CGAL::internal::intersect_get<std::pair<typename SK::Circular_arc_point_3,unsigned> >(inters[0]);
          pt[1]=CGAL::internal::intersect_get<std::pair<typename SK::Circular_arc_point_3,unsigned> >(inters[1]);
          CGAL_kernel_precondition(pt[0]!=NULL);
          CGAL_kernel_precondition(pt[1]!=NULL);
          
          typename SK::Circular_arc_3 arc1=typename SK::Circular_arc_3(arc.supporting_circle(),arc.source(),pt[0]->first);
          typename SK::Circular_arc_3 arc2=typename SK::Circular_arc_3(arc.supporting_circle(),pt[0]->first,arc.target());
          if ( has_on<SK>(arc1,pt[1]->first) ){
            *out_it++=typename SK::Circular_arc_3(arc1.supporting_circle(),arc1.source(),pt[1]->first);
            *out_it++=typename SK::Circular_arc_3(arc1.supporting_circle(),pt[1]->first,arc1.target());            
            *out_it++=arc2;
          }
          else{
            *out_it++=arc1;            
            *out_it++=typename SK::Circular_arc_3(arc2.supporting_circle(),arc2.source(),pt[1]->first);
            *out_it++=typename SK::Circular_arc_3(arc2.supporting_circle(),pt[1]->first,arc2.target());            
          }
          break;
        }
        case BIPOLAR:
          CGAL_kernel_precondition(!"This function does not accept bipolar circle as input.");
      }
      
      return out_it;
    }
    
    //this function indicates whether a theta monotone circular arc is 
    //included in the upper part of its supporting circle (wrt theta extremal pts)
    template <class SK>
    bool 
    is_upper_arc(const typename SK::Circular_arc_3& arc,
                 const typename SK::Sphere_3& sphere)
    {
      CGAL_kernel_precondition(is_theta_monotone_3<SK>(arc,sphere));
      CGAL::Circle_type type=classify_circle_3<SK>(arc.supporting_circle(),sphere);
      CGAL_kernel_precondition(type!=THREADED && type !=BIPOLAR);
      
      //case of POLAR circle
      if (type==POLAR)
        return ( ( arc.supporting_circle().center().z()-sphere.center().z() ) < 0 );
      
      //case of NORMAL circle
      typename SK::FT z_coord=extremal_points_z_coordinate<SK>(arc.supporting_circle(),sphere);
      if ( z_coord < arc.source().z() ) return true;
      if ( z_coord > arc.source().z() ) return false;
      if ( z_coord < arc.target().z() ) return true;
      if ( z_coord > arc.target().z() ) return false;
      
      //source and target are theta extremal points
      
      typename SK::Point_3 out_sphere_normal=CGAL::ORIGIN + ( arc.supporting_circle().center()-sphere.center() );
      typename SK::Point_3 ori(0,0,0);
      int wise=(out_sphere_normal > ori)? 1:-1; //check if seen ccw from the side going from the sphere center to the circle center
      
      typename SK::Vector_3 cir_center(arc.supporting_circle().center().x()-sphere.center().x(),arc.supporting_circle().center().y()-sphere.center().y(),0);
      CGAL::Comparison_result source_vs_center=compare_theta_pt_vector<SK>(arc.source(),cir_center,sphere);
      CGAL::Comparison_result target_vs_center=compare_theta_pt_vector<SK>(arc.target(),cir_center,sphere);
      
      
      if (source_vs_center==target_vs_center){//circle is cut by meridian at theta=0
        CGAL::Comparison_result source_vs_target=compare_theta_of_pts<SK>(arc.source(),arc.target(),sphere);  
        wise*=(source_vs_target==SMALLER)?-1:1;
      }
      else
        wise*=(source_vs_center==SMALLER)?1:-1;
      
      return wise==1?false:true;
    }

    //compute the intersection points of a theta monotone arc with a meridian.
    //precondition: the intersection point always exists.
    template <class SK>
    typename SK::Circular_arc_point_3 
    intersect_3( const typename SK::Circular_arc_3& arc,
                 const typename SK::Vector_3&m,
                 const typename SK::Sphere_3& sphere)
    {
      typename SK::Plane_3 plane(sphere.center(),sphere.center()+m,sphere.center()+typename SK::Vector_3(0,0,1));

      std::vector<typename SK3_Intersection_traits<SK, typename SK::Plane_3, typename SK::Circular_arc_3 >
        ::type> inters;
      intersect_3<SK>(plane,arc,std::back_inserter(inters));
      CGAL_kernel_precondition(!inters.empty());
      if (inters.size()==1){
          const typename SK::Circular_arc_point_3& pt=
            CGAL::internal::intersect_get<std::pair<typename SK::Circular_arc_point_3,unsigned> >(inters[0])->first;
          return pt;
      }
      
      CGAL_kernel_precondition(classify_circle_3<SK>(arc.supporting_circle(),sphere)!=NORMAL);
      
      const typename SK::Circular_arc_point_3& pts1 = 
        CGAL::internal::intersect_get<std::pair<typename SK::Circular_arc_point_3,unsigned> >(inters[0])->first;
      const typename SK::Circular_arc_point_3& pts2 =
        CGAL::internal::intersect_get<std::pair<typename SK::Circular_arc_point_3,unsigned> >(inters[1])->first;
      
      
      //either a polar (1 pole + 1 pt) or a threaded circle (2 pts with theta-coord = +/- pi)
      if ( (pts1.x()!=0 && pts1.y()!=0) && compare_theta_pt_vector<SK>(pts1,m,sphere)==CGAL::EQUAL)
        return pts1;
      else
        return pts2;
    }
      
    
    template <class SK>
    //Compare the z coordinates of the intersection points with the meridian defined by m of two theta monotone arcs.
    CGAL::Comparison_result 
    compare_z_at_theta_arcs( const typename SK::Circular_arc_3& arc1,
                             const typename SK::Circular_arc_3& arc2,
                             const typename SK::Vector_3& m,
                             const typename SK::Sphere_3& sphere) 
    {
      CGAL_kernel_precondition(is_theta_monotone_3<SK>(arc1,sphere));
      CGAL_kernel_precondition(is_theta_monotone_3<SK>(arc2,sphere));
      
      typename SK::Circular_arc_point_3 pt1=intersect_3<SK>(arc1,m,sphere);
      typename SK::Circular_arc_point_3 pt2=intersect_3<SK>(arc2,m,sphere);
      return CGAL::compare(pt1.z(),pt2.z());
    }

    template <class SK>
    CGAL::Comparison_result 
    compare_z_at_theta_pt_arc(const typename SK::Circular_arc_point_3& pt,
                              const typename SK::Circular_arc_3& arc,
                              const typename SK::Sphere_3& sphere)
    {
      CGAL_kernel_precondition(is_theta_monotone_3<SK>(arc,sphere));
      CGAL::Circle_type type=classify_circle_3<SK>(arc.supporting_circle(),sphere);
      CGAL_kernel_precondition(type!=BIPOLAR);
      
      switch(type){
        case THREADED:
        case POLAR:
        {
          const typename SK::Plane_3& plane=arc.supporting_circle().supporting_plane();
          typename SK::Vector_3 ortho=plane.orthogonal_vector();
          int res=typename SK::Algebraic_kernel().sign_at_object()(
              typename SK::Algebraic_kernel::Polynomial_1_3(ortho.x(),ortho.y(),ortho.z(),plane.d()),
              pt.coordinates()
          );
          res*=(plane.orthogonal_vector().z()>0)?1:-1;
          return CGAL::Comparison_result(res);
        }
        default:
        {
          typename SK::Vector_3 ortho=arc.supporting_circle().center()-sphere.center();
          typename SK::FT d=-ortho * typename SK::Vector_3(CGAL::ORIGIN,arc.supporting_circle().center());
          int res=typename SK::Algebraic_kernel().sign_at_object()(
              typename SK::Algebraic_kernel::Polynomial_1_3(ortho.x(),ortho.y(),ortho.z(),d),
              pt.coordinates()
          );
          if (res==0) return CGAL::EQUAL;
          if (res==-1) return CGAL::compare(pt.z(),extremal_points_z_coordinate<SK>(arc.supporting_circle(),sphere));
          return is_upper_arc<SK>(arc,sphere)?CGAL::SMALLER:CGAL::LARGER;
        }
      }
    }
    
  }//SphericalFunctors
}//CGAL

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_3_H
