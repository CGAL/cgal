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

#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_POINT_3_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_POINT_3_H

namespace CGAL {
  namespace SphericalFunctors {

    template< class SK>
    bool
      equal( const typename SK::Circular_arc_point_3 &p1,
             const typename SK::Circular_arc_point_3 &p2)
    {
      return p1.rep() == p2.rep();
    }
        

  //compute the half_quadrant relative to the cylindrical coordinate system defined by sphere of pt
  template<class SK>
  float half_quadrant(const typename SK::Root_of_2& Rx,
                      const typename SK::Root_of_2& Ry)
  {
    int x=CGAL_NTS sign(Rx);
    int y=CGAL_NTS sign(Ry);
    if (y>0){
      if (x>0)  switch (CGAL_NTS sign(Rx-Ry)){case -1: return 2; case 0: return 1.5; case 1: return 1; }; 
      if (x<0) switch (CGAL::opposite(CGAL_NTS sign(Rx+Ry))){case -1: return 3; case 0: return 3.5; case 1: return 4; };
      return 2.5;//OPTI : we have more information here by x_y
    }
    else{
      if (y<0){
        if (x>0) switch (CGAL_NTS sign(Rx+Ry)){case -1: return 7; case 0: return 7.5; case 1: return 8; };
        if (x<0) switch (CGAL::opposite(CGAL_NTS sign(Rx-Ry))){case -1: return 6; case 0: return 5.5; case 1: return 5; };
        return 6.5;//OPTI : we have more information here by x_y
      }
      else{
        if (x>0)  return 0.5;
        if (x<0) return 4.5;
        return 0;
      }
    }
  }    
    

  template<class SK>
  inline
  float half_quadrant(const typename SK::Circular_arc_point_3& pt,
                      const typename SK::Sphere_3& sphere)
  {
    typename SK::Root_of_2 Rx=pt.x()-sphere.center().x();
    typename SK::Root_of_2 Ry=pt.y()-sphere.center().y();
    return half_quadrant<SK>(Rx,Ry);
  }
  
  template<class SK>
  bool is_hquadrant_a_tangent(float hq){
    if (hq >7 || hq<2 || (hq>3 && hq<6))
      return true;
    return false;
  }
  
  template <class SK>
  CGAL::Comparison_result compare_theta_of_pts(const typename SK::Circular_arc_point_3& pt1,
                                               const typename SK::Circular_arc_point_3& pt2,
                                               const typename SK::Sphere_3& sphere)
  {
    CGAL_kernel_precondition( (pt1.x()!=sphere.center().x() || pt1.y()!=sphere.center().x()) && 
                              (pt2.x()!=sphere.center().x() || pt2.y()!=sphere.center().x()) );
    float hq1=half_quadrant<SK>(pt1,sphere);
    float hq2=half_quadrant<SK>(pt2,sphere);
    CGAL::Sign res=CGAL_NTS sign(hq1-hq2);
    if (res == CGAL::EQUAL){//same quadrant 
      if (floor(hq1)!=hq1)//for hquadrant boundary
        res=CGAL::EQUAL;
      else{
        //compare tan or cot, same expression due to constant sign within a hquadrant
        //and different monotonicity type of tan and cot
        if ( is_hquadrant_a_tangent<SK>(hq1) )
          res= CGAL_NTS compare( (pt1.y() - sphere.center().y()) / (pt1.x() - sphere.center().x()) ,
                              (pt2.y() - sphere.center().y()) / (pt2.x() - sphere.center().x()) ); 
        else
          res= CGAL_NTS compare( (pt2.x() - sphere.center().x()) / (pt2.y() - sphere.center().y()) ,
                              (pt1.x() - sphere.center().x()) / (pt1.y() - sphere.center().y()) );          
      }
    }
    return res;
  }

  template <class SK>
  CGAL::Comparison_result compare_theta_pt_vector(const typename SK::Circular_arc_point_3& pt,
                                                  const typename SK::Vector_3& v,
                                                  const typename SK::Sphere_3& sphere)
  {
    CGAL_kernel_precondition( (pt.x()!=sphere.center().x() || pt.y()!=sphere.center().x()) && 
                              (v.x()!=0 || v.y()!=0) &&  v.z()==0);
    float hq1=half_quadrant<SK>(pt,sphere);
    float hq2=half_quadrant<SK>(v.x(),v.y());
    CGAL::Sign res=CGAL_NTS sign(hq1-hq2);
    if (res == CGAL::EQUAL){//same quadrant 
      if (floor(hq1)!=hq1)//for hquadrant boundary
        res=CGAL::EQUAL;
      else{
        //compare tan or cot, same expression due to constant sign within a hquadrant
        //and different monotonicity type of tan and cot
        res= CGAL_NTS compare( v.x() * (pt.y() - sphere.center().y()) ,
                            v.y() * (pt.x() - sphere.center().x())  ); 
      }
    }
    return res;
  }

  template <class SK>
  CGAL::Comparison_result compare_theta_vectors(const typename SK::Vector_3& m1,
                                                const typename SK::Vector_3& m2
                                               )
  {
    CGAL_kernel_precondition( (m1.x()!=0 || m1.y()!=0) &&  m1.z()==0 &&  
                             (m2.x()!=0 || m2.y()!=0) &&  m2.z()==0 );
    float hq1=half_quadrant<SK>(m1.x(),m1.y());
    float hq2=half_quadrant<SK>(m2.x(),m2.y());
    CGAL::Sign res=CGAL_NTS sign(hq1-hq2);
    if (res == CGAL::EQUAL){//same quadrant 
      if (floor(hq1)!=hq1)//for hquadrant boundary
        res=CGAL::EQUAL;
      else{
        //compare tan or cot, same expression due to constant sign within a hquadrant
        //and different monotonicity type of tan and cot
        res= CGAL_NTS compare( m2.x() * m1.y() ,
                            m2.y() * m1.x() ); 
      }
    }
    return res;
  }
  
  template <class SK>
  CGAL::Comparison_result compare_theta_z(const typename SK::Circular_arc_point_3& pt1,
                                          const typename SK::Circular_arc_point_3& pt2,
                                          const typename SK::Sphere_3& sphere,
                                          bool decreasing_z=false)
  {
    CGAL::Comparison_result res=compare_theta_of_pts<SK>(pt1,pt2,sphere);
    if (res==CGAL::EQUAL)
      res = (decreasing_z?(CGAL::opposite(SK().compare_z_3_object()(pt1,pt2))):(SK().compare_z_3_object()(pt1,pt2)));
    return res;
  }  
    
  }//SphericalFunctors
}//CGAL

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_ON_CIRCULAR_ARC_POINT_3_H
