// Copyright (c) 2005-2009  INRIA Sophia-Antipolis (France).
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
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Sebastien Loriot


#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_TO_RIGHT_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_TO_RIGHT_H

#include <CGAL/license/Circular_kernel_3.h>


namespace CGAL {
  namespace SphericalFunctors {

//Traits providing primitives to compare two theta-monotone circular arc to the right of a common point Pt
template <class SK>
struct Trait_for_cmp_tgt{
  typename SK::Algebraic_kernel::Root_for_spheres_2_3 Pt_;
  typedef typename SK::Algebraic_kernel::Root_of_2 Tk_type;
  
  Trait_for_cmp_tgt(const typename SK::Algebraic_kernel::Root_for_spheres_2_3& Pt,
                    const typename SK::Sphere_3& sphere)
                  :Pt_(Pt.x()-sphere.center().x(),Pt.y()-sphere.center().y(),Pt.z()-sphere.center().z()){}

  typename SK::Algebraic_kernel::Root_of_2 
  unsigned_tkz_coeff_normal(const typename SK::Point_3& C,const typename SK::FT& gamma_k) const
  {
    return CGAL_NTS sign(gamma_k)*(C.x()*Pt_.y()-C.y()*Pt_.x());
  }

  Tk_type 
  unsigned_tkz_coeff_threaded(const typename SK::Point_3& C) const {
    return C.x()*Pt_.y()-C.y()*Pt_.x();
  }
  
};

// Same as Trait_for_cmp_tgt but in the case where common point is at theta=0
template <class SK>
struct Trait_for_cmp_tgt_theta_0{
  typename SK::Algebraic_kernel::Root_for_spheres_2_3 Pt_;
  typedef typename SK::FT Tk_type;
  
  Trait_for_cmp_tgt_theta_0(const typename SK::Algebraic_kernel::Root_for_spheres_2_3& Pt,
                            const typename SK::Sphere_3& sphere)
                          :Pt_(Pt.x()-sphere.center().x(),Pt.y()-sphere.center().y(),Pt.z()-sphere.center().z()){}
  
  typename SK::FT 
  unsigned_tkz_coeff_normal( const typename SK::Point_3& C,const typename SK::FT& gamma_k) const
  {
    return - CGAL_NTS sign(gamma_k)*C.y();
  }

  Tk_type
  unsigned_tkz_coeff_threaded(const typename SK::Point_3& C) const{
    return -C.y();
  }
};

// Class to compare two theta-monotone circular arc to the right of a common point.
// This class need a traits class which can be either Trait_for_cmp_tgt or Trait_for_cmp_tgt_theta_0
//depending on the properties of the point.
template<class SK,class Traits>
struct Compare_to_right_of_arcs{
  Traits traits_;
  const typename SK::Sphere_3& sphere_;
  
  Compare_to_right_of_arcs(const Traits& traits,const typename SK::Sphere_3& sphere):
    traits_(traits),sphere_(sphere){}
  
  typename SK::FT 
  gamma_k ( const typename  SK::Point_3& c,
            const typename SK::FT& rk,
            const typename SK::FT& R,
            typename SK::FT& ak2) const
  {
    ak2=c.x()*c.x()+c.y()*c.y()+c.z()*c.z();
    return (ak2+R-rk)/(2.*ak2);
  }

  typename SK::FT
  gamma_k (const typename SK::Point_3& c,
            const typename SK::FT& rk,
            const typename SK::FT& R) const
  {
    typename SK::FT ak2=c.x()*c.x()+c.y()*c.y()+c.z()*c.z();
    return (ak2+R-rk)/(2.*ak2);
  }
  
  bool is_arc_an_upper_one(const typename SK::Circular_arc_3& arc,const typename SK::Sphere_3& sphere) const {
    CGAL::Circle_type type=classify_circle_3<SK>(arc.supporting_circle(),sphere);
    if ( type==CGAL::NORMAL)
      return  extremal_points_z_coordinate<SK>(arc.supporting_circle(),sphere)-sphere.center().z() < traits_.Pt_.z();
    CGAL_kernel_precondition(type==POLAR);
    return arc.supporting_circle().center().z() < sphere.center().z() ;
  }
    
  
  typename SK::FT square_norm_tk_normal (const typename SK::Point_3& c,const typename SK::FT& rk,const typename SK::FT& R) const;
  typename SK::FT square_norm_tk_threaded (const typename SK::Point_3& c,const typename SK::FT& rk,const typename SK::FT& R) const;
  void fill_tzk_n(const typename SK::Circular_arc_3& arc,typename Traits::Tk_type& tz,typename SK::FT& n,bool is_supporting_circle_threaded) const;
  int sign_of_delta(const typename SK::Circular_arc_3& arc1,bool circle1_threaded,const typename SK::Circular_arc_3& arc2,bool circle2_threaded) const;
  typename SK::FT give_rk(const typename SK::Circular_arc_3& arc) const;
  int compare_for_delta_eq_0_threaded(const typename SK::Circular_arc_3& arc_threaded,const typename SK::Circular_arc_3& arc,bool circle_threaded) const;
  int compare_for_delta_eq_0(const typename SK::Circular_arc_3& arc1,bool circle1_threaded,const typename SK::Circular_arc_3& arc2,bool circle2_threaded) const;
  CGAL::Comparison_result operator()(const typename SK::Circular_arc_3& arc1,const typename SK::Circular_arc_3& arc2,bool do_it_to_left=false) const; 
};


template<class SK,class Traits>
typename SK::FT
Compare_to_right_of_arcs<SK,Traits>::square_norm_tk_normal(
  const typename SK::Point_3& c,
  const typename SK::FT& rk,
  const typename SK::FT& R) const
{
  typename SK::FT ak2;
  typename SK::FT gk=gamma_k(c,rk,R,ak2);
  return ak2 * (R - gk * gk * ak2);
}

template<class SK,class Traits>
typename SK::FT 
Compare_to_right_of_arcs<SK,Traits>::square_norm_tk_threaded(
  const typename SK::Point_3& c,
  const typename SK::FT& rk,
  const typename SK::FT& R) const
{
  typename SK::FT ak2;
  typename SK::FT gk=gamma_k(c,rk,R,ak2);
  return ak2 * (R - gk * gk * ak2);
}

template<class SK,class Traits>
void
Compare_to_right_of_arcs<SK,Traits>::fill_tzk_n(
  const typename SK::Circular_arc_3& arc,
  typename Traits::Tk_type& tz,
  typename SK::FT& n,
  bool is_supporting_circle_threaded) const
{
  if (is_supporting_circle_threaded){
    //use a pseudo-sphere whose center is translated by ortho from the center
    // of the circle
    //~ tz=traits_.unsigned_tkz_coeff_threaded(arc.supporting_circle().supporting_sphere_center());
    typename SK::Vector_3 ortho=arc.supporting_circle().supporting_plane().orthogonal_vector();
    tz=traits_.unsigned_tkz_coeff_threaded(CGAL::ORIGIN+ortho);
    //~ tz*=(CGAL::pole_covered_by_supporting_sphere<SK>(arc.supporting_circle())==CGAL::NPOLE)?(1):(-1);
    if (ortho.z() <= 0)
      tz=-tz;
    //~ tz*= ( (ortho.z() > 0)?(1):(-1) );
    //~ n=square_norm_tk_threaded(arc.supporting_circle().supporting_sphere_center(),
                              //~ arc.supporting_circle().supporting_sphere_squared_radius(),
                              //~ sphere_.squared_radius());
    n=square_norm_tk_threaded(CGAL::ORIGIN+ortho,
                              arc.supporting_circle().squared_radius() + ortho.squared_length(),
                              sphere_.squared_radius());
    //warning : down to here
  }
  else{
    typename SK::Point_3 tmp_pt=typename SK::Point_3(0.,0.,0.) + (arc.supporting_circle().center()-sphere_.center());
    typename SK::FT gk=
      gamma_k(tmp_pt,
              arc.supporting_circle().squared_radius(),
              sphere_.squared_radius());
    
    tz=traits_.unsigned_tkz_coeff_normal(tmp_pt,gk);
    if ( is_arc_an_upper_one(arc,sphere_) )
      tz=-tz;
    //~ tz*=is_arc_an_upper_one(arc,sphere_)?(-1):(1);
    n=square_norm_tk_normal(tmp_pt,
                            arc.supporting_circle().squared_radius(),
                            sphere_.squared_radius());
  }
}


template<class SK,class Traits>
int 
Compare_to_right_of_arcs<SK,Traits>::sign_of_delta(
  const typename SK::Circular_arc_3& arc1, bool is_circle_1_threaded,
  const typename SK::Circular_arc_3& arc2, bool is_circle_2_threaded) const
{
  typename Traits::Tk_type tz1,tz2;
  typename SK::FT n1,n2;
  fill_tzk_n(arc1,tz1,n1,is_circle_1_threaded);
  fill_tzk_n(arc2,tz2,n2,is_circle_2_threaded);
  if (tz1==0 && tz2==0)
    return 0;
  if (tz1==0)
    return CGAL_NTS sign(tz2);
  if (tz2==0)
    return -CGAL_NTS sign(tz1);
  
  if (CGAL_NTS sign(tz1)!=CGAL_NTS sign(tz2))
    return CGAL_NTS sign(tz2);
  return CGAL_NTS sign(tz2)*CGAL_NTS sign(n1*CGAL::square(tz2)-n2*CGAL::square(tz1));
}


template<class SK,class Traits>
typename SK::FT 
Compare_to_right_of_arcs<SK,Traits>::give_rk(
  const typename SK::Circular_arc_3& arc) const
{
  return (is_arc_an_upper_one(arc,sphere_)?(1):(-1))*arc.supporting_circle().squared_radius();  
}

template<class SK,class Traits>
int
Compare_to_right_of_arcs<SK,Traits>::compare_for_delta_eq_0_threaded(
  const typename SK::Circular_arc_3& arc_threaded,
  const typename SK::Circular_arc_3& arc,
  bool is_supporting_circle_threaded) const
{
  if (!is_supporting_circle_threaded)
    return ( is_upper_arc<SK>(arc,sphere_) )?(-1):(1); //IN THAT CASE WE CAN OPTIMIZE AND USE THE Z-COORDINATES OF THE POINT AND CIRCLE THETA-EXTREMAL PT
  return (-CGAL_NTS sign(arc_threaded.supporting_circle().center().z()-arc.supporting_circle().center().z()));
}
  
template<class SK,class Traits>
int 
Compare_to_right_of_arcs<SK,Traits>::compare_for_delta_eq_0(
  const typename SK::Circular_arc_3& arc1,
  bool is_circle_1_threaded,
  const typename SK::Circular_arc_3& arc2,
  bool is_circle_2_threaded) const
{
  if (is_circle_1_threaded) 
    return  compare_for_delta_eq_0_threaded(arc1,arc2,is_circle_2_threaded);
  if (is_circle_2_threaded) 
    return -compare_for_delta_eq_0_threaded(arc2,arc1,is_circle_1_threaded);
  typename SK::FT rc1=give_rk(arc1);
  typename SK::FT rc2=give_rk(arc2);
  CGAL_precondition(rc1!=0 && rc2!=0);
  if(CGAL_NTS sign(rc1)*CGAL_NTS sign(rc2)<0)
    return (CGAL_NTS sign(rc1)>0)?(1):(-1);
  return (rc1<rc2)?(1):(-1);
}


template<class SK,class Traits>
CGAL::Comparison_result Compare_to_right_of_arcs<SK,Traits>::operator()(
                          const typename SK::Circular_arc_3& arc1,
                          const typename SK::Circular_arc_3& arc2,
                          bool do_it_to_left
                        ) const
{
  bool is_circle_1_threaded = (classify_circle_3<SK>(arc1.supporting_circle(),sphere_)==CGAL::THREADED);
  bool is_circle_2_threaded = (classify_circle_3<SK>(arc2.supporting_circle(),sphere_)==CGAL::THREADED);
  
  int res=(do_it_to_left?-1:1)*sign_of_delta(arc1,is_circle_1_threaded,arc2,is_circle_2_threaded);
  if (res==0)
    res=compare_for_delta_eq_0(arc1,is_circle_1_threaded,arc2,is_circle_2_threaded);
  CGAL_precondition(res!=0);
  if (res<0) return CGAL::LARGER;
  return CGAL::SMALLER;
}

}
}

#endif //CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_TO_RIGHT_H
