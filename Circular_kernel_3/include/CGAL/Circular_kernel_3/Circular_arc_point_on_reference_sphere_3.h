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
// $URL: svn+ssh://sloriot@scm.gforge.inria.fr/svn/cgal/trunk/Circular_kernel_3/include/CGAL/Circular_kernel_3/Circular_arc_3.h $
// $Id: Circular_arc_3.h 40627 2007-10-16 15:00:59Z sloriot $
//
// Author(s) : Loriot Sebastien <Sebastien.Loriot@sophia.inria.fr>
//                   Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                   Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//                   Pedro Machado    <tashimir@gmail.com>
//                   Julien Hazebrouck
//                   Damien Leroy

#ifndef CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H
#define CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H

#include <iostream>

#include <CGAL/Circular_arc_point_3.h>
#include <CGAL/Circular_kernel_3/constant.h>

namespace CGAL {
  namespace CGALi {
  
    
  template<class SK>
  class Circular_arc_point_on_reference_sphere_3:public Circular_arc_point_3<SK>{
    typedef typename SK::Algebraic_kernel::Root_of_2 Root_of_2;
    typedef typename SK::FT FT;
    typedef typename SK::Algebraic_kernel AK;
    typedef Circular_arc_point_3<SK> T_Circular_arc_point_3;
    typedef typename SK::Theta_rep Theta_rep;
    //---------------
    Theta_rep Trep;
    public:
    Circular_arc_point_on_reference_sphere_3(const FT& ftheta,const FT& xt,const FT& yt,const FT& zt,const CGAL::HQ_NT& _hq)
      :T_Circular_arc_point_3(typename SK::Point_3(xt,yt,zt)),Trep(_hq,ftheta){};//critical point of non normal circles
        
    Circular_arc_point_on_reference_sphere_3(const CGAL::HQ_NT& _hq,const Root_of_2& ftheta,const Root_of_2& x_,const Root_of_2& y_,const Root_of_2& z_)
      :T_Circular_arc_point_3(x_,y_,z_),Trep(_hq,ftheta){};

    Circular_arc_point_on_reference_sphere_3(const CGAL::HQ_NT& _hq,const Root_of_2& ftheta,const typename SK::Algebraic_kernel::Root_for_spheres_2_3& rfs)
      :T_Circular_arc_point_3(rfs),Trep(_hq,ftheta){};
        
    Circular_arc_point_on_reference_sphere_3():T_Circular_arc_point_3(FT(0),FT(0),FT(0)),Trep(-1,FT(0)){};
      
    Circular_arc_point_on_reference_sphere_3(const CGAL::HQ_NT& hq,const typename AK::Root_for_spheres_2_3& R):T_Circular_arc_point_3(R),Trep(hq,CGAL::auto_ftype(hq)==CGAL::TAN?(R.y()/R.x()):(R.x()/R.y())){};            
    
    Circular_arc_point_on_reference_sphere_3(const CGAL::HQ_NT& hq,const T_Circular_arc_point_3& R):T_Circular_arc_point_3(R),Trep(hq,CGAL::auto_ftype(hq)==CGAL::TAN?(R.y()/R.x()):(R.x()/R.y())){};            
    //~ static inline Circular_arc_point_on_reference_sphere_3 VirtualPt_to_point_on_sphere(){
      //~ return Circular_arc_point_on_reference_sphere_3(FT(0),FT(0),FT(0),FT(0),HQ_NT(-1));
    //~ };

    //~ static inline Circular_arc_point_on_reference_sphere_3 Root_for_sphere_to_point_on_sphere(const CGAL::HQ_NT& hq,const typename AK::Root_for_spheres_2_3& R){
      //~ return Circular_arc_point_on_reference_sphere_3(hq,auto_ftype(hq)==TAN?(R.y()/R.x()):(R.x()/R.y()),R);
    //~ };      
      
      
    const Theta_rep& theta_rep() const {return Trep;};
    const Root_of_2& get_f_of_theta() const {return theta_rep().ftheta();};  
    const CGAL::HQ_NT& get_hq() const {return theta_rep().hq();}

    double get_theta_approx() const{
      double ax=CGAL::to_double(this->x());
      double ay=CGAL::to_double(this->y());
      return ( (atan2 (ay,ax)<0)?(atan2 (ay,ax)+2.*M_PI):(atan2 (ay,ax)) );
    };
    
    CGAL::Cartesian<double>::Point_3 get_point_approx() const {//just for intersection and critical points of normal circles
      return CGAL::Cartesian<double>::Point_3(CGAL::to_double(this->x()),CGAL::to_double(this->y()),CGAL::to_double(this->z()));
    }
    
    
    
  };    
    
  } // namespace CGALi
} // namespace CGAL

#endif //CGAL_SPHERICAL_KERNEL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H

