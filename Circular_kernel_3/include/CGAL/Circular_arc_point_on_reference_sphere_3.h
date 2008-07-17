// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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

#ifndef CGAL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H
#define CGAL_CIRCULAR_ARC_POINT_ON_REFERENCE_SPHERE_3_H

namespace CGAL {

template < typename SphericalKernel >
class Circular_arc_point_on_reference_sphere_3
  : public SphericalKernel::Kernel_base::Circular_arc_point_on_reference_sphere_3
{
  typedef typename SphericalKernel::Kernel_base::Circular_arc_point_on_reference_sphere_3 
                                           RCircular_arc_point_on_reference_sphere_3;

  typedef typename SphericalKernel::Root_of_2             Root_of_2;
  typedef typename SphericalKernel::Point_3               Point_3;
  typedef typename SphericalKernel::Plane_3               Plane_3;
  typedef typename SphericalKernel::Line_3                Line_3;
  typedef typename SphericalKernel::Circle_3              Circle_3;
  typedef typename SphericalKernel::Sphere_3              Sphere_3;
  typedef typename SphericalKernel::FT                    FT;
  typedef typename SphericalKernel::Algebraic_kernel      AK;


public:
  typedef typename SphericalKernel::Root_for_spheres_2_3 
    Root_for_spheres_2_3;
  typedef SphericalKernel   R; 
  typedef RCircular_arc_point_on_reference_sphere_3 Repd;

  const Repd& rep() const
  {
    return *this;
  }

  Repd& rep()
  {
    return *this;
  }

  Circular_arc_point_on_reference_sphere_3(const RCircular_arc_point_on_reference_sphere_3& p)
  :RCircular_arc_point_on_reference_sphere_3(p){}    
  
  Circular_arc_point_on_reference_sphere_3()
  : RCircular_arc_point_on_reference_sphere_3(
    typename R::Construct_circular_arc_point_on_reference_sphere_3()())
  {}

  Circular_arc_point_on_reference_sphere_3(const FT& ftheta,const FT& xt,const FT& yt,const FT& zt,const HQ_NT& _hq)
  : RCircular_arc_point_on_reference_sphere_3(
    typename R::Construct_circular_arc_point_on_reference_sphere_3()(ftheta,xt,yt,zt,_hq))
  {}
    
  Circular_arc_point_on_reference_sphere_3(const HQ_NT& _hq,const Root_of_2& ftheta,const Root_of_2& x_,const Root_of_2& y_,const Root_of_2& z_)
  : RCircular_arc_point_on_reference_sphere_3(
    typename R::Construct_circular_arc_point_on_reference_sphere_3()(_hq,ftheta,x_,y_,z_))
  {}
    
  Circular_arc_point_on_reference_sphere_3(const HQ_NT& hq,const typename AK::Root_for_spheres_2_3& R)
  : RCircular_arc_point_on_reference_sphere_3(hq,R)
  {}
    
  Circular_arc_point_on_reference_sphere_3(const HQ_NT& hq,const typename SphericalKernel::Circular_arc_point_3& R)
  : RCircular_arc_point_on_reference_sphere_3(hq,R)
  {}    

  Circular_arc_point_on_reference_sphere_3(const HQ_NT& _hq,const Root_of_2& ftheta,const typename AK::Root_for_spheres_2_3& rfs)
  : RCircular_arc_point_on_reference_sphere_3(
    typename R::Construct_circular_arc_point_on_reference_sphere_3()(_hq,ftheta,rfs))
  {}
    

  typename Qualified_result_of<typename R::Compute_circular_theta_rep_3,Circular_arc_point_on_reference_sphere_3>::type
  theta_rep() const
  { return typename R::Compute_circular_theta_rep_3()(*this);}
  
  
      
  typename Qualified_result_of<typename R::Compute_circular_theta_3,Circular_arc_point_on_reference_sphere_3>::type
  //const Root_of_2 &
  get_f_of_theta() const
  { return typename R::Compute_circular_theta_3()(*this);}  

  typename Qualified_result_of<typename R::Compute_circular_hq_3,Circular_arc_point_on_reference_sphere_3>::type
  get_hq() const
  { return typename R::Compute_circular_hq_3()(*this);}    
  
  //~ Bbox_3 bbox() const
  //~ { return typename R::Construct_bbox_3()(*this); }

};


}

#endif
