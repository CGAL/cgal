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

#ifndef CGAL_CIRCLE_ON_REFERENCE_SPHERE_3_H
#define CGAL_CIRCLE_ON_REFERENCE_SPHERE_3_H

namespace CGAL {

template < typename SphericalKernel >
class Circle_on_reference_sphere_3
  : public SphericalKernel::Kernel_base::Circle_on_reference_sphere_3
{
  typedef typename SphericalKernel::Kernel_base::Circle_on_reference_sphere_3 
                                           RCircle_on_reference_sphere_3;

  typedef typename SphericalKernel::Point_3               Point_3;
  typedef typename SphericalKernel::FT                    FT;
  typedef typename SphericalKernel::Algebraic_kernel      AK;
  typedef typename SphericalKernel::Root_for_spheres_2_3  Root_for_spheres_2_3;  

public:
  typedef SphericalKernel   R; 
  typedef RCircle_on_reference_sphere_3 Repd;

  const Repd& rep() const
  {
    return *this;
  }

  Repd& rep()
  {
    return *this;
  }

  Circle_on_reference_sphere_3()
  :RCircle_on_reference_sphere_3(){}
  
  Circle_on_reference_sphere_3(const RCircle_on_reference_sphere_3& p)
  :RCircle_on_reference_sphere_3(p){}

  Circle_on_reference_sphere_3(const FT& _r,const Point_3& _c,const typename SphericalKernel::Sphere_with_radius_3& S)
  : RCircle_on_reference_sphere_3(
    typename R::Construct_circle_on_reference_sphere_3()(_r,_c,S))
  {}

  Circle_on_reference_sphere_3(const FT& _r,const Point_3& _c,CGAL::Circle_type nat,const typename SphericalKernel::Sphere_with_radius_3& S)
  : RCircle_on_reference_sphere_3(
    typename R::Construct_circle_on_reference_sphere_3()(_r,_c,nat,S))
  {}    


  typename Qualified_result_of<typename R::Compute_type_of_circle_on_reference_sphere_3,Circle_on_reference_sphere_3>::type
  type_of_circle_on_reference_sphere() const
  { return typename R::Compute_type_of_circle_on_reference_sphere_3()(*this);}
    
  typename Qualified_result_of<typename R::Compute_supporting_sphere_radius_3,Circle_on_reference_sphere_3>::type
  supporting_sphere_radius() const
  { return typename R::Compute_supporting_sphere_radius_3()(*this);}
  
  typename Qualified_result_of<typename R::Compute_supporting_sphere_squared_radius_3,Circle_on_reference_sphere_3>::type
  supporting_sphere_squared_radius() const
  { return typename R::Compute_supporting_sphere_squared_radius_3()(*this);}
    
  typename Qualified_result_of<typename R::Compute_supporting_sphere_center_3,Circle_on_reference_sphere_3>::type
  supporting_sphere_center() const
  { return typename R::Compute_supporting_sphere_center_3()(*this);}
  
  typename Qualified_result_of<typename R::Compute_reference_sphere_3,Circle_on_reference_sphere_3>::type
  reference_sphere() const
  { return typename R::Compute_reference_sphere_3()(*this);}
  
  typename Qualified_result_of<typename R::Compute_extremal_point_z,Circle_on_reference_sphere_3>::type
  extremal_point_z() const
  { return typename R::Compute_extremal_point_z()(*this);} 
  
  typename Qualified_result_of<typename R::Compute_circle_center_coefficient_3,Circle_on_reference_sphere_3>::type
  circle_center_coefficient() const
  { return typename R::Compute_circle_center_coefficient_3()(*this);}  
  
   
    
};

}
#endif

