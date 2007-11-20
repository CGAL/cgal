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

#ifndef CGAL_HALF_CIRCLE_ON_REFERENCE_SPHERE_3_H
#define CGAL_HALF_CIRCLE_ON_REFERENCE_SPHERE_3_H

namespace CGAL {
  
template < typename SphericalKernel >
class Half_circle_on_reference_sphere_3
  : public SphericalKernel::Kernel_base::Half_circle_on_reference_sphere_3
{
  typedef typename SphericalKernel::Kernel_base::Half_circle_on_reference_sphere_3 
                                           RHalf_circle_on_reference_sphere_3;

public:
  typedef SphericalKernel   R; 
  typedef RHalf_circle_on_reference_sphere_3 Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  //~ Half_circle_on_reference_sphere_3()
  //~ :RHalf_circle_on_reference_sphere_3(){}
  
  Half_circle_on_reference_sphere_3(const RHalf_circle_on_reference_sphere_3& p)
  :RHalf_circle_on_reference_sphere_3(p){}

  Half_circle_on_reference_sphere_3(const typename SphericalKernel::Circle_on_reference_sphere_3& C,CGAL::Hcircle_type pos)    
  : RHalf_circle_on_reference_sphere_3(typename R::Construct_half_circle_on_reference_sphere_3()(C,pos))
  {}

  typename Qualified_result_of<typename R::Compute_half_circle_position_3,Half_circle_on_reference_sphere_3>::type
  get_position() const
  { return typename R::Compute_half_circle_position_3()(*this);}
    
  typename Qualified_result_of<typename R::Compute_supporting_circle_on_reference_sphere_3,Half_circle_on_reference_sphere_3>::type
  supporting_circle() const
  { return typename R::Compute_supporting_circle_on_reference_sphere_3()(*this);}
   
    
};
  
}
#endif

