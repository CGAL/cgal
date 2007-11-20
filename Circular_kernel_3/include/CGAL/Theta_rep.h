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

#ifndef CGAL_THETA_REP_H
#define CGAL_THETA_REP_H

namespace CGAL {

template < typename SphericalKernel >
class Theta_rep
  : public SphericalKernel::Kernel_base::Theta_rep
{
  typedef typename SphericalKernel::Kernel_base::Theta_rep 
                                           RTheta_rep;

  typedef typename SphericalKernel::Root_of_2             Root_of_2;
  typedef typename SphericalKernel::FT                    FT;
  typedef typename SphericalKernel::Algebraic_kernel      AK;
  typedef typename SphericalKernel::Root_for_spheres_2_3  Root_for_spheres_2_3;  

public:
  typedef SphericalKernel   R; 
  typedef RTheta_rep Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  Theta_rep(const RTheta_rep& p)
  :RTheta_rep(p){}    
    
    
  Theta_rep()
  : RTheta_rep(
    typename R::Construct_theta_rep()())
  {}

  Theta_rep(const HQ_NT& hq,const Root_of_2& r)
  : RTheta_rep(
    typename R::Construct_theta_rep()(hq,r))
  {}

    
  typename Qualified_result_of<typename R::Compute_theta_ftheta_3,Theta_rep>::type
  ftheta() const
  { return typename R::Compute_theta_ftheta_3()(*this);}
  
  typename Qualified_result_of<typename R::Compute_theta_hq_3,Theta_rep>::type
  hq() const
  { return typename R::Compute_theta_hq_3()(*this);}
  
    
};

}
#endif

