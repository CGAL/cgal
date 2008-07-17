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

#ifndef CGAL_SPHERICAL_KERNEL_THETA_REP_H
#define CGAL_SPHERICAL_KERNEL_THETA_REP_H

namespace CGAL {
  namespace CGALi {

template <class SK>
class Theta_rep{
  typedef typename SK::Algebraic_kernel::Root_of_2 Root_of_2;
  typedef typename SK::FT FT;    
  typedef std::pair<HQ_NT,Root_of_2> Rep;
  
  //
  typename SK::template Handle<Rep>::type  base;
  public:
  Theta_rep(const HQ_NT& hq,const Root_of_2& r):base(std::pair<HQ_NT,Root_of_2>(hq,r)){}
  Theta_rep(){}
  const HQ_NT& hq() const {return get(base).first;}
  const Root_of_2& ftheta() const {return CGAL::get(base).second;}
};    

  }
}

#endif
