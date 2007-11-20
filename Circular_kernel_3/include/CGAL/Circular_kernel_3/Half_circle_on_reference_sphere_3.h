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

#ifndef CGAL_SPHERICAL_HALF_CIRCLE_ON_REFERENCE_SPHERE_H
#define CGAL_SPHERICAL_HALF_CIRCLE_ON_REFERENCE_SPHERE_H

namespace CGAL {
  namespace CGALi {

  template <class SK>
  class Half_circle_on_reference_sphere_3
  {
    protected:
    const typename SK::Circle_on_reference_sphere_3& C;
    CGAL::Hcircle_type pos;
    public:
    const typename SK::Circle_on_reference_sphere_3& supporting_circle() const {return C;}
    const CGAL::Hcircle_type& get_position() const {return pos;}
    Half_circle_on_reference_sphere_3(const typename SK::Circle_on_reference_sphere_3& C,CGAL::Hcircle_type pos):C(C),pos(pos){}
  };
    
  }
}

#endif
