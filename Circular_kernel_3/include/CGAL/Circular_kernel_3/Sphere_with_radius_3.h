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


#ifndef CGAL_SPHERICAL_KERNEL_SPHERE_WITH_RADIUS_3_H
#define CGAL_SPHERICAL_KERNEL_SPHERE_WITH_RADIUS_3_H

//~ #include <CGAL/Cartesian/Sphere_3.h>
#include <CGAL/Sphere_3.h>

namespace CGAL {
  namespace CGALi {
    
  template<class SK>
  class Sphere_with_radius_3:public CGAL::Sphere_3<SK>{
    typedef typename SK::Algebraic_kernel::Root_of_2 Root_of_2;
    typedef typename SK::FT FT;
    typedef typename SK::Algebraic_kernel AK;
    typedef typename SK::Sphere_3 Sphere_3;
    //---------------
    typename SK::template Handle<FT>::type  hrad;
    public:
    typedef typename SK::Point_3 Point_3;
    Sphere_with_radius_3():Sphere_3(){};
    Sphere_with_radius_3(const FT& _r,const Point_3& _c):Sphere_3(_c,_r*_r),hrad(_r){};
    Sphere_with_radius_3(const Sphere_3& S):Sphere_3(S),hrad(-1){}
    const FT& radius() const {CGAL_precondition(CGAL::get(hrad)!=-1); return CGAL::get(hrad);}
  };
    
  } // namespace CGALi
} // namespace CGAL

#endif //CGAL_SPHERICAL_KERNEL_SPHERE_WITH_RADIUS_3_H
