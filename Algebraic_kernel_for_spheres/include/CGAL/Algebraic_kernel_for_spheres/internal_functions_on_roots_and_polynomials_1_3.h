// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado

#ifndef CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_1_3_H
#define CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_1_3_H

#include <CGAL/license/Circular_kernel_3.h>


namespace CGAL {
  namespace AlgebraicSphereFunctors {

  template < class AK, class OutputIterator >
  inline 
  OutputIterator
    solve( const typename AK::Polynomial_1_3 & e1,
	   const typename AK::Polynomial_1_3 & e2,
	   const typename AK::Polynomial_1_3 & e3,
	 OutputIterator res )
  {
    typedef typename AK::FT FT;
    typedef typename AK::Root_for_spheres_2_3 Root_for_spheres_2_3;
    CGAL_kernel_precondition(!(same_solutions<FT>(e1,e2) || same_solutions<FT>(e1,e3) ||
                               same_solutions<FT>(e2,e3)));
    const FT &a1 = e1.a();
    const FT &a2 = e2.a();
    const FT &a3 = e3.a();
    const FT &b1 = e1.b();
    const FT &b2 = e2.b();
    const FT &b3 = e3.b();
    const FT &c1 = e1.c();
    const FT &c2 = e2.c();
    const FT &c3 = e3.c();
    const FT &d1 = e1.d();
    const FT &d2 = e2.d();
    const FT &d3 = e3.d();
    FT denominateur = (a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1-a3*b2*c1+a3*b1*c2);
    //if denominateur == 0 it's because the planes are parallel
    if (denominateur == 0) return res;
    FT z = -(a2*b3*d1-a1*b3*d2+a1*b2*d3-a3*b2*d1-a2*b1*d3+a3*b1*d2)/denominateur;
    FT y = (-a1*d2*c3+a1*c2*d3-a2*c1*d3-a3*c2*d1+a3*c1*d2+a2*c3*d1)/denominateur;
    FT x = -(-b1*d2*c3+b1*c2*d3+b3*c1*d2-b2*c1*d3+b2*d1*c3-b3*c2*d1)/denominateur;
    *res++ = std::make_pair(Root_for_spheres_2_3(x,y,z), 
                            static_cast<unsigned>(1));
    return res;
  }

  }
}

#endif //CGAL_ALGEBRAIC_KERNEL_FUNCTIONS_ON_ROOTS_AND_POLYNOMIALS_1_3_H
