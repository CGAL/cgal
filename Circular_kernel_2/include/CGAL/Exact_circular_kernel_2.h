// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_EXACT_CIRCULAR_2_KERNEL_H
#define CGAL_EXACT_CIRCULAR_2_KERNEL_H

#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>


//TODO: CORRECT THE MAKE_ROOT_OF_2 of GMPq GMPz
#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpq.h>
#else

#  include <CGAL/MP_Float.h>
#  include <CGAL/Quotient.h>

#endif


// maybe it is better to change to the bbox filtered one
//#include <CGAL/Lazy_circular_kernel_2.h>
#include <CGAL/Filtered_bbox_circular_kernel_2.h>

namespace CGAL {

namespace internal {

#ifdef CGAL_USE_GMP
  typedef CGAL::Gmpq                                           NT1;
#else

  typedef Quotient<MP_Float>                       NT1;

#endif

  typedef Cartesian<NT1>                                 Linear_k1;
  typedef Algebraic_kernel_for_circles_2_2<NT1>          Algebraic_k1;
  typedef Circular_kernel_2<Linear_k1, Algebraic_k1>     CK1;

//   typedef CGAL::Interval_nt_advanced                           NT2;
//   typedef CGAL::Cartesian<NT2>                                 Linear_k2;
//   typedef CGAL::Algebraic_kernel_for_circles_2_2<NT2>          Algebraic_k2;
//   typedef CGAL::Circular_kernel_2<Linear_k2,Algebraic_k2>      CK2;

//  typedef CGAL::Lazy_circular_kernel_2<CK1,CK2>
//  Exact_circular_kernel_2;

} // namespace internal

typedef Filtered_bbox_circular_kernel_2<internal::CK1>   Exact_circular_kernel_2;

} //namespace CGAL

#endif // CGAL_EXACT_CIRCULAR_2_KERNEL_H
