// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_EXACT_CIRCULAR_2_KERNEL_H
#define CGAL_EXACT_CIRCULAR_2_KERNEL_H

#include <CGAL/license/Circular_kernel_2.h>


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


#include <CGAL/Filtered_bbox_circular_kernel_2.h>

namespace CGAL {

namespace internal {

#ifdef CGAL_USE_GMP
  typedef CGAL::Gmpq                                     NT;
#else

  typedef Quotient<MP_Float>                             NT;

#endif

  typedef Cartesian<NT>                                  Linear_k;
  typedef Algebraic_kernel_for_circles_2_2<NT>           Algebraic_k;
  typedef Circular_kernel_2<Linear_k, Algebraic_k>       CK;

} // namespace internal

typedef Filtered_bbox_circular_kernel_2<internal::CK>    Exact_circular_kernel_2;

} //namespace CGAL

#endif // CGAL_EXACT_CIRCULAR_2_KERNEL_H
