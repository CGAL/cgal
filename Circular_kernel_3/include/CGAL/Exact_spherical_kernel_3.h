// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a
// STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_EXACT_SPHERICAL_3_KERNEL_H
#define CGAL_EXACT_SPHERICAL_3_KERNEL_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Spherical_kernel_3.h>


#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpq.h>
#else
#  include <CGAL/MP_Float.h>
#  include <CGAL/Quotient.h>
#endif

namespace CGAL {


#ifdef CGAL_USE_GMP
  typedef CGAL::Gmpq                                           NT1;
#else
  typedef CGAL::Quotient<CGAL::MP_Float>                       NT1;
#endif

  typedef CGAL::Cartesian<NT1>                                 Linear_k1;
  typedef CGAL::Algebraic_kernel_for_spheres_2_3<NT1>          Algebraic_k1;
  typedef CGAL::Spherical_kernel_3<Linear_k1,Algebraic_k1>     Exact_spherical_kernel_3;

} //namespace CGAL

#endif // CGAL_EXACT_SPHERICAL_3_KERNEL_H
