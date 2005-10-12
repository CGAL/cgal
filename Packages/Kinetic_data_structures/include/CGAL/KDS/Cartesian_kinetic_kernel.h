// Copyright (c) 2005  Stanford University (USA).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_KINETIC_KERNEL_H
#define CGAL_KDS_KINETIC_KERNEL_H
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/internal/Kernel/Cartesian_kinetic_kernel_base.h>



CGAL_KDS_BEGIN_NAMESPACE

//! A kinetic kernel using cartesian coordinates
/*!  It takes a PolynomialKernel as a template parameter. The
  PolynomialKernel is used to define the Motion_function and the
  Certificate_function.
*/
template <class Polynomial_k>
class Cartesian_kinetic_kernel: 
  public internal::Cartesian_kinetic_kernel_base<Polynomial_k, 
						 Cartesian_kinetic_kernel<Polynomial_k> > {
  typedef internal::Cartesian_kinetic_kernel_base<Polynomial_k, 
						  Cartesian_kinetic_kernel<Polynomial_k> > P;
public:
  Cartesian_kinetic_kernel(Polynomial_k pk): P(pk){}
  Cartesian_kinetic_kernel(){}
};

CGAL_KDS_END_NAMESPACE

//#include <CGAL/KDS_internals/kernel_undefs.h>
#endif
