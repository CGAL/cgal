// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================
//
//    \brief provide class Arithmetic_kernel, a collection of number types. 
//

#ifndef CGAL_GMP_ARITHMETIC_KERNEL_H
#define CGAL_GMP_ARITHMETIC_KERNEL_H

#include <CGAL/config.h>

#ifdef CGAL_USE_GMP 

#include <CGAL/Arithmetic_kernel/Arithmetic_kernel_base.h>
#include <CGAL/Get_arithmetic_kernel.h>



#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>

#ifdef CGAL_USE_MPFI
#define CGAL_HAS_GMP_ARITHMETIC_KERNEL 
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpfi.h>
#endif //CGAL_USE_MPFI

namespace CGAL {

/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The GMP set of exact number types
 */
class GMP_arithmetic_kernel : public internal::Arithmetic_kernel_base {
public:
  typedef CGAL::Gmpz           Integer;
  typedef CGAL::Gmpq           Rational;
  #ifdef CGAL_USE_MPFI
  typedef CGAL::Gmpfr          Bigfloat;
  typedef CGAL::Gmpfi Bigfloat_interval;
  #endif //CGAL_USE_MPFI
};
    
template <>
struct Get_arithmetic_kernel<Gmpz> {
  typedef GMP_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<Gmpq>{
  typedef GMP_arithmetic_kernel Arithmetic_kernel;
};

#ifdef CGAL_USE_MPFI
template <>
struct Get_arithmetic_kernel<Gmpfr>{
  typedef GMP_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<Gmpfi>{
  typedef GMP_arithmetic_kernel Arithmetic_kernel;
};
#endif //CGAL_USE_MPFI

} //namespace CGAL

#endif //CGAL_USE_GMP

#endif // CGAL_ARITHMETIC_KERNEL_H
// EOF
