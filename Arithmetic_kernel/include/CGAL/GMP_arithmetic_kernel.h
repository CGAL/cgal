// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Number_types/include/CGAL/Arithmetic_kernel.h $
// $Id: Arithmetic_kernel.h 47259 2008-12-06 21:47:11Z afabri $
// 
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================
//
//    \brief provide class Arithmetic_kernel, a collection of number types. 
//

/*! \file CGAL/Arithmetic_kernel.h
 *  \brief Declarations pertaining to CGAL::Arithmetic_kernel
 */

#ifndef CGAL_GMP_ARITHMETIC_KERNEL_H
#define CGAL_GMP_ARITHMETIC_KERNEL_H

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel/Arithmetic_kernel_base.h>
#include <CGAL/Get_arithmetic_kernel.h>


#ifdef CGAL_USE_GMP 
#ifdef CGAL_USE_MPFI

#define CGAL_HAS_GMP_ARITHMETIC_KERNEL 

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpfi.h>

namespace CGAL {

/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The GMP set of exact number types
 */
class GMP_arithmetic_kernel : public internal::Arithmetic_kernel_base {
public:
  typedef CGAL::Gmpz           Integer;
  typedef CGAL::Gmpq           Rational;
  typedef CGAL::Gmpfr          Bigfloat;
  typedef CGAL::Gmpfi Bigfloat_interval;
};
    
template <>
struct Get_arithmetic_kernel<Gmpz> {
  typedef GMP_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<Gmpq>{
  typedef GMP_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<Gmpfr>{
  typedef GMP_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<Gmpfi>{
  typedef GMP_arithmetic_kernel Arithmetic_kernel;
};

} //namespace CGAL

#endif //CGAL_USE_MPFI
#endif //CGAL_USE_GMP  
#endif // CGAL_ARITHMETIC_KERNEL_H
// EOF
