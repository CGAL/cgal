// Copyright (c) 2008-2009 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL$
// $Id$
// 
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================
//
//    \brief provide class Arithmetic_kernel, a collection of number types. 
//

// 

/*! \file CGAL/Arithmetic_kernel.h
 *  \brief Declarations pertaining to CGAL::Arithmetic_kernel
 */

#ifndef CGAL_ARITHMETIC_KERNEL_H
#define CGAL_ARITHMETIC_KERNEL_H

#include <CGAL/basic.h>
#include <CGAL/CORE_arithmetic_kernel.h>
#include <CGAL/LEDA_arithmetic_kernel.h>
#include <CGAL/GMP_arithmetic_kernel.h>


// Define a default Arithmetic_kernel GMP, CORE, LEDA 

namespace CGAL{

#ifndef CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
#if defined(CGAL_HAS_LEDA_ARITHMETIC_KERNEL) 
typedef LEDA_arithmetic_kernel Arithmetic_kernel;
#define CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL 1
#endif // CGAL_USE_LEDA
#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

#ifndef CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
#if defined(CGAL_HAS_GMP_ARITHMETIC_KERNEL) 
typedef GMP_arithmetic_kernel Arithmetic_kernel;
#define CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL 1
#endif // CGAL_USE_GMP
#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

#ifndef CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
#if defined(CGAL_HAS_CORE_ARITHMETIC_KERNEL) 
typedef CORE_arithmetic_kernel Arithmetic_kernel;
#define CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL 1
#endif // CGAL_USE_CORE
#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

} // namespace CGAL 


// Macro to snap typedefs in Arithmetic_kernel
#define CGAL_SNAP_ARITHMETIC_KERNEL_TYPEDEFS(AT) \
  typedef typename AT::Integer Integer; \
  typedef typename AT::Rational Rational; \
  typedef typename AT::Field_with_sqrt Field_with_sqrt; 

#endif // CGAL_ARITHMETIC_KERNEL_H
// EOF
