// Copyright (c) 2016 Inria.
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
// Author: Marc Glisse <marc.glisse@inria.fr>

#ifndef CGAL_GMPXX_ARITHMETIC_KERNEL_H
#define CGAL_GMPXX_ARITHMETIC_KERNEL_H

#include <CGAL/config.h>

#ifdef CGAL_USE_GMPXX

#include <CGAL/Arithmetic_kernel/Arithmetic_kernel_base.h>
#include <CGAL/Get_arithmetic_kernel.h>

#include <CGAL/gmpxx.h>

namespace CGAL {
/** \ingroup CGAL_Arithmetic_kernel
 *  \brief The GMPXX set of exact number types
 */
struct GMPXX_arithmetic_kernel : internal::Arithmetic_kernel_base {
  typedef mpz_class Integer;
  typedef mpq_class Rational;
};

template <class T, class U>
struct Get_arithmetic_kernel<__gmp_expr<T, U> > {
  typedef GMPXX_arithmetic_kernel Arithmetic_kernel;
};
} //namespace CGAL
#endif //CGAL_USE_GMPXX
#endif
