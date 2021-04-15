// Copyright (c) 2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================
//
//    \brief provide class LEDA_arithmetic_kernel, a collection of number types.
//



#ifndef CGAL_CORE_ARITHMETIC_KERNEL_H
#define CGAL_CORE_ARITHMETIC_KERNEL_H

#include <CGAL/basic.h>

#ifdef CGAL_USE_CORE

#define CGAL_HAS_CORE_ARITHMETIC_KERNEL

#include <CGAL/Arithmetic_kernel/Arithmetic_kernel_base.h>
#include <CGAL/Get_arithmetic_kernel.h>

#include <CGAL/CORE_BigFloat.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_Expr.h>

namespace CGAL {

/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The CORE set of exact number types
 */
class CORE_arithmetic_kernel : public internal::Arithmetic_kernel_base {
public:
    //! exact integers
    typedef CORE::BigInt Integer;
    //! exact float nummber
    typedef CORE::BigRat Exact_float_number;
    //! exact rationals, constructible from integers
    typedef CORE::BigRat Rational;
    //! exact root expressions, constructible from integers and rationals
    typedef CORE::Expr Field_with_sqrt;
    // undocumented
    //typedef CORE::BigFloat          Bigfloat;
    typedef CORE::BigFloat          Bigfloat_interval;

};


template <>
struct Get_arithmetic_kernel<CORE::BigInt>{
  typedef CORE_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<CORE::BigRat>{
  typedef CORE_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<CORE::Expr>{
  typedef CORE_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<CORE::BigFloat>{
  typedef CORE_arithmetic_kernel Arithmetic_kernel;
};

} //namespace CGAL

#endif // CGAL_USE_CORE

#endif //  CGAL_CORE_ARITHMETIC_KERNEL_H
