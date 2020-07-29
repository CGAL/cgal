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


#ifndef CGAL_LEDA_ARITHMETIC_KERNEL_H
#define CGAL_LEDA_ARITHMETIC_KERNEL_H

#include <CGAL/basic.h>

#ifdef CGAL_USE_LEDA

#define CGAL_HAS_LEDA_ARITHMETIC_KERNEL

#include <CGAL/Arithmetic_kernel/Arithmetic_kernel_base.h>
#include <CGAL/Get_arithmetic_kernel.h>

#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_real.h>
#include <CGAL/leda_bigfloat_interval.h>


namespace CGAL {

/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The LEDA set of exact number types
 */
class LEDA_arithmetic_kernel : public internal::Arithmetic_kernel_base {
public:
  //! exact integers
  typedef leda_integer Integer;
  //! exact rationals, constructible from integers
  typedef leda_rational Rational;
  //! exact root expressions, constructible from integers and rationals
  typedef leda_real Field_with_sqrt;

  // undocumented
  typedef leda_bigfloat          Bigfloat;
  typedef leda_bigfloat_interval Bigfloat_interval;

};



template <>
struct Get_arithmetic_kernel<leda::integer> {
  typedef LEDA_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<leda::rational>{
  typedef LEDA_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<leda::real>{
  typedef LEDA_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<leda::bigfloat>{
  typedef LEDA_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<CGAL::leda_bigfloat_interval>{
  typedef LEDA_arithmetic_kernel Arithmetic_kernel;
};

} //namespace CGAL


#endif // CGAL_USE_LEDA

#endif //  CGAL_LEDA_ARITHMETIC_KERNEL_H
