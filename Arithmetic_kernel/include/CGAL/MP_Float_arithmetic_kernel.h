// Copyright (c) 2010 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot
//
// ============================================================================
//
//    \brief provide class Arithmetic_kernel, a collection of number types.
//

#ifndef CGAL_MP_FLOAT_ARITHMETIC_KERNEL_H
#define CGAL_MP_FLOAT_ARITHMETIC_KERNEL_H

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel/Arithmetic_kernel_base.h>
#include <CGAL/Get_arithmetic_kernel.h>

#define CGAL_HAS_MP_FLOAT_ARITHMETIC_KERNEL

#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

namespace CGAL {

/*! \ingroup CGAL_Arithmetic_kernel
 *  \brief  The MP_Float set of exact number types
 */
class MP_Float_arithmetic_kernel : public internal::Arithmetic_kernel_base {
public:
  typedef MP_Float                       Integer;
  typedef CGAL::Quotient<MP_Float>       Rational;
  typedef MP_Float                       Bigfloat;
  struct  Not_implemented{}              Bigfloat_interval;
};

template <>
struct Get_arithmetic_kernel<MP_Float> {
  typedef MP_Float_arithmetic_kernel Arithmetic_kernel;
};
template <>
struct Get_arithmetic_kernel<Quotient<MP_Float> > {
  typedef MP_Float_arithmetic_kernel Arithmetic_kernel;
};

} //namespace CGAL

#endif // CGAL_ARITHMETIC_KERNEL_H
// EOF
