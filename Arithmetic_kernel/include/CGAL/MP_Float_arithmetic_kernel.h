// Copyright (c) 2010 Max-Planck-Institute Saarbruecken (Germany).
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
