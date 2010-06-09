// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Number_types/include/CGAL/Sqrt_extension/Scalar_factor_traits.h $
// $Id: Scalar_factor_traits.h 52628 2009-10-20 08:59:26Z lrineau $
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_GET_ARITHMETIC_KERNEL_H
#define CGAL_SQRT_EXTENSION_GET_ARITHMETIC_KERNEL_H

#include <CGAL/basic.h>
#include <CGAL/Get_arithmetic_kernel.h>

namespace CGAL {

template <class COEFF, class ROOT>
struct Get_arithmetic_kernel<Sqrt_extension<COEFF,ROOT> >{
  typedef Get_arithmetic_kernel<COEFF> GET;
  typedef typename GET::Arithmetic_kernel Arithmetic_kernel;
};

} //namespace CGAL

#endif
