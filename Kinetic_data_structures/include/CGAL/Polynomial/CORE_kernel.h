// Copyright (c) 2005  Stanford University (USA).
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
// $Id$ $Date$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_POLYNOMIAL_CORE_KERNEL_H
#define CGAL_POLYNOMIAL_POLYNOMIAL_CORE_KERNEL_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Polynomial/Kernel.h>
#include <CGAL/Polynomial/internal/CORE_polynomial.h>
#include <CGAL/Polynomial/CORE_Expr_root_stack.h>

namespace CGAL { namespace POLYNOMIAL {



// CORE_Expr_root_stack::FT
class CORE_kernel: public Kernel<internal::CORE_polynomial, CORE_Expr_root_stack,
				 CORE::Expr >
{
  typedef CORE_kernel This;
  typedef Kernel<internal::CORE_polynomial, CORE_Expr_root_stack, CORE::Expr> P;
public:

  CORE_kernel(const CORE_Expr_root_stack::Traits &tr=CORE_Expr_root_stack::Traits()):
    P(tr){}

};

} } //namespace CGAL::POLYNOMIAL
#endif
