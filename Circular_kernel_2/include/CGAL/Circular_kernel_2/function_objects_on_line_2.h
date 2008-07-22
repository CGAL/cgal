// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_FUNCTION_OBJECTS_ON_LINE_2_H
#define CGAL_CIRCULAR_KERNEL_FUNCTION_OBJECTS_ON_LINE_2_H

#include <CGAL/Circular_kernel_2/internal_functions_on_line_2.h>

namespace CGAL {
namespace LinearFunctors {

  template < class CK >
  class Construct_line_2 : public  CK::Linear_kernel::Construct_line_2
  {
    public:

    typedef typename CK::Line_2 result_type;

    result_type
    operator() ( const typename CK::Polynomial_1_2 &eq )
      {
	return construct_line_2<CK>(eq);
      }
  };

  template < class CK >
  class Get_equation
  {
    public:

    typedef typename CK::Polynomial_1_2 result_type;

    result_type
    operator() ( const typename CK::Line_2 & l )
      {
	return LinearFunctors::get_equation<CK>(l);
      }
  };

} // namespace LinearFunctors
} // namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_FUNCTION_OBJECTS_ON_LINE_2_H
