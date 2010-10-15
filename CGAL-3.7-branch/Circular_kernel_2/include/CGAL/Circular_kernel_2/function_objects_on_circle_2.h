// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)


#ifndef CGAL_CIRCULAR_KERNEL_FUNCTION_OBJECTS_ON_CIRCLE_2_H
#define CGAL_CIRCULAR_KERNEL_FUNCTION_OBJECTS_ON_CIRCLE_2_H

#include <CGAL/Circular_kernel_2/internal_functions_on_circle_2.h>

#include <CGAL/Circular_kernel_2/function_objects_on_line_2.h> 
// to be removed when CGAL::Kernel has a Get_equation

namespace CGAL {
namespace CircularFunctors {

  template < class CK >
  class Construct_circle_2 : public  CK::Linear_kernel::Construct_circle_2
  {
  public:
    
    typedef typename CK::Circular_arc_2 Circular_arc_2;
    typedef typename CK::Linear_kernel::Construct_circle_2::result_type
      result_type;

    using CK::Linear_kernel::Construct_circle_2::operator();

    result_type
    operator() ( const typename CK::Polynomial_for_circles_2_2 &eq )
      {
	      return construct_circle_2<CK>(eq);
      }

	  const result_type& 
	  operator() (const Circular_arc_2 & a) const
	    {
	      return (a.rep().supporting_circle());
	    }

  };

} // namespace CircularFunctors

#ifndef CGAL_CFG_DONT_OVERLOAD_TOO_MUCH
  template < typename K>
  struct Qualified_result_of<CircularFunctors::Construct_circle_2<K>,
                           typename K::Circular_arc_2>
  {
    typedef const typename K::Circle_2 &   type;
  };
#endif

} // namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_FUNCTION_OBJECTS_ON_CIRCLE_2_H
