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

#ifndef CGAL_CIRCULAR_KERNEL_FUNCTIONS_ON_LINE_2_H
#define CGAL_CIRCULAR_KERNEL_FUNCTIONS_ON_LINE_2_H

namespace CGAL {
namespace LinearFunctors {

  template < class CK >
  typename CK::Polynomial_1_2
  get_equation( const typename CK::Line_2 & L )
  {
    typedef typename CK::RT RT;
 
    return typename CK::Polynomial_1_2(L.a(),L.b(),L.c());
  }
  
  template < class CK >
  typename CK::Line_2  
  construct_line_2 ( const typename CK::Polynomial_1_2 &eq )
  {
    return typename CK::Line_2(eq[2],eq[1],eq[0]); 
  }
  
  template < class CK >
  bool
  has_on(const typename CK::Line_2 & l,
	 const typename CK::Circular_arc_point_2 &p)
  {
    typedef typename CK::Algebraic_kernel            AK;
    typedef typename CK::Polynomial_1_2 Polynomial_1_2;
    Polynomial_1_2 equation = CGAL::LinearFunctors::get_equation<CK>(l);

    return(AK().sign_at_object()(equation,p.coordinates())== ZERO);
  }

  template < class CK >
  inline bool
  non_oriented_equal(const typename CK::Line_2 & a1,
	             const typename CK::Line_2 & a2) {
    if(identical(a1,a2)) return true;
    const typename CK::RT &a1c = a1.a(); 
    const typename CK::RT &b1c = a1.b();
    const typename CK::RT &c1c = a1.c();
    const typename CK::RT &a2c = a2.a(); 
    const typename CK::RT &b2c = a2.b();
    const typename CK::RT &c2c = a2.c();
    return (a1c*b2c == a2c*b1c) &&
           (a1c*c2c == a2c*c1c) &&
           (b1c*c2c == b2c*c1c);
  }

  template< class CK, class OutputIterator>
  OutputIterator
  intersect_2( const typename CK::Line_2 & l,
	       const typename CK::Circle_2 & c,
	       OutputIterator res )
  {
    typedef typename CK::Algebraic_kernel            AK;
    typedef typename CK::Polynomial_1_2              Equation_line;
    typedef typename CK::Polynomial_for_circles_2_2  Equation_circle; 
    typedef typename CK::Root_for_circles_2_2        Root_for_circles_2_2;
    
    Equation_line e1 = CGAL::get_equation<CK>(l);
    Equation_circle e2 = CGAL::get_equation<CK>(c);
    
    typedef std::vector< std::pair < Root_for_circles_2_2, unsigned > > 
      solutions_container;
    solutions_container solutions;

    AK().solve_object()(e1, e2, std::back_inserter(solutions)); 
    // to be optimized

    typedef typename CK::Circular_arc_point_2 Circular_arc_point_2;

    for ( typename solutions_container::iterator it = solutions.begin(); 
	  it != solutions.end(); ++it )
      {
	*res++ = make_object
	  (std::make_pair(Circular_arc_point_2(it->first), it->second ));
      }

    return res;
  }

} // namespace LinearFunctors
} // namespace CGAL
#endif // CGAL_CIRCULAR_KERNEL_FUNCTIONS_ON_LINE_2_H
