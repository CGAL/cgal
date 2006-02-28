// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_2_H
#define CGAL_ALGEBRAIC_KERNEL_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_2_H

#include <CGAL/Algebraic_kernel/internal_functions_on_roots_and_polynomials_2_2.h>
#include <CGAL/Algebraic_kernel/internal_functions_on_roots_and_polynomial_1_2_and_2_2.h>
#include <CGAL/Algebraic_kernel/internal_functions_comparison_root_for_circles_2_2.h>

CGAL_BEGIN_NAMESPACE

namespace AlgebraicFunctors {
  
  template < class AK >
  class Solve
  {
    typedef typename AK::Polynomial_for_circles_2_2 Equation_Circle;
    typedef typename AK::Polynomial_1_2             Equation_Line;

  public:
    
    template < class OutputIterator >
    OutputIterator
    operator()(const Equation_Circle & e1, 
	       const Equation_Circle & e2, 
	       OutputIterator res) const
    { return AlgebraicFunctors::solve<AK> ( e1, e2, res); }
    
    template < class OutputIterator >
    OutputIterator
    operator()(const Equation_Line & e1, 
	       const Equation_Circle & e2, 
	       OutputIterator res) const
    { return AlgebraicFunctors::solve<AK> ( e1, e2, res); }

    template < class OutputIterator >
    OutputIterator
    operator()(const Equation_Line & e1, 
	       const Equation_Line & e2, 
	       OutputIterator res) const
    { return AlgebraicFunctors::solve<AK> ( e1, e2, res); }

  };

  template < class AK >
  class Construct_polynomial_for_circles_2_2
  {
    typedef typename AK::RT                         RT;
    typedef typename AK::Polynomial_for_circles_2_2 Polynomial_for_circles_2_2;
    
  public:
    
    Polynomial_for_circles_2_2
    operator()(const RT& xc, const RT& yc, const RT& r_sq) const
    { return Polynomial_for_circles_2_2(xc, yc, r_sq); }

  };

  template < class AK >
  class Construct_polynomial_1_2
  {
    typedef typename AK::RT             RT;
    typedef typename AK::Polynomial_1_2 Polynomial_1_2;
    
  public:

    Polynomial_1_2
    operator()( const RT& a, const RT& b, const RT& c) const
    { return Polynomial_1_2(a, b, c); }

  };

  template < class AK >
  class Sign_at
  {    
  public:
    
    Sign 
    operator()( const typename AK::Polynomial_for_circles_2_2 & equation,
		const typename AK::Root_for_circles_2_2 & r ) const
    { return AlgebraicFunctors::sign_at<AK>(equation, r); }

    Sign 
    operator()( const typename AK::Polynomial_1_2 & equation,
		const typename AK::Root_for_circles_2_2 & r ) const
    { return AlgebraicFunctors::sign_at<AK>(equation, r); }

  };
    
  template < class AK >
  class X_critical_points
  {    
  public:

    typename AK::Root_for_circles_2_2
    operator()(const typename AK::Polynomial_for_circles_2_2 & c, 
	       bool i) const
    { return AlgebraicFunctors::x_critical_point<AK>(c,i); }

    template <class OutputIterator>
    OutputIterator
    operator()(const typename AK::Polynomial_for_circles_2_2 & c, 
	       OutputIterator res) const
    { return AlgebraicFunctors::x_critical_points<AK>(c,res); }

  };

  template < class AK >
  class Y_critical_points
  {    
  public:
    
    typename AK::Root_for_circles_2_2
    operator()(const typename AK::Polynomial_for_circles_2_2 & c, 
	       bool i) const
    { return AlgebraicFunctors::y_critical_point<AK>(c,i); }

    template <class OutputIterator>
    OutputIterator
    operator()(const typename AK::Polynomial_for_circles_2_2 & c, 
	       OutputIterator res) const
    { return AlgebraicFunctors::y_critical_points<AK>(c,res); }

  };

  template <typename RT>
  class Compare_x
  {
  public:

    Comparison_result 
    operator()(const Root_for_circles_2_2<RT>& r1, 
	       const Root_for_circles_2_2<RT>& r2) const
    { return compare_x<RT>(r1, r2); }

  };

  template <typename RT>
  class Compare_y
  {
  public:

    Comparison_result 
    operator()(const Root_for_circles_2_2<RT>& r1, 
	       const Root_for_circles_2_2<RT>& r2) const
    { return compare_y<RT>(r1, r2); }

  };

  template <typename RT>
  class Compare_xy
  {
  public:

    Comparison_result 
    operator()(const Root_for_circles_2_2<RT>& r1, 
	       const Root_for_circles_2_2<RT>& r2) const
    { return compare_xy<RT>(r1, r2); }

  };

} // namespace AlgebraicFunctors

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_2_H
