// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_2_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_2_H

#include <CGAL/Algebraic_kernel_for_circles/internal_functions_on_roots_and_polynomials_2_2.h>
#include <CGAL/Algebraic_kernel_for_circles/internal_functions_on_roots_and_polynomial_1_2_and_2_2.h>
#include <CGAL/Algebraic_kernel_for_circles/internal_functions_comparison_root_for_circles_2_2.h>

namespace CGAL {

namespace AlgebraicFunctors {
  
  template < class AK >
  class Solve
  {
    typedef typename AK::Polynomial_for_circles_2_2 Equation_Circle;
    typedef typename AK::Polynomial_1_2             Equation_Line;

  public:
    typedef void         result_type; 

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
    operator()(const Equation_Circle & e1, 
	       const Equation_Line & e2, 
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
    
    typedef Polynomial_for_circles_2_2 result_type;

    result_type
    operator()(const RT& xc, const RT& yc, const RT& r_sq) const
    { return Polynomial_for_circles_2_2(xc, yc, r_sq); }

  };

  template < class AK >
  class Construct_polynomial_1_2
  {
    typedef typename AK::RT             RT;
    typedef typename AK::Polynomial_1_2 Polynomial_1_2;
    
  public:
    
    typedef Polynomial_1_2 result_type;

    result_type
    operator()( const RT& a, const RT& b, const RT& c) const
    { return Polynomial_1_2(a, b, c); }

  };

  template < class AK >
  class Sign_at
  {    
    typedef typename AK::Polynomial_1_2             Polynomial_1_2;
    typedef typename AK::Polynomial_for_circles_2_2 Polynomial_for_circles_2_2;
    typedef typename AK::Root_for_circles_2_2       Root_for_circles_2_2;

  public:
    typedef CGAL::Sign   result_type;

    result_type
    operator()( const Polynomial_for_circles_2_2 & equation,
		const Root_for_circles_2_2 & r ) const
    { return AlgebraicFunctors::sign_at<AK>(equation, r); }

    result_type
    operator()( const Polynomial_1_2 & equation,
		const Root_for_circles_2_2 & r ) const
    { return AlgebraicFunctors::sign_at<AK>(equation, r); }

  };
    
  template < class AK >
  class X_critical_points
  {    
    typedef typename AK::Root_for_circles_2_2       Root_for_circles_2_2;
    typedef typename AK::Polynomial_for_circles_2_2 Polynomial_for_circles_2_2;

  public:
    typedef void         result_type;

    Root_for_circles_2_2
    operator()(const Polynomial_for_circles_2_2 & c, 
	       bool i) const
    { return AlgebraicFunctors::x_critical_point<AK>(c,i); }

    template <class OutputIterator>
    OutputIterator
    operator()(const Polynomial_for_circles_2_2 & c, 
	       OutputIterator res) const
    { return AlgebraicFunctors::x_critical_points<AK>(c,res); }

  };

  template < class AK >
  class Y_critical_points
  {    
    typedef typename AK::Root_for_circles_2_2       Root_for_circles_2_2;
    typedef typename AK::Polynomial_for_circles_2_2 Polynomial_for_circles_2_2;

  public:
    typedef void         result_type;

    Root_for_circles_2_2
    operator()(const Polynomial_for_circles_2_2 & c, 
	       bool i) const
    { return AlgebraicFunctors::y_critical_point<AK>(c,i); }

    template <class OutputIterator>
    OutputIterator
    operator()(const Polynomial_for_circles_2_2 & c, 
	       OutputIterator res) const
    { return AlgebraicFunctors::y_critical_points<AK>(c,res); }

  };

  template < class AK >
  class Compare_x
  {
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
    typedef typename AK::RT                   RT;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator()(const Root_for_circles_2_2& r1, 
	       const Root_for_circles_2_2& r2) const
    { return AlgebraicFunctors::compare_x<RT>(r1, r2); }

  };

  template < class AK >
  class Compare_y
  {
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
    typedef typename AK::RT                   RT;
    
  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator()(const Root_for_circles_2_2& r1, 
	       const Root_for_circles_2_2& r2) const
    { return AlgebraicFunctors::compare_y<RT>(r1, r2); }

  };

  template < class AK >
  class Compare_xy
  {
    typedef typename AK::Root_for_circles_2_2 Root_for_circles_2_2;
    typedef typename AK::RT                   RT;

  public:
    typedef CGAL::Comparison_result result_type;

    result_type
    operator()(const Root_for_circles_2_2& r1, 
	       const Root_for_circles_2_2& r2) const
    { return AlgebraicFunctors::compare_xy<RT>(r1, r2); }

  };

} // namespace AlgebraicFunctors

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_2_H
