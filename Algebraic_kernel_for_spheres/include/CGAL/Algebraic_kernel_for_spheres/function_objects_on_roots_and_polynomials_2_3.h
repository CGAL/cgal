// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL$
// $Id$
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_SPHERES_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_3_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_SPHERES_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/Algebraic_kernel_for_spheres/internal_functions_comparison_root_for_spheres_2_3.h>
#include <CGAL/Algebraic_kernel_for_spheres/internal_functions_on_roots_and_polynomials_2_3.h>
#include <CGAL/Algebraic_kernel_for_spheres/internal_functions_on_roots_and_polynomial_1_3_and_2_3.h>
#include <CGAL/Algebraic_kernel_for_spheres/internal_functions_on_roots_and_polynomials_1_3.h>
namespace CGAL {

namespace AlgebraicSphereFunctors {
  
  template < class AK >
  class Solve
  {

    typedef typename AK::Polynomial_for_spheres_2_3
      Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3
      Polynomial_1_3;
    typedef std::pair<
      Polynomial_for_spheres_2_3,
      Polynomial_1_3>       Equation_Circle;
    typedef typename AK::Polynomials_for_line_3 Polynomials_for_line_3;
    
  public:
    template < class OutputIterator >
      OutputIterator
      operator()
      (const Equation_Circle & e1,
       const Equation_Circle & e2,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Equation_Circle & e1,
       const Polynomials_for_line_3 & e2,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomials_for_line_3 & e1,
       const Equation_Circle & e2,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Equation_Circle & e1,
       const Polynomial_1_3 & e2,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomial_1_3 & e1,
       const Equation_Circle & e2,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> (e1, e2, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Equation_Circle & e1,
       const Polynomial_for_spheres_2_3 & e2,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomial_for_spheres_2_3 & e1,
       const Equation_Circle & e2,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomial_for_spheres_2_3 & e1,
       const Polynomial_1_3 & e2,
       const Polynomial_1_3 & e3,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, e3, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomial_for_spheres_2_3 & e1,
       const Polynomials_for_line_3 & e2,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomials_for_line_3 & e1,
       const Polynomial_for_spheres_2_3 & e2,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomial_for_spheres_2_3 & e1,
       const Polynomial_for_spheres_2_3 & e2,
       const Polynomial_1_3 & e3,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, e3, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomial_1_3 & e1,
       const Polynomial_for_spheres_2_3 & e2,
       const Polynomial_for_spheres_2_3 & e3,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, e3, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomial_1_3 & e1,
       const Polynomial_1_3 & e2,
       const Polynomial_for_spheres_2_3 & e3,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, e3, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomial_for_spheres_2_3 & e1,
       const Polynomial_for_spheres_2_3 & e2,
       const Polynomial_for_spheres_2_3 & e3,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, e3, res); }

    template < class OutputIterator >
      OutputIterator
      operator()
      (const Polynomial_1_3 & e1,
       const Polynomial_1_3 & e2,
       const Polynomial_1_3 & e3,
       OutputIterator res) const
      { return AlgebraicSphereFunctors::solve<AK> ( e1, e2, e3, res); }

    
    

  };

  template < class AK >
  class Construct_polynomial_for_spheres_2_3
  {
    typedef typename AK::RT                                        RT;
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    
  public:
    Polynomial_for_spheres_2_3
    operator()(const RT& xc, const RT& yc,const RT& zc, const RT& r_sq) const
    { return Polynomial_for_spheres_2_3(xc, yc, zc, r_sq); }
  };


  template < class AK >
  class Construct_polynomial_1_3
  {
    typedef typename AK::RT                                        RT;
    typedef typename AK::Polynomial_1_3 Polynomial_1_3;
    
  public:
    Polynomial_1_3
    operator()(const RT& a, const RT& b,const RT& c, const RT& d) const
    { return Polynomial_1_3(a, b, c, d); }
  };

  template < class AK >
  class Construct_polynomials_for_line_3
  {
    typedef typename AK::FT                                        FT;
    typedef typename AK::Polynomials_for_line_3 Polynomials_for_line_3;
    
  public:
    Polynomials_for_line_3
    operator()(const FT& a1, const FT& b1,
               const FT& a2, const FT& b2,
               const FT& a3, const FT& b3) const
    { return Polynomials_for_line_3(a1, b1, a2, b2, a3, b3); }
  };

  template < class AK >
  class Sign_at
  {
    typedef typename AK::Polynomial_1_3             Polynomial_1_3;
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Root_for_spheres_2_3       Root_for_spheres_2_3;

  public:
    typedef CGAL::Sign   result_type;

    result_type
    operator()( const Polynomial_for_spheres_2_3 & equation,
		const Root_for_spheres_2_3 & r ) const
    { return AlgebraicSphereFunctors::sign_at<AK>(equation, r); }

    result_type
    operator()( const Polynomial_1_3 & equation,
		const Root_for_spheres_2_3 & r ) const
    { return AlgebraicSphereFunctors::sign_at<AK>(equation, r); }

  };
    

  template < class AK >
  class X_critical_points
  {
    typedef typename AK::Root_of_2                  Root_of_2;
    typedef typename AK::Root_for_spheres_2_3       Root_for_spheres_2_3;
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3             Polynomial_1_3;

  public:
    typedef void         result_type;

    Root_for_spheres_2_3
    operator()(const Polynomial_for_spheres_2_3 & c, 
	       bool i) const
    { return AlgebraicSphereFunctors::x_critical_point<AK>(c,i); }

    template <class OutputIterator>
    OutputIterator
    operator()(const Polynomial_for_spheres_2_3 & c, 
	       OutputIterator res) const
    { return AlgebraicSphereFunctors::x_critical_points<AK>(c,res); }

    Root_for_spheres_2_3
    operator()(const std::pair< Polynomial_for_spheres_2_3, Polynomial_1_3 > & c, 
	       bool i) const
    { return AlgebraicSphereFunctors::x_critical_point<AK>(c,i); }

    template <class OutputIterator>
    OutputIterator
    operator()(const std::pair< Polynomial_for_spheres_2_3, Polynomial_1_3 > & c, 
	       OutputIterator res) const
    { return AlgebraicSphereFunctors::x_critical_points<AK>(c,res); }
  };

  template < class AK >
  class Y_critical_points
  {
    typedef typename AK::Root_of_2                  Root_of_2;
    typedef typename AK::Root_for_spheres_2_3       Root_for_spheres_2_3;
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3             Polynomial_1_3;

  public:
    typedef void         result_type;

    Root_for_spheres_2_3
    operator()(const Polynomial_for_spheres_2_3 & c, 
	       bool i) const
    { return AlgebraicSphereFunctors::y_critical_point<AK>(c,i); }

    template <class OutputIterator>
    OutputIterator
    operator()(const Polynomial_for_spheres_2_3 & c, 
	       OutputIterator res) const
    { return AlgebraicSphereFunctors::y_critical_points<AK>(c,res); }

    Root_for_spheres_2_3
    operator()(const std::pair< Polynomial_for_spheres_2_3, Polynomial_1_3 > & c, 
	       bool i) const
    { return AlgebraicSphereFunctors::y_critical_point<AK>(c,i); }

    template <class OutputIterator>
    OutputIterator
    operator()(const std::pair< Polynomial_for_spheres_2_3, Polynomial_1_3 > & c, 
	       OutputIterator res) const
    { return AlgebraicSphereFunctors::y_critical_points<AK>(c,res); }
  };

  template < class AK >
  class Z_critical_points
  {
    typedef typename AK::Root_of_2                  Root_of_2;
    typedef typename AK::Root_for_spheres_2_3       Root_for_spheres_2_3;
    typedef typename AK::Polynomial_for_spheres_2_3 Polynomial_for_spheres_2_3;
    typedef typename AK::Polynomial_1_3             Polynomial_1_3;

  public:
    typedef void         result_type;

    Root_for_spheres_2_3
    operator()(const Polynomial_for_spheres_2_3 & c, 
	       bool i) const
    { return AlgebraicSphereFunctors::z_critical_point<AK>(c,i); }

    template <class OutputIterator>
    OutputIterator
    operator()(const Polynomial_for_spheres_2_3 & c, 
	       OutputIterator res) const
    { return AlgebraicSphereFunctors::z_critical_points<AK>(c,res); }

    Root_for_spheres_2_3
    operator()(const std::pair< Polynomial_for_spheres_2_3, Polynomial_1_3 > & c, 
	       bool i) const
    { return AlgebraicSphereFunctors::z_critical_point<AK>(c,i); }

    template <class OutputIterator>
    OutputIterator
    operator()(const std::pair< Polynomial_for_spheres_2_3, Polynomial_1_3 > & c, 
	       OutputIterator res) const
    { return AlgebraicSphereFunctors::z_critical_points<AK>(c,res); }

  };
  
  template <typename RT>
  class Compare_x
  {
  public:
    Comparison_result 
    operator()(const Root_for_spheres_2_3<RT>& r1, 
	     const Root_for_spheres_2_3<RT>& r2) const
    { return AlgebraicSphereFunctors::compare_x<RT>(r1, r2); }

  };

  template <typename RT>
  class Compare_y
  {
  public:
     Comparison_result 
    operator()(const Root_for_spheres_2_3<RT>& r1, 
	     const Root_for_spheres_2_3<RT>& r2) const
    { return AlgebraicSphereFunctors::compare_y<RT>(r1, r2); }
  };

  template <typename RT>
  class Compare_z
  {
  public:
     Comparison_result 
    operator()(const Root_for_spheres_2_3<RT>& r1, 
	     const Root_for_spheres_2_3<RT>& r2) const
    { return AlgebraicSphereFunctors::compare_z<RT>(r1, r2); }
  };

  template <typename RT>
  class Compare_xy
  {
  public:
    Comparison_result 
    operator()(const Root_for_spheres_2_3<RT>& r1, 
	     const Root_for_spheres_2_3<RT>& r2) const
    { return AlgebraicSphereFunctors::compare_xy<RT>(r1, r2); }
  };

   template <typename RT>
  class Compare_xyz
  {
  public:
    Comparison_result 
    operator()(const Root_for_spheres_2_3<RT>& r1, 
	     const Root_for_spheres_2_3<RT>& r2) const
    { return AlgebraicSphereFunctors::compare_xyz<RT>(r1, r2); }
  };

} // namespace AlgebraicSphereFunctors

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_FOR_SPHERES_FUNCTION_OBJECTS_ON_ROOTS_AND_POLYNOMIALS_3_H
