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
// $URL:svn+ssh://pmachado@scm.gforge.inria.fr/svn/cgal/trunk/Algebraic_kernel_for_circles/include/CGAL/Algebraic_kernel_for_circles_2_2.h $
// $Id:Algebraic_kernel_for_circles_2_2.h 33723 2006-08-25 08:36:44 +0000 (Fri, 25 Aug 2006) pmachado $
//
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_2_2_H

#include <CGAL/Quotient.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Polynomials_2_2.h>
#include <CGAL/Polynomials_1_2.h>
#include <CGAL/Root_for_circles_2_2.h>

#include <CGAL/Algebraic_kernel_for_circles/function_objects_on_roots_and_polynomials_2_2.h>

namespace CGAL {

  template< class RT_ >
  struct Algebraic_kernel_for_circles_2_2
  {
    typedef Algebraic_kernel_for_circles_2_2<RT_>      Self;

    typedef RT_                                        RT;
    typedef typename Root_of_traits< RT >::RootOf_1    FT;

    typedef CGAL::Polynomial_1_2<RT>               Polynomial_1_2; 
    typedef CGAL::Polynomial_for_circles_2_2<RT>   Polynomial_for_circles_2_2; 
    // problem RT / FT ?

    typedef typename Root_of_traits< RT >::RootOf_2  Root_of_2;
    typedef CGAL::Root_for_circles_2_2< RT >         Root_for_circles_2_2;

    typedef AlgebraicFunctors::Construct_polynomial_1_2<Self>
                                         Construct_polynomial_1_2;
    typedef AlgebraicFunctors::Construct_polynomial_for_circles_2_2<Self>
                                         Construct_polynomial_for_circles_2_2;

    typedef AlgebraicFunctors::Solve<Self>             Solve;
    typedef AlgebraicFunctors::Sign_at<Self>           Sign_at;
    typedef AlgebraicFunctors::X_critical_points<Self> X_critical_points;
    typedef AlgebraicFunctors::Y_critical_points<Self> Y_critical_points;
    typedef AlgebraicFunctors::Compare_x<Self>         Compare_x;
    typedef AlgebraicFunctors::Compare_y<Self>         Compare_y;
    typedef AlgebraicFunctors::Compare_xy<Self>        Compare_xy;

    Construct_polynomial_1_2
                 construct_polynomial_1_2_object() const
    { return Construct_polynomial_1_2(); }

    Construct_polynomial_for_circles_2_2
                 construct_polynomial_for_circles_2_2_object() const
    { return Construct_polynomial_for_circles_2_2(); }

    Solve solve_object() const
    { return Solve(); }

    Sign_at sign_at_object() const
    { return Sign_at(); }

    X_critical_points x_critical_points_object() const
    { return X_critical_points(); }

    Y_critical_points y_critical_points_object() const
    { return Y_critical_points(); }

    Compare_x compare_x_object() const
    { return Compare_x(); }
    
    Compare_y compare_y_object() const
    { return Compare_y(); }
    
    Compare_xy compare_xy_object() const
    { return Compare_xy(); }

  };

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_2_2_H
