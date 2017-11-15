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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//             Sylvain Pion
//             Pedro Machado
//             Julien Hazebrouck
//             Damien Leroy

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_SPHERES_2_3_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_SPHERES_2_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/Root_of_traits.h>
#include <CGAL/Polynomials_2_3.h>
#include <CGAL/Polynomials_1_3.h>
#include <CGAL/Polynomials_for_line_3.h>
#include <CGAL/Root_for_spheres_2_3.h>

#include <CGAL/Algebraic_kernel_for_spheres/function_objects_on_roots_and_polynomials_2_3.h>
#include <CGAL/global_functions_on_roots_and_polynomials_1_3.h>
#include <CGAL/global_functions_on_roots_and_polynomials_2_3.h>
#include <CGAL/global_functions_on_root_for_sphere_2_3.h>

namespace CGAL {

  template< class RT_ >
  struct Algebraic_kernel_for_spheres_2_3
  {
    typedef Algebraic_kernel_for_spheres_2_3<RT_>      Self;

    typedef RT_                                        RT;
    typedef typename Root_of_traits< RT >::RootOf_1    FT;
            
    typedef CGAL::Polynomials_for_line_3<FT>           Polynomials_for_line_3;
    typedef CGAL::Polynomial_for_spheres_2_3<FT>       Polynomial_for_spheres_2_3;
    typedef CGAL::Polynomial_1_3<FT>                   Polynomial_1_3;
    // problem RT / FT ?

    typedef typename Root_of_traits< RT >::RootOf_2    Root_of_2;
    typedef CGAL::Root_for_spheres_2_3< RT >           Root_for_spheres_2_3;

    typedef AlgebraicSphereFunctors::Construct_polynomial_for_spheres_2_3<Self>
                                         Construct_polynomial_for_spheres_2_3;
    typedef AlgebraicSphereFunctors::Construct_polynomial_1_3<Self>
                                         Construct_polynomial_1_3;
    typedef AlgebraicSphereFunctors::Construct_polynomials_for_line_3<Self>
                                         Construct_polynomials_for_line_3;

    typedef AlgebraicSphereFunctors::Solve<Self>             Solve;
    typedef AlgebraicSphereFunctors::Sign_at<Self>           Sign_at;
    typedef AlgebraicSphereFunctors::X_critical_points<Self> X_critical_points;
    typedef AlgebraicSphereFunctors::Y_critical_points<Self> Y_critical_points;
    typedef AlgebraicSphereFunctors::Z_critical_points<Self> Z_critical_points;
    typedef AlgebraicSphereFunctors::Compare_x<RT>           Compare_x;
    typedef AlgebraicSphereFunctors::Compare_y<RT>           Compare_y;
    typedef AlgebraicSphereFunctors::Compare_z<RT>           Compare_z;
    typedef AlgebraicSphereFunctors::Compare_xy<RT>          Compare_xy;
    typedef AlgebraicSphereFunctors::Compare_xyz<RT>         Compare_xyz;

    Construct_polynomial_for_spheres_2_3
                 construct_polynomial_for_spheres_2_3_object() const
    { return Construct_polynomial_for_spheres_2_3(); }

    Construct_polynomial_1_3
                 construct_polynomial_1_3_object() const
    { return Construct_polynomial_1_3(); }

    Construct_polynomials_for_line_3
                 construct_polynomials_for_line_3_object() const
    { return Construct_polynomials_for_line_3(); }

    Solve solve_object() const
    { return Solve(); }

    Sign_at sign_at_object() const
    { return Sign_at(); }

    X_critical_points x_critical_points_object() const
    { return X_critical_points(); }

    Y_critical_points y_critical_points_object() const
    { return Y_critical_points(); }

    Z_critical_points z_critical_points_object() const
    { return Z_critical_points(); }

    Compare_x compare_x_object() const
    { return Compare_x(); }
    
    Compare_y compare_y_object() const
    { return Compare_y(); }

    Compare_z compare_z_object() const
    { return Compare_z(); }
    
    Compare_xy compare_xy_object() const
    { return Compare_xy(); }
    
    Compare_xyz compare_xyz_object() const
    { return Compare_xyz(); }

  };

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_FOR_SPHERES_2_3_H
