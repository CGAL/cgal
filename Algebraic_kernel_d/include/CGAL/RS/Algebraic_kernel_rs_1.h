// Copyright (c) 2006-2010 Inria Lorraine (France). All rights reserved.
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
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_ALGEBRAIC_KERNEL_RS_1
#define CGAL_RS_ALGEBRAIC_KERNEL_RS_1

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/RS/functors_1.h>

template <class _C,
          class _G=
#ifdef CGAL_RS_USE_UGCD
                   CGAL::Modgcd_1
#else
                   CGAL::Rsgcd_1
#endif
>
struct Algebraic_kernel_rs_1{

        typedef _C                                      Coefficient;
        typedef _G                                      Gcd;
        typedef CGAL::Polynomial<Coefficient>           Polynomial_1;
        typedef RSFunctors::Algebraic                   Algebraic_real_1;
        typedef RSFunctors::Bound                       Bound;
        typedef RSFunctors::Multiplicity                Multiplicity_type;

        // constructor: we must initialize RS just a time, so this is a
        // good time to do it
        Algebraic_kernel_rs_1(){CGAL::init_solver();};
        ~Algebraic_kernel_rs_1(){CGAL::reset_solver();};

        typedef RSFunctors::Construct_alg_1<Polynomial_1,Coefficient,Gcd>
                                                Construct_algebraic_real_1;
        typedef RSFunctors::Compute_polynomial_1<Polynomial_1>
                                                        Compute_polynomial_1;
        typedef RSFunctors::Isolate_1<Polynomial_1,Gcd> Isolate_1;
        typedef RSFunctors::Is_square_free_1<Polynomial_1,Gcd>
                                                        Is_square_free_1;
        typedef RSFunctors::Make_square_free_1<Polynomial_1,Gcd>
                                                        Make_square_free_1;
        typedef RSFunctors::Square_free_factorize_1<Polynomial_1,Gcd>
                                                        Square_free_factorize_1;
        typedef RSFunctors::Is_coprime_1<Polynomial_1,Gcd>
                                                        Is_coprime_1;
        typedef RSFunctors::Make_coprime_1<Polynomial_1,Gcd>
                                                        Make_coprime_1;
        typedef RSFunctors::Solve_1<Polynomial_1,Gcd>   Solve_1;
        typedef RSFunctors::Number_of_solutions_1<Polynomial_1>
                                                        Number_of_solutions_1;
        typedef RSFunctors::Sign_at_1<Polynomial_1,Gcd> Sign_at_1;
        typedef RSFunctors::Is_zero_at_1<Polynomial_1,Gcd>
                                                        Is_zero_at_1;
        typedef RSFunctors::Compare_1<Gcd>              Compare_1;
        typedef RSFunctors::Bound_between_1<Gcd>        Bound_between_1;
        typedef RSFunctors::Approximate_absolute_1      Approximate_absolute_1;
        typedef RSFunctors::Approximate_relative_1      Approximate_relative_1;

#define CGALRS_CREATE_FUNCTION_OBJECT(T,N) \
        T N##_object()const{return T();}
        CGALRS_CREATE_FUNCTION_OBJECT(Construct_algebraic_real_1,
                                      construct_algebraic_real_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Compute_polynomial_1,compute_polynomial_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Isolate_1,isolate_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Is_square_free_1,is_square_free_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Make_square_free_1,make_square_free_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Square_free_factorize_1,
                                      square_free_factorize_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Is_coprime_1,is_coprime_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Make_coprime_1,make_coprime_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Solve_1,solve_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Number_of_solutions_1,
                                      number_of_solutions_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Sign_at_1,sign_at_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Is_zero_at_1,is_zero_at_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Compare_1,compare_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Bound_between_1,bound_between_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Approximate_absolute_1,
                                      approximate_absolute_1)
        CGALRS_CREATE_FUNCTION_OBJECT(Approximate_relative_1,
                                      approximate_relative_1)
#undef CGALRS_CREATE_FUNCTION_OBJECT
};  // Algebraic_kernel_d_1_RS

#endif  // CGAL_RS_ALGEBRAIC_KERNEL_RS_1
