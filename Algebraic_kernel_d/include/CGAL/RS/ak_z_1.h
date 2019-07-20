// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.

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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_AK_Z_1_H
#define CGAL_RS_AK_Z_1_H

#include <cstddef> // included only to define size_t
#include <CGAL/Polynomial_traits_d.h>
#include "algebraic_z_1.h"
#include "comparator_1.h"
#include "signat_1.h"
#include "functors_z_1.h"

// This file defines the "Z-algebraic kernel". This kind of kernel performs
// all the internal operations using an integer polynomial type (the name
// "Z" comes from there). For this, a converter functor (passed as a
// template parameter) is used, which converts the input polynomial to the
// integer representation.

namespace CGAL{
namespace RS_AK1{

template <class ExtPolynomial_,
          class IntPolynomial_,
          class PolConverter_,
          class Bound_,
          class Isolator_,
          class Refiner_,
          class Ptraits_=CGAL::Polynomial_traits_d<ExtPolynomial_>,
          class ZPtraits_=CGAL::Polynomial_traits_d<IntPolynomial_> >
class Algebraic_kernel_z_1{
        public:
        typedef ExtPolynomial_                          Polynomial_1;
        typedef IntPolynomial_                          ZPolynomial_1;
        typedef PolConverter_                           PolConverter;
        typedef typename Polynomial_1::NT               Coefficient;
        typedef Bound_                                  Bound;
        private:
        typedef Isolator_                               Isolator;
        typedef Refiner_                                Refiner;
        typedef Ptraits_                                Ptraits;
        typedef ZPtraits_                               ZPtraits;
        typedef CGAL::RS_AK1::Signat_1<ZPolynomial_1,Bound>
                                                        Signat;
        typedef CGAL::RS_AK1::Simple_comparator_1<ZPolynomial_1,
                                                  Bound,
                                                  Refiner,
                                                  Signat,
                                                  ZPtraits>
                                                        Comparator;
        public:
        typedef CGAL::RS_AK1::Algebraic_z_1<Polynomial_1,
                                            ZPolynomial_1,
                                            Bound,
                                            Refiner,
                                            Comparator,
                                            Ptraits,
                                            ZPtraits>
                                                        Algebraic_real_1;
        typedef size_t                                  size_type;
        typedef unsigned                                Multiplicity_type;

        // default constructor and destructor
        public:
        Algebraic_kernel_z_1(){};
        ~Algebraic_kernel_z_1(){};

        // functors from the CGAL concept
        public:
        typedef CGAL::RS_AK1::Construct_algebraic_real_z_1<Polynomial_1,
                                                           ZPolynomial_1,
                                                           PolConverter,
                                                           Algebraic_real_1,
                                                           Bound,
                                                           Coefficient,
                                                           Isolator>
                                                Construct_algebraic_real_1;
        typedef CGAL::RS_AK1::Compute_polynomial_z_1<Polynomial_1,
                                                     Algebraic_real_1>
                                                        Compute_polynomial_1;
        typedef CGAL::RS_AK1::Isolate_z_1<Polynomial_1,
                                          ZPolynomial_1,
                                          PolConverter,
                                          Bound,
                                          Algebraic_real_1,
                                          Isolator,
                                          Comparator,
                                          Signat,
                                          Ptraits,
                                          ZPtraits>
                                                        Isolate_1;
        typedef typename Ptraits::Is_square_free        Is_square_free_1;
        typedef typename Ptraits::Make_square_free      Make_square_free_1;
        typedef typename Ptraits::Square_free_factorize Square_free_factorize_1;
        typedef CGAL::RS_AK1::Is_coprime_z_1<Polynomial_1,Ptraits>
                                                        Is_coprime_1;
        typedef CGAL::RS_AK1::Make_coprime_z_1<Polynomial_1,Ptraits>
                                                        Make_coprime_1;
        typedef CGAL::RS_AK1::Solve_z_1<Polynomial_1,
                                        ZPolynomial_1,
                                        PolConverter,
                                        Bound,
                                        Algebraic_real_1,
                                        Isolator,
                                        Signat,
                                        Ptraits,
                                        ZPtraits>
                                                        Solve_1;
        typedef CGAL::RS_AK1::Number_of_solutions_z_1<Polynomial_1,
                                                      ZPolynomial_1,
                                                      PolConverter,
                                                      Isolator>
                                                        Number_of_solutions_1;

        typedef CGAL::RS_AK1::Sign_at_z_1<Polynomial_1,
                                          ZPolynomial_1,
                                          PolConverter,
                                          Bound,
                                          Algebraic_real_1,
                                          Refiner,
                                          Signat,
                                          Ptraits,
                                          ZPtraits>
                                                        Sign_at_1;
        typedef CGAL::RS_AK1::Is_zero_at_z_1<Polynomial_1,
                                             ZPolynomial_1,
                                             PolConverter,
                                             Bound,
                                             Algebraic_real_1,
                                             Refiner,
                                             Signat,
                                             Ptraits,
                                             ZPtraits>
                                                        Is_zero_at_1;
        typedef CGAL::RS_AK1::Compare_z_1<Algebraic_real_1,
                                          Bound,
                                          Comparator>
                                                        Compare_1;
        typedef CGAL::RS_AK1::Bound_between_z_1<Algebraic_real_1,
                                                Bound,
                                                Comparator>
                                                        Bound_between_1;
        typedef CGAL::RS_AK1::Approximate_absolute_z_1<Polynomial_1,
                                                       Bound,
                                                       Algebraic_real_1,
                                                       Refiner>
                                                Approximate_absolute_1;
        typedef CGAL::RS_AK1::Approximate_relative_z_1<Polynomial_1,
                                                       Bound,
                                                       Algebraic_real_1,
                                                       Refiner>
                                                Approximate_relative_1;

#define CREATE_FUNCTION_OBJECT(T,N) \
        T N##_object()const{return T();}
        CREATE_FUNCTION_OBJECT(Construct_algebraic_real_1,
                               construct_algebraic_real_1)
        CREATE_FUNCTION_OBJECT(Compute_polynomial_1,
                               compute_polynomial_1)
        CREATE_FUNCTION_OBJECT(Isolate_1,
                               isolate_1)
        CREATE_FUNCTION_OBJECT(Is_square_free_1,
                               is_square_free_1)
        CREATE_FUNCTION_OBJECT(Make_square_free_1,
                               make_square_free_1)
        CREATE_FUNCTION_OBJECT(Square_free_factorize_1,
                               square_free_factorize_1)
        CREATE_FUNCTION_OBJECT(Is_coprime_1,
                               is_coprime_1)
        CREATE_FUNCTION_OBJECT(Make_coprime_1,
                               make_coprime_1)
        CREATE_FUNCTION_OBJECT(Solve_1,
                               solve_1)
        CREATE_FUNCTION_OBJECT(Number_of_solutions_1,
                               number_of_solutions_1)
        CREATE_FUNCTION_OBJECT(Sign_at_1,
                               sign_at_1)
        CREATE_FUNCTION_OBJECT(Is_zero_at_1,
                               is_zero_at_1)
        CREATE_FUNCTION_OBJECT(Compare_1,
                               compare_1)
        CREATE_FUNCTION_OBJECT(Bound_between_1,
                               bound_between_1)
        CREATE_FUNCTION_OBJECT(Approximate_absolute_1,
                               approximate_absolute_1)
        CREATE_FUNCTION_OBJECT(Approximate_relative_1,
                               approximate_relative_1)
#undef CREATE_FUNCTION_OBJECT

}; // class Algebraic_kernel_z_1

} // namespace RS_AK1
} // namespace CGAL

#endif // CGAL_RS_AK_Z_1_H
