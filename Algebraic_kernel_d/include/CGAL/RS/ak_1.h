// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_AK_1_H
#define CGAL_RS_AK_1_H

#include <cstddef> // included only to define size_t
#include <CGAL/Polynomial_traits_d.h>
#include "algebraic_1.h"
#include "comparator_1.h"
#include "signat_1.h"
#include "functors_1.h"

namespace CGAL{
namespace RS_AK1{

template <class Polynomial_,
          class Bound_,
          class Isolator_,
          class Refiner_,
          class Ptraits_=CGAL::Polynomial_traits_d<Polynomial_> >
class Algebraic_kernel_1{
        public:
        typedef Polynomial_                             Polynomial_1;
        typedef typename Polynomial_1::NT               Coefficient;
        typedef Bound_                                  Bound;
        private:
        typedef Isolator_                               Isolator;
        typedef Refiner_                                Refiner;
        typedef Ptraits_                                Ptraits;
        typedef CGAL::RS_AK1::Signat_1<Polynomial_1,Bound>
                                                        Signat;
        typedef CGAL::RS_AK1::Simple_comparator_1<Polynomial_1,
                                                  Bound,
                                                  Refiner,
                                                  Signat,
                                                  Ptraits>
                                                        Comparator;
        public:
        typedef CGAL::RS_AK1::Algebraic_1<Polynomial_1,
                                          Bound,
                                          Refiner,
                                          Comparator,
                                          Ptraits>
                                                        Algebraic_real_1;
        typedef size_t                                  size_type;
        typedef unsigned                                Multiplicity_type;

        // default constructor and destructor
        public:
        Algebraic_kernel_1(){};
        ~Algebraic_kernel_1(){};

        // functors from the CGAL concept
        public:
        typedef CGAL::RS_AK1::Construct_algebraic_real_1<Polynomial_1,
                                                         Algebraic_real_1,
                                                         Bound,
                                                         Coefficient,
                                                         Isolator>
                                                Construct_algebraic_real_1;
        typedef CGAL::RS_AK1::Compute_polynomial_1<Polynomial_1,
                                                   Algebraic_real_1>
                                                        Compute_polynomial_1;
        typedef CGAL::RS_AK1::Isolate_1<Polynomial_1,
                                        Bound,
                                        Algebraic_real_1,
                                        Isolator,
                                        Comparator,
                                        Signat,
                                        Ptraits>
                                                        Isolate_1;
        typedef typename Ptraits::Is_square_free        Is_square_free_1;
        typedef typename Ptraits::Make_square_free      Make_square_free_1;
        typedef typename Ptraits::Square_free_factorize
                                                Square_free_factorize_1;
        typedef CGAL::RS_AK1::Is_coprime_1<Polynomial_1,Ptraits>
                                                        Is_coprime_1;
        typedef CGAL::RS_AK1::Make_coprime_1<Polynomial_1,Ptraits>
                                                        Make_coprime_1;
        typedef CGAL::RS_AK1::Solve_1<Polynomial_1,
                                         Bound,
                                         Algebraic_real_1,
                                         Isolator,
                                         Signat,
                                         Ptraits>
                                                        Solve_1;
        typedef CGAL::RS_AK1::Number_of_solutions_1<Polynomial_1,Isolator>
                                                        Number_of_solutions_1;

        typedef CGAL::RS_AK1::Sign_at_1<Polynomial_1,
                                        Bound,
                                        Algebraic_real_1,
                                        Refiner,
                                        Signat,
                                        Ptraits>
                                                        Sign_at_1;
        typedef CGAL::RS_AK1::Is_zero_at_1<Polynomial_1,
                                           Bound,
                                           Algebraic_real_1,
                                           Refiner,
                                           Signat,
                                           Ptraits>
                                                        Is_zero_at_1;
        typedef CGAL::RS_AK1::Compare_1<Algebraic_real_1,
                                        Bound,
                                        Comparator>
                                                        Compare_1;
        typedef CGAL::RS_AK1::Bound_between_1<Algebraic_real_1,
                                              Bound,
                                              Comparator>
                                                        Bound_between_1;
        typedef CGAL::RS_AK1::Approximate_absolute_1<Polynomial_1,
                                                     Bound,
                                                     Algebraic_real_1,
                                                     Refiner>
                                                Approximate_absolute_1;
        typedef CGAL::RS_AK1::Approximate_relative_1<Polynomial_1,
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

}; // class Algebraic_kernel_1

} // namespace RS_AK1
} // namespace CGAL

#endif // CGAL_RS_AK_1_H
