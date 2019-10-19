// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_RS2_ISOLATOR_1_H
#define CGAL_RS_RS2_ISOLATOR_1_H

#include "rs2_calls.h"
#include "polynomial_converter_1.h"
#include <CGAL/Gmpfi.h>
#include <vector>

namespace CGAL{
namespace RS2{

template <class Polynomial_,class Bound_>
class RS2_isolator_1{
        public:
        typedef Polynomial_                             Polynomial;
        typedef Bound_                                  Bound;
        private:
        typedef Gmpfi                                   Interval;
        public:
        RS2_isolator_1(const Polynomial&);
        Polynomial polynomial()const;
        int number_of_real_roots()const;
        bool is_exact_root(int i)const;
        Bound left_bound(int i)const;
        Bound right_bound(int i)const;
        private:
        Polynomial _polynomial;
        std::vector<Interval> _real_roots;
};

template <class Polynomial_,class Bound_>
RS2_isolator_1<Polynomial_,Bound_>::
RS2_isolator_1(const Polynomial_ &){
        CGAL_error_msg("not implemented for these polynomial/bound types");
}

// Isolator of Gmpz polynomials and Gmpfr bounds (the best case for RS).
template <>
RS2_isolator_1<CGAL::Polynomial<CGAL::Gmpz>,Gmpfr>::
RS2_isolator_1(const CGAL::Polynomial<CGAL::Gmpz> &p):_polynomial(p){
        //mpz_t *coeffs=(mpz_t*)malloc((degree+1)*sizeof(mpz_t));
        //mpfi_ptr *intervals_mpfi=(mpfi_ptr*)malloc(degree*sizeof(mpfi_ptr));
        //for(unsigned i=0;i<=degree;++i)
        //        coeffs[i][0]=*(p[i].mpz());
        RS2::RS2_calls::init_solver();
        RS2::RS2_calls::create_rs_upoly(p,rs_get_default_up());
        //free(coeffs);
        set_rs_precisol(0);
        set_rs_verbose(0);
        rs_run_algo((char*)"UISOLE");
        RS2::RS2_calls::insert_roots(std::back_inserter(_real_roots));
}

// Isolator of Gmpq polynomials and Gmpfr bounds. The polynomial must be
// converted to a Gmpz one with the same roots.
template <>
RS2_isolator_1<CGAL::Polynomial<CGAL::Gmpq>,Gmpfr>::
RS2_isolator_1(const CGAL::Polynomial<CGAL::Gmpq> &p):_polynomial(p){
        RS2::RS2_calls::init_solver();
        RS2::RS2_calls::create_rs_upoly(
                        CGAL::RS_AK1::Polynomial_converter_1<
                                        CGAL::Polynomial<Gmpq>,
                                        CGAL::Polynomial<Gmpz> >()(p),
                        rs_get_default_up());
        set_rs_precisol(0);
        set_rs_verbose(0);
        rs_run_algo((char*)"UISOLE");
        RS2::RS2_calls::insert_roots(std::back_inserter(_real_roots));
}

template <class Polynomial_,class Bound_>
Polynomial_
RS2_isolator_1<Polynomial_,Bound_>::polynomial()const{
        return _polynomial;
}

template <class Polynomial_,class Bound_>
int
RS2_isolator_1<Polynomial_,Bound_>::number_of_real_roots()const{
        return _real_roots.size();
}

template <class Polynomial_,class Bound_>
bool
RS2_isolator_1<Polynomial_,Bound_>::is_exact_root(int i)const{
        return _real_roots[i].inf()==_real_roots[i].sup();
}

template <class Polynomial_,class Bound_>
Bound_
RS2_isolator_1<Polynomial_,Bound_>::left_bound(int i)const{
        return _real_roots[i].inf();
}

template <class Polynomial_,class Bound_>
Bound_
RS2_isolator_1<Polynomial_,Bound_>::right_bound(int i)const{
        return _real_roots[i].sup();
}

} // namespace RS2
} // namespace CGAL

#endif // CGAL_RS_RS2_ISOLATOR_1_H
