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
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_RS3_REFINER_1_H
#define CGAL_RS_RS3_REFINER_1_H

#include <CGAL/Polynomial_traits_d.h>
#include "rs2_calls.h"
#include <rs3_fncts.h>
#include "Gmpfr_make_unique.h"

namespace CGAL{
namespace RS3{

template <class Polynomial_,class Bound_>
struct RS3_refiner_1{
        void operator()(const Polynomial_&,Bound_&,Bound_&,int,int=0);
}; // class RS3_refiner_1

template <class Polynomial_,class Bound_>
void
RS3_refiner_1<Polynomial_,Bound_>::
operator()(const Polynomial_&,Bound_&,Bound_&,int,int){
        CGAL_error_msg("RS3 refiner not implemented for these types");
        return;
}

template<>
void
RS3_refiner_1<Polynomial<Gmpz>,Gmpfr>::
operator()
(const Polynomial<Gmpz> &pol,Gmpfr &left,Gmpfr &right,int prec,int k){
        typedef Polynomial<Gmpz>                        Polynomial;
        typedef Polynomial_traits_d<Polynomial>         Ptraits;
        typedef Ptraits::Degree                         Degree;
        CGAL_precondition(left<=right);
        // TODO: add precondition to check whether the interval is a point
        // or the evaluations on its endpoints have different signs
        //std::cout<<"refining ["<<left<<","<<right<<"]"<<std::endl;
        int deg=Degree()(pol);
        mpz_t* coefficients=(mpz_t*)malloc((deg+1)*sizeof(mpz_t));
        __mpfi_struct interval;
        // Make sure the endpoints do not share references.
        CGAL_RS_GMPFR_MAKE_UNIQUE(left,temp_left);
        CGAL_RS_GMPFR_MAKE_UNIQUE(right,temp_right);
        // Construct the mpfi_t interval which will be refined.
        interval.left=*(left.fr());
        interval.right=*(right.fr());
        // Construct the polynomial which will be refined (copy pointers).
        for(int i=0;i<=deg;++i)
                coefficients[i][0]=*(pol[i].mpz());
        // Call RS.
        RS2::RS2_calls::init_solver();
        rs3_refine_u_root(&interval,
                          coefficients,
                          deg,
                          prec+CGAL::max(left.get_precision(),
                                         right.get_precision()),
                          k,
                          k);
        // Clear variables.
        free(coefficients);
        mpfr_clear(left.fr());
        mpfr_clear(right.fr());
        // Copy results back to the Gmpfr endpoints.
        mpfr_custom_init_set(left.fr(),
                             mpfr_custom_get_kind(&interval.left),
                             mpfr_custom_get_exp(&interval.left),
                             mpfr_get_prec(&interval.left),
                             mpfr_custom_get_mantissa(&interval.left));
        mpfr_custom_init_set(right.fr(),
                             mpfr_custom_get_kind(&interval.right),
                             mpfr_custom_get_exp(&interval.right),
                             mpfr_get_prec(&interval.right),
                             mpfr_custom_get_mantissa(&interval.right));
        CGAL_postcondition(left<=right);
        //std::cout<<"ref root is ["<<left<<","<<right<<"]"<<std::endl;
        return;
}

template<>
void
RS3_refiner_1<Polynomial<Gmpq>,Gmpfr>::
operator()
(const Polynomial<Gmpq> &qpol,Gmpfr &left,Gmpfr &right,int prec,int k){
        typedef Polynomial<Gmpz>                        ZPolynomial;
        typedef Polynomial_traits_d<ZPolynomial>        ZPtraits;
        typedef ZPtraits::Degree                        ZDegree;
        CGAL_precondition(left<=right);
        // TODO: add precondition to check whether the interval is a point
        // or the evaluations on its endpoints have different signs
        //std::cout<<"refining ["<<left<<","<<right<<"]"<<std::endl;
        // Construct a Gmpz polynomial from the original Gmpq polynomial.
        Polynomial<Gmpz> zpol=CGAL::RS_AK1::Polynomial_converter_1<
                                        CGAL::Polynomial<Gmpq>,
                                        CGAL::Polynomial<Gmpz> >()(qpol);
        int deg=ZDegree()(zpol);
        mpz_t* coefficients=(mpz_t*)malloc((deg+1)*sizeof(mpz_t));
        __mpfi_struct interval;
        // Make sure the endpoints do not share references.
        CGAL_RS_GMPFR_MAKE_UNIQUE(left,temp_left);
        CGAL_RS_GMPFR_MAKE_UNIQUE(right,temp_right);
        // Construct the mpfi_t interval which will be refined.
        interval.left=*(left.fr());
        interval.right=*(right.fr());
        // Construct the polynomial which will be refined (copy pointers).
        for(int i=0;i<=deg;++i)
                coefficients[i][0]=*(zpol[i].mpz());
        // Call RS.
        RS2::RS2_calls::init_solver();
        rs3_refine_u_root(&interval,
                          coefficients,
                          deg,
                          prec+CGAL::max(left.get_precision(),
                                         right.get_precision()),
                          k,
                          k);
        // Clear variables.
        free(coefficients);
        mpfr_clear(left.fr());
        mpfr_clear(right.fr());
        // Copy results back to the Gmpfr endpoints.
        mpfr_custom_init_set(left.fr(),
                             mpfr_custom_get_kind(&interval.left),
                             mpfr_custom_get_exp(&interval.left),
                             mpfr_get_prec(&interval.left),
                             mpfr_custom_get_mantissa(&interval.left));
        mpfr_custom_init_set(right.fr(),
                             mpfr_custom_get_kind(&interval.right),
                             mpfr_custom_get_exp(&interval.right),
                             mpfr_get_prec(&interval.right),
                             mpfr_custom_get_mantissa(&interval.right));
        CGAL_postcondition(left<=right);
        //std::cout<<"ref root is ["<<left<<","<<right<<"]"<<std::endl;
        return;
}

} // namespace RS3
} // namespace CGAL

#endif // CGAL_RS_RS3_REFINER_1_H
