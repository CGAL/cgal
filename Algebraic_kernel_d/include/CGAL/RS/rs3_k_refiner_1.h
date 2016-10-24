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

#ifndef CGAL_RS_RS3_K_REFINER_1_H
#define CGAL_RS_RS3_K_REFINER_1_H

#include <CGAL/Polynomial_traits_d.h>
#include "polynomial_converter_1.h"
#include "rs2_calls.h"
#include <rs3_fncts.h>
#include "Gmpfr_make_unique.h"

// If we want assertions, we need to evaluate.
#ifndef CGAL_NO_PRECONDITIONS
#include "signat_1.h"
#endif

namespace CGAL{
namespace RS3{

template <class Polynomial_,class Bound_>
struct RS3_k_refiner_1{
        void operator()(const Polynomial_&,Bound_&,Bound_&,int);
}; // class RS3_k_refiner_1

template <class Polynomial_,class Bound_>
void
RS3_k_refiner_1<Polynomial_,Bound_>::
operator()(const Polynomial_&,Bound_&,Bound_&,int){
        CGAL_error_msg("RS3 k-refiner not implemented for these types");
        return;
}

template<>
void
RS3_k_refiner_1<Polynomial<Gmpz>,Gmpfr>::
operator()
(const Polynomial<Gmpz> &pol,Gmpfr &left,Gmpfr &right,int prec){
        typedef Polynomial<Gmpz>                        Polynomial;
        typedef Polynomial_traits_d<Polynomial>         Ptraits;
        typedef Ptraits::Degree                         Degree;
        CGAL_precondition(left<=right);
#ifndef CGAL_NO_PRECONDITIONS
        typedef Ptraits::Make_square_free               Sfpart;
        typedef CGAL::RS_AK1::Signat_1<Polynomial,Gmpfr>
                                                        Signat;
        Polynomial sfpp=Sfpart()(pol);
        Signat signof(sfpp);
        CGAL::Sign sl=signof(left);
        CGAL_precondition(sl!=signof(right)||(left==right&&sl==ZERO));
#endif
        //std::cout<<"refining ["<<left<<","<<right<<"]"<<std::endl;
        int deg=Degree()(pol);
        mpz_t* coefficients=(mpz_t*)malloc((deg+1)*sizeof(mpz_t));
        __mpfi_struct interval;
        CGAL_RS_GMPFR_MAKE_UNIQUE(left,temp_left);
        CGAL_RS_GMPFR_MAKE_UNIQUE(right,temp_right);
        interval.left=*(left.fr());
        interval.right=*(right.fr());
        for(int i=0;i<=deg;++i)
                coefficients[i][0]=*(pol[i].mpz());
        RS2::RS2_calls::init_solver();
        rs3_refine_u_root(&interval,
                          coefficients,
                          deg,
                          prec+CGAL::max(left.get_precision(),
                                         right.get_precision()),
                          1,
                          1);
        free(coefficients);
        mpfr_clear(left.fr());
        mpfr_clear(right.fr());
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
RS3_k_refiner_1<Polynomial<Gmpq>,Gmpfr>::
operator()
(const Polynomial<Gmpq> &qpol,Gmpfr &left,Gmpfr &right,int prec){
        typedef Polynomial<Gmpz>                        ZPolynomial;
        typedef Polynomial_traits_d<ZPolynomial>        ZPtraits;
        typedef ZPtraits::Degree                        ZDegree;
        CGAL_precondition(left<=right);
#ifndef CGAL_NO_PRECONDITIONS
        typedef ZPtraits::Make_square_free              ZSfpart;
        typedef CGAL::RS_AK1::Signat_1<ZPolynomial,Gmpfr>
                                                        Signat;
#endif
        //std::cout<<"refining ["<<left<<","<<right<<"]"<<std::endl;
        Polynomial<Gmpz> zpol=CGAL::RS_AK1::Polynomial_converter_1<
                                        CGAL::Polynomial<Gmpq>,
                                        CGAL::Polynomial<Gmpz> >()(qpol);
#ifndef CGAL_NO_PRECONDITIONS
        ZPolynomial zsfpp=ZSfpart()(zpol);
        Signat signof(zsfpp);
        CGAL::Sign sl=signof(left);
        CGAL_precondition(sl!=signof(right)||(left==right&&sl==ZERO));
#endif
        int deg=ZDegree()(zpol);
        mpz_t* coefficients=(mpz_t*)malloc((deg+1)*sizeof(mpz_t));
        __mpfi_struct interval;
        CGAL_RS_GMPFR_MAKE_UNIQUE(left,temp_left);
        CGAL_RS_GMPFR_MAKE_UNIQUE(right,temp_right);
        interval.left=*(left.fr());
        interval.right=*(right.fr());
        for(int i=0;i<=deg;++i)
                coefficients[i][0]=*(zpol[i].mpz());
        RS2::RS2_calls::init_solver();
        rs3_refine_u_root(&interval,
                          coefficients,
                          deg,
                          prec+CGAL::max(left.get_precision(),
                                         right.get_precision()),
                          1,
                          1);
        free(coefficients);
        mpfr_clear(left.fr());
        mpfr_clear(right.fr());
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

#endif // CGAL_RS_RS3_K_REFINER_1_H
