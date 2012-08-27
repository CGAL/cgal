// Copyright (c) 2011 National and Kapodistrian University of Athens (Greece).
// All rights reserved.
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

#ifndef CGAL_RS_ISOLE_1_H
#define CGAL_RS_ISOLE_1_H

#include <CGAL/RS/basic.h>
#include <CGAL/RS/rs_calls_1.h>
#include <CGAL/Gmpfi.h>
#include <vector>

namespace CGAL{

namespace RS{

// CGAL::RS::isolator<P> solves a polynomial of type P and returns a vector
// of Gmpfi's containing the solutions. Currently, the only implemented
// isolator solves polynomials of type CGAL::Polynomial<CGAL::Gmpz>.

template <class Polynomial_>
struct isolator{
        inline std::vector<Gmpfi>
        operator()(const Polynomial_&,unsigned int prec=CGAL_RS_DEF_PREC);
};

template <class Polynomial_>
inline std::vector<Gmpfi>
isolator<Polynomial_>::operator()(const Polynomial_ &p,unsigned int prec){
        CGAL_error_msg(
                "isolator not implemented for this type of polynomials");
        return std::vector<Gmpfi>();
}

template <>
inline std::vector<Gmpfi>
isolator<Polynomial<Gmpz> >::operator()(const Polynomial<Gmpz> &p,
                                        unsigned int prec){
        int numsols;
        unsigned int degree=p.degree();
        mpz_t *coeffs=(mpz_t*)malloc((degree+1)*sizeof(mpz_t));
        mpfi_ptr *intervals_mpfi=(mpfi_ptr*)malloc(degree*sizeof(mpfi_ptr));
        std::vector<Gmpfi> intervals;
        for(unsigned int i=0;i<=degree;++i)
                coeffs[i][0]=*(p[i].mpz());
        init_solver();
        create_rs_upoly(coeffs,degree,rs_get_default_up());
        free(coeffs);
        set_rs_precisol(prec);
        set_rs_verbose(CGAL_RS_VERB);
        rs_run_algo(CGALRS_CSTR("UISOLE"));
        numsols=affiche_sols_eqs(intervals_mpfi);
        for(int j=0;j<numsols;++j)
                intervals.push_back(Gmpfi(intervals_mpfi[j]));
        free(intervals_mpfi);
        return intervals;
}

} // namespace RS

} // namespace CGAL

#endif  // CGAL_RS_ISOLE_1_H
