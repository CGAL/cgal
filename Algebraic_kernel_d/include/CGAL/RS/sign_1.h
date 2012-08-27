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

#ifndef CGAL_RS_SIGN_1_H
#define CGAL_RS_SIGN_1_H

#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <CGAL/RS/basic.h>
#include <CGAL/RS/dyadic.h>
#include <CGAL/RS/polynomial_1.h>
#include <CGAL/RS/algebraic_1.h>
#include <CGAL/RS/solve_1.h>
#include <CGAL/assertions.h>

namespace CGAL{

struct RSSign{

        // This function calculates the sign of the evaluation of a polynomial
        // at a given algebraic number. If it is impossible to know the sign
        // evaluating the interval, it calls sign_1_rs, which uses RS to do it.
        static CGAL::Sign sign_1(const RS_polynomial_1 &p,const Algebraic_1 &x){
                RS::rs_sign s=p.sign_mpfi(x.mpfi());
                if(s!=RS::RS_UNKNOWN)
                        return RS::convert_rs_sign(s);
                return sign_1_rs(p,x);
        }

        // compute the sign of the polynomial at a given dyadic
        static inline CGAL::Sign
        exactsignatdyadic(const RS_polynomial_1 &p,CGALRS_dyadic_srcptr d){
                return p.sign_dyadic(d);
        }

        // compute the sign of the polynomial at a given mpfr
        static inline CGAL::Sign
        exactsignat(const RS_polynomial_1 &p,mpfr_srcptr m){
                return p.sign_mpfr(m);
        }

        // This function uses interval arithmetic to calculate the sign of the
        // evaluation. First, it converts the given mpfr number to an mpfi
        // interval with the given precision. Later, it evaluates the
        // polynomial in the interval an looks for the sign. If it is not able
        // to determine it, it calls the exact signat function.
        static CGAL::Sign quicksignat
        (const RS_polynomial_1 &p,mpfr_srcptr xcoord,mp_prec_t prec){
                mpfi_t x_approx;
                RS::rs_sign s;
                mpfi_init2(x_approx,prec);
                mpfi_set_fr(x_approx,xcoord);
                s=p.sign_mpfi(x_approx);
                mpfi_clear(x_approx);
                if(s!=RS::RS_UNKNOWN)
                        return RS::convert_rs_sign(s);
                return p.sign_mpfr(xcoord);
        }

        // This function is supposed to return the signat calculated in
        // the best way.
        static CGAL::Sign signat(const RS_polynomial_1 &p,mpfr_srcptr xcoord){
                mp_prec_t xprec=mpfr_get_prec(xcoord);
                if(xprec>>(mp_prec_t)12) // switch to exact if xprec>=8192
                        return p.sign_mpfr(xcoord);
                return quicksignat(p,xcoord,xprec);
        }

}; // struct RSSign

} // namespace CGAL

#endif  // CGAL_RS_SIGN_1_H
