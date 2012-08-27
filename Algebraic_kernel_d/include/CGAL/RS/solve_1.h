// Copyright (c) 2006-2009 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_SOLVE_1_H
#define CGAL_RS_SOLVE_1_H

#include <CGAL/basic.h>
#include <CGAL/RS/basic.h>
#include <CGAL/RS/dyadic.h>
#include <CGAL/RS/polynomial_1.h>
#include <CGAL/RS/algebraic_1.h>
#include <CGAL/RS/rs_calls_1.h>
#include <CGAL/Gmpfi.h>
#include <vector>

namespace CGAL{

class RS_polynomial_1;
class Algebraic_1;

// solve given the precision, returns the number of roots
inline int solve_1(mpfi_ptr *x,
                   const RS_polynomial_1 &p1,
                   unsigned int prec=CGAL_RS_DEF_PREC){
        rs_reset_all();
        create_rs_upoly(p1.get_coefs(),p1.get_degree(),rs_get_default_up());
        set_rs_precisol(prec);
        set_rs_verbose(CGAL_RS_VERB);
        rs_run_algo(CGALRS_CSTR("UISOLE"));
        return affiche_sols_eqs(x);
}

// calculate the sign of a polynomial evaluated at the root of another
inline Sign sign_1_rs(const RS_polynomial_1 &p1,
                      const Algebraic_1 &a,
                      unsigned int prec=CGAL_RS_MIN_PREC){
        mpz_t **constr;
        int *degs;
        CGAL_assertion(a.is_consistent());
        rs_reset_all ();
        // tell RS to find the roots of this polynomial
        create_rs_upoly (a.pol().get_coefs (), a.pol().get_degree (),
                        rs_get_default_up ());
        // the constraint vector will have just one element
        constr = (mpz_t**)malloc(sizeof(mpz_t*));
        *constr = p1.get_coefs ();
        degs = (int*)malloc(sizeof(int));
        *degs = p1.get_degree ();
        create_rs_uconstr (constr, degs, rs_get_default_ineqs_u ());
        set_rs_precisol (prec);
        set_rs_verbose (CGAL_RS_VERB);
        rs_run_algo(CGALRS_CSTR("UISOLES"));
        return affiche_signs_constr(a.nr());
}

} // namespace CGAL

#endif  // CGAL_RS_SOLVE_1_H
