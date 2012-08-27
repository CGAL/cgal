// Copyright (c) 2009-2010 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_SIGN_1_RS_H
#define CGAL_RS_SIGN_1_RS_H

#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <CGAL/RS/algebraic_1.h>
#include <CGAL/assertions.h>

namespace RS3{

inline CGAL::Sign sign_1(const CGAL::RS_polynomial_1 &p,
                         const CGAL::Algebraic_1 &a){

        mpz_t* list_constr[]={p.get_coefs()};
        int list_degs[]={p.get_degree()};
        mpfi_t* sols_u=(mpfi_t*)malloc(sizeof(mpfi_t));
        *(sols_u[0])=*(a.mpfi());
        mpfi_t **vals_constr=(mpfi_t**)malloc(sizeof(mpfi_t*));
        vals_constr[0]=(mpfi_t*)malloc(sizeof(mpfi_t));
        mpfi_init(vals_constr[0][0]);

        rs3_refine_eval_u(a.pol().get_coefs(),a.pol().get_degree(),NULL,0,
                          (const mpz_t **)list_constr,list_degs,
                          1,sols_u,1,
                          vals_constr,NULL,MPFR_PREC_MIN,1,1,1);

        if(mpfr_zero_p(&vals_constr[0][0]->left)!=0 &&
           mpfr_zero_p(&vals_constr[0][0]->right)!=0){
                return CGAL::ZERO;
        }
        CGAL_assertion(mpfi_has_zero(vals_constr[0][0])<=0);
        if(mpfi_is_pos(vals_constr[0][0])){
                return CGAL::POSITIVE;
        }
        return CGAL::NEGATIVE;
}

} // namespace RS3

#endif  // CGAL_RS_SIGN_1_RS_H
