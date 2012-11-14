// Copyright (c) 2009 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_REFINE_1_RS_H
#define CGAL_RS_REFINE_1_RS_H

#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>
#include <CGAL/RS/algebraic_1.h>
#include <CGAL/assertions.h>

namespace RS3{

inline void refine_1(const CGAL::Algebraic_1 &a,unsigned int s=10000){
        CGAL_precondition(a.inf()<=a.sup());
        // If the algebraic endpoints have room for a refinement bigger
        // than desired, refine as much is possible. This would ensure that
        // small refinements are never made. Other possibility to avoid
        // small refinements is to refine the number just after it was
        // isolated (this early refinement would also have the advantage
        // that one can choose to refine until reaching a certain criterion
        // and assume it is true for later refinements, but this can slow
        // down the Solve_1 functor). And a last option would be to
        // determine a lower bound on the precisions RS likes to refine.
        if(s<mpfr_get_prec(&((mpfi_ptr)(a.mpfi()))->left)-2)
                s=mpfr_get_prec(&((mpfi_ptr)(a.mpfi()))->left)-2;
        if(s<mpfr_get_prec(&((mpfi_ptr)(a.mpfi()))->right)-2)
                s=mpfr_get_prec(&((mpfi_ptr)(a.mpfi()))->right)-2;
        // If the precision of the endpoints is not enough for the desired
        // refinement, allocate a new mpfi and swap later the result.
        if(mpfr_get_prec(&((mpfi_ptr)(a.mpfi()))->left)<s+2||
           mpfr_get_prec(&((mpfi_ptr)(a.mpfi()))->right)<s+2){
                mpfi_t n;
                mpfi_init2(n,s+2);
                mpfi_set(n,a.mpfi());
                rs3_refine_u_root(n,
                                a.pol().get_coefs(),
                                a.pol().get_degree(),
                                s+1,
                                0,
                                0);
                mpfi_swap(n,(mpfi_ptr)a.mpfi());
                // It is not possible to clear the original mpfi because it
                // could be allocated inside RS memory.
                // TODO: clear the mpfi except the first time it is refined
        }else{
                rs3_refine_u_root((mpfi_ptr)a.mpfi(),
                                a.pol().get_coefs(),
                                a.pol().get_degree(),
                                s+1,
                                0,
                                0);
        }
        CGAL_postcondition(a.inf()<=a.sup());
}

} // namespace RS3

#endif  // CGAL_RS_REFINE_1_RS_H
