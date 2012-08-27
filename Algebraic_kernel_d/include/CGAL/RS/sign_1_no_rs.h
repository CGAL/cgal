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

#ifndef CGAL_RS_SIGN_1_NO_RS_H
#define CGAL_RS_SIGN_1_NO_RS_H

#include <mpfi.h>
#include <CGAL/RS/basic.h>
#include <CGAL/RS/sign_1_no_rs.h>
#include <CGAL/RS/polynomial_1.h>
#include <CGAL/RS/polynomial_1_utils.h>
#include <CGAL/RS/algebraic_1.h>
#include <CGAL/RS/sign_1.h>
#include <CGAL/RS/refine_1_rs.h>

namespace CGAL{

// compute the sign of the polynomial at a given algebraic number
template <class _Gcd_policy>
CGAL::Sign sign_1_no_rs(const RS_polynomial_1 &p,const Algebraic_1 &x){
        typedef _Gcd_policy     Gcd;

        unsigned bisect_steps=/*4*/1000;
        rs_sign s;
        CGAL::Sign sleft,sright;
        RS_polynomial_1 *gcd,*deriv;
        // we check first by evaluating intervals
        s=p.sign_mpfi(x.mpfi());
        if(s!=RS::RS_UNKNOWN)
                return convert_rs_sign(s);

        // we are not sure, we calculate the gcd in order to know if both
        // polynomials have common roots
        gcd=&(Gcd()(sfpart_1<Gcd>()(p),sfpart_1<Gcd>()(x.pol())));
        if(!gcd->get_degree_static())
                goto refineandreturn;

        // at this point, gcd is not 1; we proceed as follows:
        // -we refine x until having p monotonic in x's interval (to be
        // sure that p has at most one root on that interval)
        // -if the gcd has a root on this interval, both roots are equal
        // (we return 0), otherwise, we refine until having a result

        // how to assure that p is monotonic: when its derivative is never zero
        deriv=&(sfpart_1<Gcd>()(p).derive());
        while(deriv->sign_mpfi(x.mpfi())==RS::RS_UNKNOWN)
                RS3::refine_1(x,(bisect_steps*=2));
                //bisect_n<Gcd>(x,(bisect_steps*=2));

        // how to know that the gcd has a root: just evaluating endpoints
        if((sleft=RSSign::exactsignat(*gcd,x.left()))==CGAL::ZERO||
           (sright=RSSign::exactsignat(*gcd,x.right()))==CGAL::ZERO||
           (sleft!=sright)){
                return CGAL::ZERO;
        }

refineandreturn:
        // the sign is not zero, we refine until having a result
        for(;;){
                RS3::refine_1(x,bisect_steps);
                s=p.sign_mpfi(x.mpfi());
                if(s==RS::RS_POSITIVE)
                        return CGAL::POSITIVE;
                if(s==RS::RS_NEGATIVE)
                        return CGAL::NEGATIVE;
                bisect_steps*=2;
        }
}

} // namespace CGAL

#endif  // CGAL_RS_SIGN_1_NO_RS_H
