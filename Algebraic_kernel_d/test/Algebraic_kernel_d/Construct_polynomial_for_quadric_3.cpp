// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
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
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>

#include <CGAL/config.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Algebraic_kernel_d/Construct_polynomial_for_quadric_3.h>

template < class NT >
void test_Construct_polynomial_for_quadric_3() {

    typedef CGAL::Polynomial< NT >       Poly_nt_1;
    typedef CGAL::Polynomial< Poly_nt_1 > Poly_nt_2;
    typedef CGAL::Polynomial< Poly_nt_2 > Poly_nt_3;
    
    NT pa, pb, pc, pd, pe, pf, pg, ph, pk, pl;
    
    // q_p := -85*x^2-55*x*y-37*x*z-35*x+97*y^2+50*y*z+79*y+56*z^2+49*z+63;
    pa = -85;
    pb = -55;
    pc = -37;
    pd = 97;
    pe = 50;
    pf = 56;
    pg = -35;
    ph = 79;
    pk = 49;
    pl = 63;
    
    Poly_nt_1 p00(pl, pg, pa);
    Poly_nt_1 p01(ph, pb);
    Poly_nt_1 p02(pd);

    Poly_nt_1 p10(pk,pc);
    Poly_nt_1 p11(pe);    

    Poly_nt_1 p20(pf);

    Poly_nt_2 p0(p00, p01, p02);
    Poly_nt_2 p1(p10, p11);
    Poly_nt_2 p2(p20);

    Poly_nt_3 p(p0,p1,p2);
    
    assert(CGAL::Construct_polynomial_for_quadric_3< NT >()
           (pa,pb,pc,pd,pe,pf,pg,ph,pk,pl) == p);
    
    typedef CGAL::Exponent_vector EV;
    typename CGAL::Polynomial_traits_d< Poly_nt_3 >::Get_innermost_coefficient
        icoefficient;

    // second test
    assert(icoefficient(p,EV(2,0,0)) == -85);
    assert(icoefficient(p,EV(1,1,0)) == -55);
    assert(icoefficient(p,EV(1,0,1)) == -37);
    assert(icoefficient(p,EV(0,2,0)) ==  97);
    assert(icoefficient(p,EV(0,1,1)) ==  50);
    assert(icoefficient(p,EV(0,0,2)) ==  56);
    assert(icoefficient(p,EV(1,0,0)) == -35);
    assert(icoefficient(p,EV(0,1,0)) ==  79);
    assert(icoefficient(p,EV(0,0,1)) ==  49);
    assert(icoefficient(p,EV(0,0,0)) ==  63);
}

int main() {
#if CGAL_USE_LEDA
    {
        typedef CGAL::LEDA_arithmetic_kernel AK;
        test_Construct_polynomial_for_quadric_3< AK::Integer >();
    }
#endif
#if LiS_HAVE_CORE
    {
        typedef CGAL::CORE_arithmetic_kernel AK;
        test_Construct_polynomial_for_quadric_3< AK::Integer >();
    }
#endif
    return 0;
}



