// Copyright (c) 2002-2008 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//                 Dominik Huelse  <dominik.huelse@gmx.de>
//
// ============================================================================

/*! \file CGAL/Polynomial/modular_gcd_utils.h
 *  \brief Provides additional utils for the modular GCD calculation
 */

#ifndef CGAL_POLYNOMIAL_MODULAR_GCD_UTILS_H
#define CGAL_POLYNOMIAL_MODULAR_GCD_UTILS_H

#include <CGAL/basic.h>
#include <vector>
#include <CGAL/Polynomial.h>

#include <CGAL/Timer.h>

namespace CGAL{

namespace internal {

template <class NT>
void euclidean_division_obstinate(const NT& F1, const NT& F2,
        NT& Q, NT& R){

    CGAL_precondition(F2 != 0);

    CGAL::div_mod(F1, F2, Q, R);
    CGAL_postcondition(F1 == F2*Q + R);
}


template <class NT>
void euclidean_division_obstinate(const Polynomial<NT>& F1,
        const Polynomial<NT>& F2,
        Polynomial<NT>& Q, Polynomial<NT>& R){

//    std::cout<<" my_modular_gcd_utils "<<std::endl;
    CGAL_precondition(!F2.is_zero());
    int d1 = F1.degree();
    int d2 = F2.degree();
    if ( d1 < d2 ) {
        Q = Polynomial<NT>(NT(0)); R = F1;
        CGAL_postcondition( !(boost::is_same< typename Algebraic_structure_traits<NT>::Is_exact,
                        CGAL::Tag_true >::value) ||  F1 == Q*F2 + R); return;
    }

    typedef std::vector<NT> Vector;
    Vector V_R, V_Q;
    V_Q.reserve(d1);
    if(d2==0){
        for(int i=d1;i>=0;--i){
            V_Q.push_back(F1[i]/F2[0]);
        }
        V_R.push_back(NT(0));
    }
    else{
        V_R.reserve(d1);
        V_R=Vector(F1.begin(),F1.end());
        Vector tmp1;
        tmp1.reserve(d2);
        for(int k=0; k<=d1-d2; ++k){
            V_Q.push_back(V_R[d1-k]/F2[d2]);
            for(int j=0;j<d2;++j){
                tmp1.push_back(F2[j]*V_Q[k]);
            }
            V_R[d1-k]=0;
            for(int i=d1-d2-k;i<=d1-k-1;++i){
                V_R[i]=V_R[i]-tmp1[i-(d1-d2-k)];
            }
            tmp1.clear();
        }


    }
    Q = Polynomial<NT>(V_Q.rbegin(),V_Q.rend());
    R = Polynomial<NT>(V_R.begin(),V_R.end());
    CGAL_postcondition(F1 == F2*Q + R);
}

} // namespace internal
} // namespace CGAL

#endif //#ifnedef CGAL_POLYNOMIAL_MODULAR_GCD_UTILS_H 1

// EOF
