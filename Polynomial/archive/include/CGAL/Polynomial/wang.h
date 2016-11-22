// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
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
//
// Author(s)     : Dominik Huelse <dominik.huelse@gmx.de>
//                 Michael Hemmer <hemmer@mpi-inf.mpg.de>
//                 
//
// ============================================================================

/*! \file CGAL/Polynomial/wang.h
 *  \brief Wang's algorithm for Rational Reconstruction. 
 */

#ifndef CGAL_WANG_H
#define CGAL_WANG_H 1


#include <cmath>
#include <CGAL/basic.h>
#include <cstdlib>

namespace CGAL {

namespace internal{

template<typename Integer>
inline
bool wang_general(const Integer& u, const Integer& m, 
        Integer& n, Integer& d,
        const Integer& N, const Integer& D) {
    Integer r0,r1,t0,t1,q,hilf;

//  std::cout<<" wang general "<<std::endl;
    Integer u1 = u;
    if(u1<0){
        u1=u1+m;
    }
    CGAL_precondition(u1>=0);
    CGAL_precondition((m>u) && (2*N*D<m));
    r0=m;
    t0=0;
    r1=u1;
    t1=1;
    while(r1>N){
        q = CGAL::div(r0,r1);
        hilf=r0;
        r0=r1;            
        r1=hilf-q*r1;
        hilf=t0;
        t0=t1;
        t1=hilf-q*t1; 
    }
    n=r1;
    d=t1;
    if(d<0){
        n=-n;
        d=-d; 
    }
 
    if(d<=D && (CGAL::gcd(n,d))==1)
        return true;
    else{
        return false;
    }

} // wang_general


/*! 
 *  \brief Wang's algorithm for Rational Reconstruction
 */
template<typename Integer>
bool wang( const Integer& u, const Integer& m, 
        Integer& n, Integer& d ){
    
    typename CGAL::Algebraic_structure_traits<Integer>::Sqrt sqrt;
    // set N and D to wang's default values
    Integer N = sqrt(CGAL::div(m,Integer(2))); 
    Integer D = N-Integer(1);
    
    return internal::wang_general(u, m, n, d, N, D);  
    
}// wang

} // namespace internal
} // namespace CGAL

#endif // CGAL_WANG_H

// EOF
