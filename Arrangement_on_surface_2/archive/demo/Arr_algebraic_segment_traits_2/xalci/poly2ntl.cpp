// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : NiX
// File          : NiX/poly2ntl.C
// NiX_release   : $Name:  $
// Revision      : $Revision: 1.1 $
// Revision_date : $Date: 2009-06-30 13:14:58 $
//
// Author(s)     : ... ?
//
// ============================================================================

/*! \file NiX/poly2ntl.C
 *  
 *  Conversion of polynomials to NTL format & back
 */


#include <CGAL/config.h>

#define NDEBUG 1
// #define CGAL_POLYNOMIAL_USE_NTL_MUL

// #define CGAL_ACK_BENCHMARK_RES

#ifdef CGAL_POLYNOMIAL_USE_NTL_MUL

//#include <CGAL/leda_integer.h>
#include <CGAL/CORE/BigInt.h>
#include <NTL/ZZX.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Timer.h>

#define POLY_NTL_MIN_DEGREE 10 // minimal degree of polynomials to multiply
                                // using NTL

#endif // CGAL_POLYNOMIAL_USE_NTL_MUL

namespace CGAL {

#ifdef CGAL_POLYNOMIAL_USE_NTL_MUL

// namespace internal {
// int *primes = CGAL::CGALi::primes;
// }  

struct NTL_bigint_rep {
    long alloc;
    long size;
    mp_limb_t data;
};

#endif

#ifdef CGAL_POLYNOMIAL_USE_NTL_MUL

namespace {

// typedef leda_integer CCC_int;
typedef CORE::BigInt Integer;

typedef Polynomial< Integer > Poly_1;

void poly2ntl(const Poly_1& p, NTL::ZZX& q) {

//     if(p.is_zero()) { // TODO should be handled earlier
//         //special handling
//         //a is zero if a.rep.length() == 0;
//         q = NTL::ZZX();
//     } 
    q.rep.SetLength(p.degree() + 1);
    
    int i;
    Poly_1::const_iterator pit;
    for(i = 0, pit = p.begin(); pit != p.end(); pit++, i++) {

        NTL::ZZ& zz = q.rep[i];
        mpz_srcptr tmp = pit->get_mp();
        int sz = tmp->_mp_size;
        if(sz == 0)
            continue;
        if(sz < 0)
            sz = -sz;
        
        zz.SetSize(sz);
        NTL_bigint_rep *rep = (NTL_bigint_rep *)zz.rep;
        rep->size = tmp->_mp_size;
        // copy limbs directly
        memcpy(&rep->data, tmp->_mp_d, sz*sizeof(mp_limb_t));
         //std::cerr << "pit = " << *pit << "; and " <<
            //zz << "\n\n";
    }
}

} // anonymous namespace

#endif

#ifdef CGAL_ACK_BENCHMARK_RES
#warning timing resultants
extern Timer res_tm;
#endif

#ifdef CGAL_POLYNOMIAL_USE_NTL_MUL

#warning using NTL

template <>
Poly_1& Poly_1::operator *= (const Poly_1& p2) {

    Poly_1 p1 = *this;

    if(p1.is_zero() || p2.is_zero()) {
//         std::cout << "mul NTL: zero poly\n";
        return (*this) = Poly_1(Integer(0));
    }
//  TODO: use this if poly size is small..
    if(p1.degree() <= POLY_NTL_MIN_DEGREE &&
        p2.degree() <= POLY_NTL_MIN_DEGREE) {

        internal::Creation_tag TAG;
        Poly_1 p(TAG, p1.degree() + p2.degree() + 1);
        for (int i=0; i <= p1.degree(); ++i)
          for (int j=0; j <= p2.degree(); ++j)
            p.coeff(i+j) += (p1[i]*p2[j]); 
        p.reduce();
//         std::cout << "mul usual: " << p << "\n\n";;
        return (*this) = p ;
    }
    NTL::ZZX q, q2;
    poly2ntl(p1, q);
    poly2ntl(p2, q2);

    q *= q2;

    int d = NTL::deg(q);
//     if(d == -1) {
//         std::cout << "Fatal: zero poly\n";
//         throw -1;
//     }

    this->copy_on_write(); // ??
    this->coeffs().resize(d + 1);

    // TODO: use reduce ??
    mpz_t tmp;
     mpz_init(tmp);
    for(int i = 0; i <= d; i++) {
        
        const NTL::ZZ& zz = q.rep[i];
        if(NTL::IsZero(zz)) {
            coeff(i) = Integer(0);
            continue;
        } 

        NTL_bigint_rep *rep = (NTL_bigint_rep *)zz.rep;
        int sz = rep->size;
        if(sz < 0)
            sz = -sz;
         
        mpz_realloc2(tmp, sz * GMP_NUMB_BITS);
        tmp->_mp_size = rep->size;
        memcpy(tmp->_mp_d, &rep->data, sz*sizeof(mp_limb_t));
         
//         coeff(i).makeCopy();
//         mpz_ptr mpd = coeff(i).get_mp();
//         mpd->_mp_size = rep->size;
//         mpz_realloc2(mpd, sz * GMP_NUMB_BITS);
//         memcpy(mpd->_mp_d, &rep->data, sz*sizeof(mp_limb_t));
        coeff(i) = Integer(tmp);
    
//          mpz_init_set(coeff(i).get_mp(), tmp);
    }
    mpz_clear(tmp);
   
//         CGALi::Creation_tag TAG;
//         Poly_1 p(TAG, p1.degree() + p2.degree() + 1);
//         for (int i=0; i <= p1.degree(); ++i)
//           for (int j=0; j <= p2.degree(); ++j)
//             p.coeff(i+j) += (p1[i]*p2[j]);
//         //p.reduce();
// //         std::cout << "mul usual: " << p << "\n\n";;
//
//     if(*this != p) {
//
//         std::cout << "------------ p1: " << p1 << "---------- p2: " <<
//                 p2 << "\n";
//         std::cout << "FATAL: " << *this << "----------- and " << p << "\n";
//
//     }

//     Poly_1 pp(vec.begin(), vec.end());
    //p.reduce();
//       std::cout << "mul NTL: " << *this << "\n";
    return (*this);// = pp;
}

template <> 
Integer prs_resultant_ufd< Integer >(Poly_1 A, Poly_1 B) {

#ifdef CGAL_ACK_BENCHMARK_RES
// std::cout << "start res " << "\n";
res_tm.start();
#endif

    // implemented using the subresultant algorithm for resultant computation
    // see [Cohen, 1993], algorithm 3.3.7

    typedef Integer NT;

    if (A.is_zero() || B.is_zero()) return NT(0);

    NTL::ZZX q1, q2;
    poly2ntl(A, q1);
    poly2ntl(B, q2);

    NTL::ZZ zz;
    NTL::resultant(zz, q1, q2);

    if(NTL::IsZero(zz))
        return Integer(0);
    
    Integer res;
    NTL_bigint_rep *rep = (NTL_bigint_rep *)zz.rep;
    int sz = rep->size;
    if(sz < 0)
        sz = -sz;
       
    mpz_ptr tmp = res.get_mp();  
    mpz_realloc2(tmp, sz * GMP_NUMB_BITS);
    tmp->_mp_size = rep->size;
    memcpy(tmp->_mp_d, &rep->data, sz*sizeof(mp_limb_t));


//     std::cout << "stop res " << "\n";
#ifdef CGAL_ACK_BENCHMARK_RES
// std::cout << "stop res " << "\n";
res_tm.stop();
#endif

    return res;
}
#else // CGAL_POLYNOMIAL_USE_NTL_MUL

#if 0
template <> 
Integer prs_resultant_ufd< Integer >(Poly_1 A, Poly_1 B) {

#ifdef CGAL_ACK_BENCHMARK_RES
res_tm.start();
#endif

    // implemented using the subresultant algorithm for resultant computation
    // see [Cohen, 1993], algorithm 3.3.7

    typedef Integer NT;

    if (A.is_zero() || B.is_zero()) return NT(0);

    int signflip;
    if (A.degree() <123 B.degree()) {
        Polynomial<NT> T = A; A = B; B = T;
        signflip = (A.degree() & B.degree() & 1);
    } else {
        signflip = 0;
    }

    NT a = A.content(), b = B.content();
    NT g(1), h(1), t = CGAL::ipower(a, B.degree()) * CGAL::ipower(b, A.degree());
    Polynomial<NT> Q, R; NT d;
    int delta;

    A /= a; B /= b;
    do {
        signflip ^= (A.degree() & B.degree() & 1);
        Polynomial<NT>::pseudo_division(A, B, Q, R, d);
        delta = A.degree() - B.degree();
        typedef CGAL::Algebraic_structure_traits<NT>::Is_exact
          Is_exact;
    
        A = B;
        B = R / (g * CGAL::ipower(h, delta));
        g = A.lcoeff();
        // h = h^(1-delta) * g^delta
        CGALi::hgdelta_update(h, g, delta);
    } while (B.degree() > 0);
    // h = h^(1-deg(A)) * lcoeff(B)^deg(A)
    delta = A.degree();
    g = B.lcoeff();
    CGALi::hgdelta_update(h, g, delta);
    h = signflip ? -(t*h) : t*h;
    Algebraic_structure_traits<NT>::Simplify simplify;
    simplify(h);

   return h;
}
#endif

#endif // CGAL_POLYNOMIAL_USE_NTL_MUL

} // namespace CGAL

