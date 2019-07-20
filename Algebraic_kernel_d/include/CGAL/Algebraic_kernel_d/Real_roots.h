// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.


/*! \file NiX/Real_roots.h
    \brief This file defines the class NiX::Real_roots. 
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_REAL_ROOTS_H
#define CGAL_ALGEBRAIC_KERNEL_D_REAL_ROOTS_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>

#include <list>
#include <vector>
#include <queue>


namespace CGAL { 

namespace internal {
    
/*! \ingroup NiX_Real_roots
 *  \brief This class provides operators for a comfortable construction of
 *  AlgebraicReal from Polynomials using a specific RealRootIsolator.  
 *
 *  A valid template argument for AlgebraicReal is NiX::Algebraic_real. 
 */
template < class AlgebraicReal , 
           class RealRootIsolator  >
class Real_roots{
public:
    //! The Real_roots type it self.
    typedef Real_roots<AlgebraicReal,RealRootIsolator> Self;
    //! First template argument.
    typedef AlgebraicReal Algebraic_real;
    //! Second template argument. 
    typedef RealRootIsolator  Real_root_isolator;
    //! The Polnomial type used by Algebraic_real and the Real_root_isolator.
    typedef typename Algebraic_real::Polynomial_1 Polynomial;

private:
    typedef typename AlgebraicReal::Coefficient Coefficient;    
    typedef typename AlgebraicReal::Rational    Rational; 
private: 
template < class PolynomialIterator,
           class IntIterator,
           class AlgebraicRealOutputIterator,              
           class IntOutputIterator>
int gen_agebraic_reals_with_mults( PolynomialIterator           fac,
                                   PolynomialIterator           fac_end,
                                   IntIterator                  mul,
                                   IntIterator                  CGAL_precondition_code(mul_end),
                                   AlgebraicRealOutputIterator  oi_root,
                                   IntOutputIterator            oi_mult){   
    
    Self real_roots;
    typedef std::pair<Algebraic_real, int> PAIR;
    
    // find zeroes of each factor and sort them in ascending order
    std::priority_queue< PAIR,std::vector<PAIR>,std::greater<PAIR> > pqueue;
    std::vector<Algebraic_real>                    tmp;

    while(fac !=  fac_end){
        CGAL_assertion(mul != mul_end);
        tmp.clear();
        real_roots(*fac, std::back_inserter(tmp));
        for (int j = 0; j < static_cast<int>(tmp.size()); j++) {
            pqueue.push(PAIR(tmp[j], *mul));
        }
        fac++;
        mul++;
    }
    
    // output factors and multiplicities
    int n = 0;
    while (!pqueue.empty()) {
        *oi_root++ = pqueue.top().first;
        *oi_mult++ = pqueue.top().second;
        n++;
        pqueue.pop();
    }
    return n;
}

private:
template < class PolynomialConstIterator,
           class IntConstIterator,
           class PolynomialOutputIterator>
void write_factors_by_multiplicity(PolynomialConstIterator  fac,
                                   PolynomialConstIterator  fac_end,
                                   IntConstIterator         mul,
                                   IntConstIterator         CGAL_assertion_code(mul_end),
                                   PolynomialOutputIterator oi_poly){
    // output table such that table[m] contains square-free factor of
    // multiplicity m (or else constant poly 1)
    int m = 0;
    while (fac != fac_end) {
        CGAL_assertion(mul != mul_end);
        while (m < *mul) {
            *oi_poly++ = Polynomial(Coefficient(1)); m++;
        }
        *oi_poly++ = *fac; m++;
        ++fac; ++mul;
    }
}

public:
/*! \brief computes all roots of the square free polynomial P in
 * ascending order and returns the number of real roots.
 */   
template <class AlgebraicRealOutputIterator>
int operator()(const Polynomial&        poly , 
               AlgebraicRealOutputIterator it){        
    
    CGAL_precondition_msg( typename CGAL::Polynomial_traits_d< Polynomial >::Is_square_free()(poly), "P not square free.");
    return (*this)(Real_root_isolator(poly),it);
}
   

public:
/*! \brief computes all roots of the polynomial P in ascending order
 * and their multiplicity and returns the number of real roots.
 *
 *  This operator returns the number \e n of distinct real roots of \c poly.
 *  Each root is represented as an object of an instance AlgebraicReal. 
 *  The operator writes these \e n real zeroes in ascending order to \c
 *  oi_root. 
 *  It writes the multiplicities of the zeroes in the same order to
 *  \c oi_mult .   
 */
template <class AlgebraicRealOutputIterator, 
          class IntOutputIterator>
int operator()(const Polynomial&           poly,
               AlgebraicRealOutputIterator oi_root,
               IntOutputIterator           oi_mult){  
    CGAL_precondition(CGAL::degree(poly) >= 0);
    // fast exit 
    if (CGAL::degree(poly) == 0) 
        return (poly.is_zero())?-1:0;
          
    std::list<Polynomial>       sqffac;
    std::list<int>              facmul;
    
    filtered_square_free_factorize_utcf(poly,
					    std::back_inserter(sqffac), 
					    std::back_inserter(facmul));
    
    
    int number_of_real_roots=
        gen_agebraic_reals_with_mults(sqffac.begin(),sqffac.end(),
                                      facmul.begin(),facmul.end(),
                                      oi_root,
                                      oi_mult);
    return number_of_real_roots;
}
         
public:     
/*! \brief computes all roots defined by the Real_root_isolator object in
 * ascending order 
 */
template <class AlgebraicRealOutputIterator>
int operator()(const Real_root_isolator&    isolator, 
               AlgebraicRealOutputIterator  it){ 
    
    Polynomial poly = isolator.polynomial();
    //cout << "P: "<<poly<<endl;
    // take out exact known roots out of Poly
    for(int j = 0 ; j < isolator.number_of_real_roots(); j++){
        if(isolator.is_exact_root(j)){
            Rational root(isolator.left_bound(j));
            CGAL::simplify(root);
            typedef CGAL::Fraction_traits<Rational> FT;
            typename FT::Numerator_type num;
            typename FT::Denominator_type denom;
            typename FT::Decompose decomp;
            decomp(root,num,denom);
            Polynomial linear_factor(Coefficient(-num),
                                     Coefficient(denom));
            poly=CGAL::integral_division(poly,linear_factor);                
        }
    }   
    //cout << "P_without_exact: "<<poly<<endl;        
    std::vector<AlgebraicReal> conjugated_roots;
    std::back_insert_iterator<std::vector<AlgebraicReal> > con_it 
        = std::back_inserter(conjugated_roots);
    // construct AlgebraicReal
    for(int j = 0 ; j < isolator.number_of_real_roots(); j++){
        if(isolator.is_exact_root(j)){
            // exact roots (Rational
            Rational root=isolator.left_bound(j);
            CGAL::simplify(root);
            *it++=AlgebraicReal(root);
        }else{
            // other roots 
            Rational left = isolator.left_bound(j);
            Rational right= isolator.right_bound(j);
            CGAL::simplify(left);
            CGAL::simplify(right);
            AlgebraicReal tmp(poly,left,right);
            *it++=tmp;
            *con_it++=tmp;
        } 
    }   
    AlgebraicReal::conjugate(conjugated_roots.begin(),
                             conjugated_roots.end());
    return isolator.number_of_real_roots();
}

/*! \brief factor \c p by multiplicities, return both factors and their roots
 *
 *  This operator is for those users who have an
 *  Polynomial \c poly which is not necessarily square-free, and who
 *  want to get both its square-free factorization and its real roots with
 *  their respective multiplicities.
 *
 *  This operator returns the number \e n of distinct real roots of \c poly .
 *  Each root is represented as an object of an instance AlgebraicReal. 
 *  The operator writes these \e n real zeroes in ascending order to \c
 *  oi_root. 
 *  It writes the multiplicities of the roots in the same order to
 *  \c oi_mult .                                                              
 *  Finally, it writes the square-free factors of \c p
 *  to \c oi_poly such that the factor #<I>k</I> written is the factor
 *  with exponent \e k in the square-free factorization of \c p , or
 *  the constant polynomial 1 if a factor of multiplicity \e k does
 *  not occur. Yes, this means that the first factor written, which is
 *  factor #0, will always be 1.
 *
 *  The data types involved are determined the \c AlgebraicReal
 *  template argument. In particular, \c p must be of type
 *  \c NiX::Polynomial<AlgebraicReal::Coefficient>, as are its factors; 
 *  the roots are of type \c AlgebraicReal ; and the multiplicities are 
 *  of type \c int .
 *
 */
template <class AlgebraicRealOutputIterator,
          class IntOutputIterator, 
          class PolynomialOutputIterator>
int operator()( const Polynomial&           poly, 
                AlgebraicRealOutputIterator oi_root,
                IntOutputIterator           oi_mult, 
                PolynomialOutputIterator    oi_poly) {
    
    CGAL_precondition(CGAL::degree(poly) >= 0);
    
    // fast exit 
    if (CGAL::degree(poly) == 0) 
        return (poly.is_zero())?-1:0;
          
    
    std::list<Polynomial>       sqffac;
    std::list<int>              facmul;
    
    filtered_square_free_factorize_utcf(poly,
					    std::back_inserter(sqffac),
					    std::back_inserter(facmul));
    
    write_factors_by_multiplicity(sqffac.begin(),sqffac.end(),
                                  facmul.begin(),facmul.end(),
                                  oi_poly);
    
    int numer_of_real_roots =
      gen_agebraic_reals_with_mults(sqffac.begin(),sqffac.end(),
				    facmul.begin(),facmul.end(),
				    oi_root,
				    oi_mult);
    return numer_of_real_roots;
	  }
	   };

} // namespace internal

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_REAL_ROOTS_ROOTS_H
