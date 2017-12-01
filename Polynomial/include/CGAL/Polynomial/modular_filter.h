// Copyright (c) 2002-2008 Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     : Michael Hemmer 
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_POLYNOMIAL_MODULAR_FILTER_H
#define CGAL_POLYNOMIAL_MODULAR_FILTER_H

#include <CGAL/basic.h>

#include <CGAL/Polynomial.h>
#include <CGAL/polynomial_utils.h>
#include <CGAL/Polynomial/prs_resultant.h>
#include <CGAL/Modular_traits.h>

namespace CGAL {

namespace internal {
    template <class NT> inline
    bool
    may_have_common_factor_(
        const Polynomial<NT>& p1,
        const Polynomial<NT>& p2,
        ::CGAL::Tag_true){
      
      // Enforce IEEE double precision and to nearest before using modular arithmetic
      CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);
      
        CGAL_precondition(p1.degree()!=-1);
        CGAL_precondition(p2.degree()!=-1);
        
        if(CGAL::total_degree(p1) == 0){ return p1.is_zero();}
        if(CGAL::total_degree(p2) == 0){ return p2.is_zero();}
        
        typedef Polynomial<NT>                Polynomial_nt;
        typedef Modular_traits<Polynomial_nt> MT;
        typedef typename MT::Residue_type       Polynomial_mt;

        typename MT::Modular_image    modular_image;

        typename CGAL::Polynomial_traits_d<Polynomial_nt>::Degree_vector
            degree_vector_nt;
            
        typename CGAL::Polynomial_traits_d<Polynomial_mt>::Degree_vector
            degree_vector_mt;

        Polynomial_mt m1=modular_image(p1);
        Polynomial_mt m2=modular_image(p2);

        CGAL::Exponent_vector exp_vec_nt_1 = degree_vector_nt( p1 );
        CGAL::Exponent_vector exp_vec_nt_2 = degree_vector_nt( p2 );
        CGAL::Exponent_vector exp_vec_mt_1 = degree_vector_mt( m1 );
        CGAL::Exponent_vector exp_vec_mt_2 = degree_vector_mt( m2 );
                
        // return true if the exponent vector changes (degree loss)
        if( (exp_vec_nt_1 != exp_vec_mt_1) ||
            (exp_vec_nt_2 != exp_vec_mt_2 ) )
            return true;

        typename CGAL::Polynomial_traits_d<Polynomial_mt>::Total_degree 
            tdegree_mt;

        if( tdegree_mt( CGAL::gcd( m1, m2 ) ) > 0 )
            return true;
        else
            return false;                
    }
   
    template <class NT> inline
    bool may_have_common_factor_(const Polynomial<NT>& ,
                                 const Polynomial<NT>& ,
                                 ::CGAL::Tag_false) {return true;}
    
/*! \ingroup CGAL_polynomial_utils
 *  \brief Test whether \c P and \c Q may have a common factor. 
 *
 *  This function is based on a fast modular arithmetic and serves as a
 *  filter to avoid expensive exact computations to determine whether \c P
 *  and \c Q have a common factor.. 
 *  
 *  If the function return false, then \c P and \c Q have no common factor for
 *  sure. 
 *  \c P and \c Q have with a high probability a common factor, if it returns
 *  true. 
 */
template <class NT> inline 
bool may_have_common_factor(const Polynomial<NT>& P,
                            const Polynomial<NT>& Q){
// TODO: Should this compiler switch be renamed?
#ifdef CGAL_MODULAR_FILTER_OFF
    return true;
#endif

    CGAL_precondition( Residue::get_current_prime()!=0 );
    typedef Polynomial<NT> POLY;
    typedef Modular_traits<POLY> Mtr;
    typename Mtr::Is_modularizable is_modularizable;
    return internal::may_have_common_factor_(P,Q,is_modularizable);   
}

/*! \ingroup CGAL_polynomial_utils
 *  \brief Test whether the polynomial \c P may has a multiple root.
 *
 *  This function is based on a fast modular arithmetic and serves as a
 *  filter to avoid expensive exact computations to determine whether \c P
 *  has a multiple root or not. 
 *
 *  If the function return false, then \c P has no multiple root for sure.
 *  \c P has with a high probability a multiple root, if it returns true.
 */
template <class NT> inline
bool may_have_multiple_factor_(const Polynomial<NT>& P, CGAL::Tag_true ){

  // Enforce IEEE double precision and to nearest before using modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);

    // Create modular images of p
    typedef Polynomial<NT>                Polynomial_nt;
    typedef Modular_traits<Polynomial_nt> MT;
    typedef typename MT::Residue_type       Polynomial_mt;

    typename MT::Modular_image    modular_image;

    typename CGAL::Polynomial_traits_d<Polynomial_nt>::Degree_vector
        degree_vector_nt;
            
    typename CGAL::Polynomial_traits_d<Polynomial_mt>::Degree_vector
        degree_vector_mt;
            
    Polynomial_mt m = modular_image( P );
    
    CGAL::Exponent_vector exp_vec_nt = degree_vector_nt( P );
    CGAL::Exponent_vector exp_vec_mt = degree_vector_mt( m );
    
    if( exp_vec_nt != exp_vec_mt )
        return true;
    
    // Check modular image to be square free
    typename CGAL::Polynomial_traits_d< Polynomial_mt >::Is_square_free 
        is_square_free;

    return( !is_square_free( m ) ); 
}

template< class NT > inline
bool may_have_multiple_factor_( const Polynomial<NT>&, CGAL::Tag_false ) {
    return true;
}

template< class NT > inline
bool may_have_multiple_factor( const Polynomial<NT>& P ) {
  if(CGAL::total_degree(P) <= 1)
        return false;

    // Modular filter
    CGAL_precondition( Residue::get_current_prime()!=0 );
    typedef Polynomial<NT> POLY;
    typedef Modular_traits<POLY> Mtr;
    typename Mtr::Is_modularizable is_modularizable;
    return internal::may_have_multiple_factor_(P, is_modularizable);       
}

} //namespace internal
} //namespace CGAL

#endif //CGAL_MODULAR_FILTER_H
