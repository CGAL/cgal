// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// Author(s)     :  Michael Hemmer
//
// ============================================================================
#ifndef CGAL_POLYNOMIAL_SQUARE_FREE_FACTORIZATION_H
#define CGAL_POLYNOMIAL_SQUARE_FREE_FACTORIZATION_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial/fwd.h>
#include <CGAL/Polynomial/misc.h>
#include <CGAL/Polynomial/Polynomial_type.h>

namespace CGAL {
namespace internal {

// square-free factorization
// 
// the implementation uses two dispatches:
//   a) first look at the coefficient's algebra type
//   b) if the algebra type is of the concept field, try to decompose
// the same holds for square-free factorization up to constant factors (utcf)
//
// sqff -------> algebra type ? ----field-----> decomposable ?
//                     |                   A       |     |
//                    UFD                  |       no   yes
//                     |                   |       |     |
//                     V                           |     |
//               Yun's algo <----------------------      |
//                     A                                 |
//                     |                   |             |
//                    UFD                  |             V
//                     |                   |         decompose and use
// sqff_utcf --> algebra_type ? ----field--          sqff_utcf with numerator
//                     |
//               integral domain
//                     |
//                     V
//               Yun's algo (utcf)

template <class IC,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize(const IC&, OutputIterator1, OutputIterator2){
    return 0;
}
template <class IC,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial(const IC&, OutputIterator1, OutputIterator2){
    return 0;
}

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2);
template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2);
template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial_(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2, CGAL::Tag_true);
template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial_(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2, CGAL::Tag_false);
template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial_(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2, Integral_domain_tag);
template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial_(const Polynomial<Coeff>&, OutputIterator1, OutputIterator2, Unique_factorization_domain_tag);




template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize
(const Polynomial<Coeff>&  poly, OutputIterator1 factors, OutputIterator2 multiplicities)
{
    typedef Polynomial<Coeff> POLY;
    typedef Polynomial_traits_d< POLY > PT;
    typedef typename PT::Construct_polynomial Construct_polynomial;
    typedef typename PT::Univariate_content_up_to_constant_factor Ucont_utcf;
    typedef typename PT::Integral_division_up_to_constant_factor  Idiv_utcf;
    
    if (typename PT::Total_degree()(poly) == 0){return 0;}
   
    Coeff ucont_utcf = Ucont_utcf()(poly); 
    POLY  regular_poly = Idiv_utcf()(poly,Construct_polynomial()(ucont_utcf));

    int result = square_free_factorize_for_regular_polynomial( 
            regular_poly, factors, multiplicities);

    if (CGAL::total_degree(ucont_utcf) > 0){
        typedef std::vector< Coeff > Factors_uc;
        typedef std::vector< int > Multiplicities_uc;
        Factors_uc factors_uc;
        Multiplicities_uc multiplicities_uc;
        result += square_free_factorize( ucont_utcf,
                std::back_inserter(factors_uc),
                std::back_inserter(multiplicities_uc) );
        
        for( typename Factors_uc::iterator it = factors_uc.begin(); 
             it != factors_uc.end(); ++it ){
          *factors++ = Construct_polynomial()(*it);
        }
        for( Multiplicities_uc::iterator it = multiplicities_uc.begin();
             it != multiplicities_uc.end(); ++it ){
            *multiplicities++ = (*it);
        }
    }
    return result;                                
}

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial
(const Polynomial<Coeff>&  p, OutputIterator1 factors, OutputIterator2 multiplicities){
    typedef Polynomial<Coeff> POLY;
    typedef typename CGAL::Fraction_traits<POLY>::Is_fraction Is_fraction;
    return square_free_factorize_for_regular_polynomial_(p,factors,multiplicities,Is_fraction());
}

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial_
(const Polynomial<Coeff>& p, OutputIterator1 factors, OutputIterator2 multiplicities, CGAL::Tag_true){
    
    typedef Polynomial<Coeff> POLY;
    typedef Polynomial_traits_d< POLY > PT;
    typedef Fraction_traits<POLY> FT; 
    
    typename FT::Numerator_type num; 
    typename FT::Denominator_type denom; 
    typename FT::Decompose()(p,num,denom);
    
    std::vector<typename FT::Numerator_type> ifacs;
    int result =  square_free_factorize_for_regular_polynomial(num,std::back_inserter(ifacs),multiplicities);
    
    typedef typename std::vector<typename FT::Numerator_type>::iterator Iterator;
    denom = typename FT::Denominator_type(1);
    for ( Iterator it = ifacs.begin(); it != ifacs.end(); ++it) {
        POLY q = typename FT::Compose()(*it, denom);
        *factors++ = typename PT::Canonicalize()(q);
    }
    
    return result; 
}

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial_
(const Polynomial<Coeff>& p, OutputIterator1 factors, OutputIterator2 multiplicities, CGAL::Tag_false){
    typedef Polynomial<Coeff> POLY;
    typedef typename Algebraic_structure_traits<POLY>::Algebraic_category Algebraic_category;
    return square_free_factorize_for_regular_polynomial_(p,factors,multiplicities,Algebraic_category());
}

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial_
(const Polynomial<Coeff>& p, OutputIterator1 factors, OutputIterator2 multiplicities, Integral_domain_tag){
    // Yun's Square-Free Factorization
    // see [Geddes et al, 1992], Algorithm 8.2

    typedef Polynomial<Coeff> POLY;
    typedef Polynomial_traits_d<POLY> PT;
    typedef typename PT::Innermost_coefficient_type IC;
    typename PT::Innermost_leading_coefficient ilcoeff;
    //typename PT::Innermost_coefficient_to_polynomial ictp;
    typename PT::Construct_innermost_coefficient_const_iterator_range range;
    typename Algebraic_extension_traits<IC>::Denominator_for_algebraic_integers dfai;
    typename Algebraic_extension_traits<IC>::Normalization_factor nfac;
    typename Scalar_factor_traits<POLY>::Scalar_factor sfac;  
    typename Scalar_factor_traits<POLY>::Scalar_div sdiv;
    typedef typename Scalar_factor_traits<POLY>::Scalar Scalar;

    
    if (typename PT::Total_degree()(p) == 0) return 0;

    POLY a = CGAL::canonicalize(p);
    POLY b = CGAL::differentiate(a);
    POLY c = CGAL::internal::gcd_utcf_(a, b);

    if (c == Coeff(1)) {
        *factors = a;
        *multiplicities = 1;
        return 1;
    }

    int i = 1, n = 0;

    // extending both polynomials a and b by the denominator for algebraic 
    // integers, which comes out from c=gcd(a,b), such that a and b are 
    // divisible by c
    IC lcoeff = ilcoeff(c);
    IC denom = dfai(range(c).first, range(c).second);
    lcoeff *= denom * nfac(denom);
    POLY w = (a * POLY(lcoeff)) / c;
    POLY y = (b * POLY(lcoeff)) / c;

    // extracting a common scalar factor out of w=a/c and y=b/c simultaneously,
    // such that the parts in z=y-w' are canceled out as they should
    Scalar sfactor = sfac(y,sfac(w));
    sdiv(w, sfactor); 
    sdiv(y, sfactor);

    POLY  z = y - CGAL::differentiate(w);
    POLY g;

    while (!z.is_zero()) {
        g = CGAL::internal::gcd_utcf_(w, z);
        if (g.degree() > 0) {
            *factors++ = g;
            *multiplicities++ = i;
            n++;
        }
        i++;
        lcoeff = ilcoeff(g); // same as above
        denom =dfai(range(c).first, range(c).second); 
        lcoeff *= denom * nfac(denom);
        w = (w * POLY(lcoeff)) / g;
        y = (z * POLY(lcoeff)) / g;
        Scalar sfactor = sfac(y,sfac(w));
        sdiv(w, sfactor); 
        sdiv(y, sfactor);
       
        z = y - CGAL::differentiate(w);
    }

    *factors = w;
    *multiplicities++ = i;
    n++;

    return n;
}

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_for_regular_polynomial_
(const Polynomial<Coeff>& p, OutputIterator1 factors, OutputIterator2 multiplicities, Unique_factorization_domain_tag){
    // Yun's Square-Free Factorization
    // see [Geddes et al, 1992], Algorithm 8.2
    /* 
       @inproceedings{y-osfda-76,
       author = {David Y.Y. Yun},
       title = {On square-free decomposition algorithms},
       booktitle = {SYMSAC '76: Proceedings of the third ACM symposium on Symbolic 
                    and algebraic computation},
       year = {1976},
       pages = {26--35},
       location = {Yorktown Heights, New York, United States},
       doi = {http://doi.acm.org/10.1145/800205.806320},
       publisher = {ACM Press},
       address = {New York, NY, USA},
       }
    */
    
    typedef Polynomial<Coeff> POLY;
    typedef Polynomial_traits_d<POLY> PT;

    if (typename PT::Total_degree()(p) == 0) return 0;
    POLY a = CGAL::canonicalize(p);

    POLY b = CGAL::differentiate(a);  
    POLY c = CGAL::gcd(a, b);
   
    if (c == Coeff(1)) {
        *factors = a;
        *multiplicities = 1;
        return 1;
    }

    int i = 1, n = 0;
    POLY w = a/c, y = b/c, z = y - CGAL::differentiate(w), g;
    while (!z.is_zero()) {
        g = CGAL::gcd(w, z);
        if (g.degree() > 0) {
            *factors++ = g;
            *multiplicities++ = i;
            n++;
        }
        i++;
        w /= g;
        y = z/g;
        z = y - CGAL::differentiate(w);
    }
    *factors = w;
    *multiplicities++ = i;
    n++;

    return n;
}




// square-free factorization utcf 

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_utcf
(const Polynomial<Coeff>&  p, OutputIterator1 factors, OutputIterator2 multiplicities)
{
    return square_free_factorize(p,factors,multiplicities);
}

template <class Coeff,  class OutputIterator1, class OutputIterator2>
inline int square_free_factorize_utcf_for_regular_polynomial
(const Polynomial<Coeff>&  p, OutputIterator1 factors, OutputIterator2 multiplicities)
{
    return square_free_factorize_for_regular_polynomial(p,factors,multiplicities);
}




// filtered square-free factorization

// ### filtered versions #### 

/*! \brief As NiX::square_free_factorize, but filtered by  
 *  NiX::may_have_multiple_root 
 *  
 *  Use this function if the polynomial might be square free. 
 */  
template <class Coeff, class OutputIterator1, class OutputIterator2> 
inline
int filtered_square_free_factorize(
                                       Polynomial<Coeff> p,
                                       OutputIterator1 factors,
                                       OutputIterator2 multiplicities)
{
  if(CGAL::internal::may_have_multiple_factor(p)){
      return internal::square_free_factorize(p, factors, multiplicities);
  }else{
      *factors++        = CGAL::canonicalize(p);
      *multiplicities++ = 1;
      return 1;
  }
}

/*! \brief As NiX::square_free_factorize_utcf, but filtered by  
 *  NiX::may_have_multiple_root 
 *  
 *  Use this function if the polynomial might be square free. 
 */  
template <class Coeff,  class OutputIterator1, class OutputIterator2> 
inline
int filtered_square_free_factorize_utcf( const Polynomial<Coeff>& p, 
                                         OutputIterator1 factors,
                                         OutputIterator2 multiplicities)
{
    return filtered_square_free_factorize(p,factors,multiplicities);
}

} // namespace internal
} //namespace CGAL

#endif // CGAL_POLYNOMIAL_SQUARE_FREE_FACTORIZATION_H
