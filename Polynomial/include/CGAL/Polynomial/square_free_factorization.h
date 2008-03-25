// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================
#ifndef CGAL_POLYNOMIAL_SQUARE_FREE_FACTORIZATION_H
#define CGAL_POLYNOMIAL_SQUARE_FREE_FACTORIZATION_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>

CGAL_BEGIN_NAMESPACE
;

namespace POLYNOMIAL {

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

template <class NT, class OutputIterator1, class OutputIterator2> inline
int square_free_factorization_for_regular_polynomial_( Polynomial<NT> a, OutputIterator1 factors,
        OutputIterator2 multiplicities, Field_tag ) {
    
    typedef typename Fraction_traits< Polynomial<NT> >::Is_fraction
        Is_fraction;

    return square_free_factorization_for_regular_polynomial_(a, factors, multiplicities, Is_fraction());
}

template <class NT, class OutputIterator1, class OutputIterator2> inline
int square_free_factorization_for_regular_polynomial_(Polynomial<NT> a, OutputIterator1 factors,
        OutputIterator2 multiplicities, Unique_factorization_domain_tag ) {
            
    return square_free_factorization_for_regular_polynomial_(a, factors, multiplicities);
}

template <class NT, class OutputIterator1, class OutputIterator2> inline
int square_free_factorization_for_regular_polynomial_(Polynomial<NT> a, OutputIterator1 factors,
        OutputIterator2 multiplicities, ::CGAL::Tag_false ) {
    
    return square_free_factorization_for_regular_polynomial_(a, factors, multiplicities);
}

template <class NT, class OutputIterator1, class OutputIterator2> inline
int square_free_factorization_for_regular_polynomial_(Polynomial<NT> a, OutputIterator1 factors,
        OutputIterator2 multiplicities, ::CGAL::Tag_true ) {
            
    typedef Polynomial<NT> POLY;
    typedef typename Fraction_traits<POLY>::Numerator_type INTPOLY;
    typedef typename Fraction_traits<POLY>::Denominator_type DENOM;
    typename Fraction_traits<POLY>::Decompose decompose;
    typename Fraction_traits<POLY>::Compose compose;
    
    typedef typename INTPOLY::NT INTNT;
    typedef std::vector<INTPOLY> PVEC;
    typedef typename PVEC::iterator Iterator;
    typedef typename Algebraic_structure_traits<INTNT>::Algebraic_category Algebraic_category;

    DENOM dummy;
    PVEC fac;
    std::back_insert_iterator<PVEC> fac_bi(fac);

    a.simplify_coefficients();
    INTPOLY p;
    decompose(a,p, dummy);
    int d = square_free_factorization_utcf_(p, fac_bi, multiplicities, Algebraic_category());
    for (Iterator it = fac.begin(); it != fac.end(); ++it) {
        POLY q = compose(*it, DENOM(1));
        q /= q.lcoeff();
        q.simplify_coefficients();
        *factors++ = q;
    }
    return d;
}

template <class NT,  class OutputIterator1, class OutputIterator2>
int square_free_factorization_for_regular_polynomial_(Polynomial<NT> a, OutputIterator1 factors,
        OutputIterator2 multiplicities ) {
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
    
    typedef Polynomial<NT> POLY;
    if (a.degree() == 0) return 0;
    POLY b = diff(a);  
    POLY c = gcd(a, b);
   
    if (c == NT(1)) {
        *factors = a;
        *multiplicities = 1;
        return 1;
    }

    int i = 1, n = 0;
    POLY w = a/c, y = b/c, z = y - diff(w), g;
    while (!z.is_zero()) {
        g = gcd(w, z);
        if (g.degree() > 0) {
            *factors++ = g;
            *multiplicities++ = i;
            n++;
        }
        i++;
        w /= g;
        y = z/g;
        z = y - diff(w);
    }
    *factors = w;
    *multiplicities++ = i;
    n++;

    return n;
}


/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
 *  \brief factor the univariate polynomial \c p by multiplicities.
 *
 *  That means: factor it into square-free and pairwise coprime non-constant
 *  factors <I>g<SUB>i</SUB></I> with multiplicities <I>m<SUB>i</SUB></I>
 *  such that <I>p = alpha * g<SUB>1</SUB><SUP>m<SUB>1</SUB></SUP> *
 *  ... * g<SUB>n</SUB><SUP>m<SUB>n</SUB></SUP> </I>.
 *  This is known as square-free factorization in the literature.
 *  The number \e n is returned. The factors \e gi and multiplicities
 *  \e mi are written through the respective output iterators.
 *
 *  \pre The coefficient domain \c NT must be a \c field or \c UFDomain
 *  of characteristic 0.
 *  \c OutputIterator1 must allow the value type \c Polynomial<NT>.
 *  \c OutputIterator2 must allow the value type \c int.
 *
 *  Use this function if you are sure, that the polynomial has multiple 
 *  factors, otherwise use NiX::filtered_square_free_factorization.
 */

template< class ICoeff, class OutputIterator1, class OutputIterator2 > 
int square_free_factorization( ICoeff, OutputIterator1, OutputIterator2 ) {
    return 0;
}

template< class Coeff, class OutputIterator1, class OutputIterator2 > 
int square_free_factorization( Polynomial<Coeff> polynomial, 
                               OutputIterator1 factors,
                               OutputIterator2 multiplicities ) {
    typedef typename Polynomial_traits_d< Polynomial< Coeff > >::Univariate_content Univariate_content;
    typedef typename Algebraic_structure_traits<Coeff>::Algebraic_category Algebraic_category;
    
    Coeff univariate_content = Univariate_content()( polynomial );    

    if( typename Polynomial_traits_d< Coeff >::Total_degree()(univariate_content) > 0 ) {
        Polynomial< Coeff > regular_polynomial = polynomial / univariate_content;        
        int multiplicity = square_free_factorization_for_regular_polynomial_( polynomial / univariate_content,
                                                          factors, multiplicities, Algebraic_category() );

        typedef std::vector< Coeff > Factors_uc;
        typedef std::vector< int > Multiplicities_uc;
        Factors_uc factors_uc;
        Multiplicities_uc multiplicities_uc;
        multiplicity += square_free_factorization( univariate_content, 
                                        std::back_inserter(factors_uc),
                                        std::back_inserter(multiplicities_uc) );
        
        for( typename Factors_uc::iterator it = factors_uc.begin(); it != factors_uc.end(); ++it )
            *factors++ = Polynomial<Coeff>((*it) / CGAL::unit_part((*it)) );
         
        for( Multiplicities_uc::iterator it = multiplicities_uc.begin();
             it != multiplicities_uc.end(); ++it )
            *multiplicities++ = (*it);
        
        return multiplicity; 
    } else {
        return square_free_factorization_for_regular_polynomial_( polynomial,
                                 factors, multiplicities, Algebraic_category() );    
    }                                 
} 

template <class NT, class OutputIterator1, class OutputIterator2> inline
int square_free_factorization_utcf_(
                                    Polynomial<NT> a,
                                    OutputIterator1 factors,
                                    OutputIterator2 multiplicities,
                                    Field_tag)
{
    typedef typename Fraction_traits< Polynomial<NT> >::Is_fraction
        Is_fraction;

    return square_free_factorization_for_regular_polynomial_(a, factors, multiplicities, Is_fraction());
}

template <class NT, class OutputIterator1, class OutputIterator2> inline
int square_free_factorization_utcf_(
                                    Polynomial<NT> a,
                                    OutputIterator1 factors,
                                    OutputIterator2 multiplicities,
                                    Unique_factorization_domain_tag)
{
    return square_free_factorization_for_regular_polynomial_(a, factors, multiplicities);
}

template <class NT,  class OutputIterator1, class OutputIterator2>
int square_free_factorization_utcf_(
                                Polynomial<NT> a, 
                                OutputIterator1 factors,
                                OutputIterator2 multiplicities,
                                Integral_domain_tag)
{
    // Yun's Square-Free Factorization
    // see [Geddes et al, 1992], Algorithm 8.2

    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits_d<POLY>::Innermost_leading_coefficient ilcoeff;
    //typename Polynomial_traits_d<POLY>::Innermost_coefficient_to_polynomial ictp;
    typename Polynomial_traits_d<POLY>::Innermost_coefficient_begin begin;
    typename Polynomial_traits_d<POLY>::Innermost_coefficient_end end;
    typename Algebraic_extension_traits<IC>::Denominator_for_algebraic_integers dfai;
    typename Algebraic_extension_traits<IC>::Normalization_factor nfac;
    typename Scalar_factor_traits<POLY>::Scalar_factor sfac;  
    typename Scalar_factor_traits<POLY>::Scalar_div sdiv;
    typedef typename Scalar_factor_traits<POLY>::Scalar Scalar;

    if (a.degree() == 0) return 0;

    a = canonicalize_polynomial(a);
    POLY b = diff(a);
    POLY c = gcd_utcf(a, b);

    if (c == NT(1)) {
        *factors = a;
        *multiplicities = 1;
        return 1;
    }

    int i = 1, n = 0;

    // extending both polynomials a and b by the denominator for algebraic 
    // integers, which comes out from c=gcd(a,b), such that a and b are 
    // divisible by c
    IC lcoeff = ilcoeff(c);
    IC denom = dfai(begin(c), end(c));
    lcoeff *= denom * nfac(denom);
    POLY w = (a * POLY(lcoeff)) / c;
    POLY y = (b * POLY(lcoeff)) / c;

    // extracting a common scalar factor out of w=a/c and y=b/c simultaneously,
    // such that the parts in z=y-w' are canceled out as they should
    Scalar sfactor = sfac(y,sfac(w));
    sdiv(w, sfactor); 
    sdiv(y, sfactor);

    POLY  z = y - diff(w);
    POLY g;

    while (!z.is_zero()) {
        g = gcd_utcf(w, z);
        if (g.degree() > 0) {
            *factors++ = g;
            *multiplicities++ = i;
            n++;
        }
        i++;
        lcoeff = ilcoeff(g); // same as above
        denom =dfai(begin(c), end(c)); 
        lcoeff *= denom * nfac(denom);
        w = (w * POLY(lcoeff)) / g;
        y = (z * POLY(lcoeff)) / g;
        Scalar sfactor = sfac(y,sfac(w));
        sdiv(w, sfactor); 
        sdiv(y, sfactor);
       
        z = y - diff(w);
    }
    *factors = w;
    *multiplicities++ = i;
    n++;

    return n;
}

/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
 *  \brief factor the univariate polynomial \c p by multiplicities with respect
 *  to constant factors
 *
 *  Same functionality as square_free_factorization, but
 *    a) no prefactor \c alpha is returned due to the respect to constant factors,
 *    b) the coefficient domain \c NT may also be \c IntegralDomain now.
 *
 */
template <class NT,  class OutputIterator1, class OutputIterator2> inline
int square_free_factorization_utcf(
                                   const Polynomial<NT>& p, 
                                   OutputIterator1 factors,
                                   OutputIterator2 multiplicities)
{
    typedef typename Algebraic_structure_traits<NT>::Algebraic_category Algebraic_category;

    NT alpha = typename Polynomial_traits_d< Polynomial< NT > >::Univariate_content_up_to_constant_factor()(p);
    return square_free_factorization_utcf_
        (p/alpha, factors, multiplicities, Algebraic_category());
}





} // namespace POLYNOMIAL

CGAL_END_NAMESPACE

#endif // CGAL_POLYNOMIAL_SQUARE_FREE_FACTORIZATION_H
