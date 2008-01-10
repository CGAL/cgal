// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Tobias Reithmann <treith@mpi-inf.mpg.de>
//                 Michael Hemmer   <hemmer@informatik.uni-mainz.de>
//
// ============================================================================

/*! \file NiX/polynomial_gcd.h
 *   \brief Greatest common divisors and related operations on polynomials. 
 */

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Fraction_traits.h>

#include <CGAL/Polynomial/ipower.h>
#include <CGAL/Polynomial/hgdelta_update.h>
#include <CGAL/Polynomial/square_free_factorization.h>
#include <CGAL/Polynomial/modular_filter.h>

#ifndef CGAL_POLYNOMIAL_POLYNOMIAL_GCD_H
#define CGAL_POLYNOMIAL_POLYNOMIAL_GCD_H

#ifndef CGAL_POLYNOMIAL_GCD_AVOID_CANONICALIZE
#define CGAL_POLYNOMIAL_GCD_AVOID_CANONICALIZE 1
#endif 



CGAL_BEGIN_NAMESPACE

// Forward declarations of ALL functions of this file to avoid problems with
//  some compilers.
//template< class NT > inline bool may_have_multiple_root( const Polynomial<NT>& P );

namespace POLYNOMIAL {

template <class NT> inline 
Polynomial<NT> gcd_( const Polynomial<NT>& p1, const Polynomial<NT>& p2 );
template <class NT> inline
Polynomial<NT> gcd_innermost_coefficient_dispatch( const Polynomial<NT>& p1, 
                                         const Polynomial<NT>& p2, Field_tag );
template <class NT> inline
Polynomial<NT> gcd_innermost_coefficient_dispatch( const Polynomial<NT>& p1, 
               const Polynomial<NT>& p2, Unique_factorization_domain_tag tag );

template <class NT> inline
Polynomial<NT> gcd_( const Polynomial<NT>& p1, const Polynomial<NT>& p2, 
                                                           ::CGAL::Tag_false ); 
template <class NT> inline
Polynomial<NT> gcd_( const Polynomial<NT>& p1, const Polynomial<NT>& p2, 
                                                        ::CGAL::Tag_true tag );
#ifndef NiX_POLY_USE_PRIMITIVE_GCD
    template <class NT>
    Polynomial<NT> gcd_( Polynomial<NT> p1, Polynomial<NT> p2, 
                                             Unique_factorization_domain_tag ); 
#else
    template <class NT>
    Polynomial<NT> gcd_( Polynomial<NT> p1, Polynomial<NT> p2, UFDomain_tag );
#endif // NiX_POLY_USE_PRIMITIVE_GCD

template <class NT>
Polynomial<NT> gcd_( Polynomial<NT> p1, Polynomial<NT> p2, Field_tag );
template <class NT> inline
Polynomial<NT> gcd( const Polynomial<NT>& p1, const Polynomial<NT>& p2 );
template <class NT> inline
NT gcd_utcf_( const NT& a, const NT& b );
template <class NT> inline
Polynomial<NT> gcd_utcf_( const Polynomial<NT>& p1, const Polynomial<NT>& p2 );
template <class NT> inline
Polynomial<NT> gcd_utcf_( const Polynomial<NT>& p1, const Polynomial<NT>& p2, 
                                                           ::CGAL::Tag_false );
template <class NT>
Polynomial<NT> gcd_utcf_( Polynomial<NT> p1, Polynomial<NT> p2, 
                                                            ::CGAL::Tag_true );
template <class NT> inline
Polynomial<NT> gcd_utcf_( const Polynomial<NT>& p1, const Polynomial<NT>& p2, 
                                                               Field_tag tag );
template <class NT>
NT content_utcf_( const Polynomial<NT>& p );
template <class NT>
Polynomial<NT> gcd_utcf_( Polynomial<NT> p1, Polynomial<NT> p2, 
                                                         Integral_domain_tag );

} // namespace POLYNOMIAL

template <class NT> inline
Polynomial<NT> gcd_utcf( const Polynomial<NT>& p1, const Polynomial<NT>& p2 );

namespace INTERN_POLYNOMIAL_GCD {

template <class NT> inline
Polynomial<NT> gcdex_( Polynomial<NT> x, Polynomial<NT> y, Polynomial<NT>& xf, 
                                       Polynomial<NT>& yf, ::CGAL::Tag_false );
template <class NT>
Polynomial<NT> gcdex_( Polynomial<NT> x, Polynomial<NT> y, Polynomial<NT>& xf, 
                                               Polynomial<NT>& yf, Field_tag );
template <class NT> inline
Polynomial<NT> gcdex_( Polynomial<NT> x, Polynomial<NT> y, Polynomial<NT>& xf, 
                                        Polynomial<NT>& yf, ::CGAL::Tag_true );

} // namespace POLYNOMIAL

template <class NT> inline
Polynomial<NT> gcdex( Polynomial<NT> p1, Polynomial<NT> p2, Polynomial<NT>& f1, 
                                                          Polynomial<NT>& f2 );
template <class NT>
Polynomial<NT> pseudo_gcdex(
#ifdef DOXYGEN_RUNNING
    Polynomial<NT> p1, Polynomial<NT> p2,
    Polynomial<NT>& f2, Polynomial<NT>& f2, NT& v
#else
    Polynomial<NT> x, Polynomial<NT> y,
    Polynomial<NT>& xf, Polynomial<NT>& yf, NT& vf
#endif // DOXYGEN_RUNNING
);

template <class NT, class OutputIterator1, class OutputIterator2> inline
int filtered_square_free_factorization( Polynomial<NT> p, 
                     OutputIterator1 factors, OutputIterator2 multiplicities );
template <class NT,  class OutputIterator1, class OutputIterator2> inline
int filtered_square_free_factorization_utcf( const Polynomial<NT>& p, 
                     OutputIterator1 factors, OutputIterator2 multiplicities );


// end of prototypes
////////////////////////////////////////////////////////////////////////////////


// TODO: This is a dummy-version of may_have_multiple_root. The original
//       EXACUS function is defined in polynomial_utils.h, but does currently
//       not work because of the missing modular traits (and other stuff).
/*template <class NT> inline
bool may_have_multiple_root(const Polynomial<NT>& P){
  return true;
}*/


// 1) gcd (basic form without cofactors)
//     uses three level of dispatch on tag types:
//     a) if the algebra type of the innermost coefficient is a field,
//         ask for decomposability. for UFDs compute the gcd directly
//     b) if NT supports integralization, the gcd is computed on
//         integralized polynomials
//     c) over a field (unless integralized), use the Euclidean algorithm;
//         over a UFD, use the subresultant algorithm

namespace POLYNOMIAL {

template <class NT> inline
Polynomial<NT> gcd_(
    const Polynomial<NT>& p1, const Polynomial<NT>& p2
) {
    typedef typename Polynomial_traits_d<NT>::Innermost_coefficient Innermost_coefficient;
    typedef typename Algebraic_structure_traits<Innermost_coefficient>::Algebraic_category Algebraic_category;

    // dispatch 1a) innermost coefficients algebra type is field or UFD
    return gcd_innermost_coefficient_dispatch(p1, p2, Algebraic_category());
}

template <class NT> inline
Polynomial<NT> gcd_innermost_coefficient_dispatch(
    const Polynomial<NT>& p1, const Polynomial<NT>& p2, Field_tag
) {
    typedef typename Fraction_traits< Polynomial<NT> >
        ::Is_fraction Is_fraction;

    // dispatch 1b) NT is decomposable or not
    return gcd_(p1, p2, Is_fraction());
}

template <class NT> inline
Polynomial<NT> gcd_innermost_coefficient_dispatch(
    const Polynomial<NT>& p1, const Polynomial<NT>& p2, Unique_factorization_domain_tag tag
) {
    return gcd_(p1, p2, tag);
}

template <class NT> inline
Polynomial<NT> gcd_(
    const Polynomial<NT>& p1, const Polynomial<NT>& p2, ::CGAL::Tag_false
) {
    typedef typename Algebraic_structure_traits<NT>::Algebraic_category Algebraic_category;

    // dispatch 1c) choose the right algorithm for a field or a UFD
    return gcd_(p1, p2, Algebraic_category());
}

template <class NT> inline
Polynomial<NT> gcd_(
    const Polynomial<NT>& p1, const Polynomial<NT>& p2, ::CGAL::Tag_true tag
) {
    return gcd_utcf_(p1, p2, tag);
}

#ifndef NiX_POLY_USE_PRIMITIVE_GCD
template <class NT>
Polynomial<NT> gcd_(
    Polynomial<NT> p1, Polynomial<NT> p2, Unique_factorization_domain_tag
) {
    // implemented using the subresultant algorithm for gcd computation
    // see [Cohen, 1993], algorithm 3.3.1

    // handle trivial cases
    if (p1.is_zero()){
        if (p2.is_zero()) return Polynomial<NT>(NT(1));
        else return p2 / p2.unit_part();
    }
    if (p2.is_zero())
        return p1 / p1.unit_part();
    if (p2.degree() > p1.degree()) {
        Polynomial<NT> p3 = p1; p1 = p2; p2 = p3;
    }

    // compute gcd of content
    NT p1c = p1.content(), p2c = p2.content();
    NT gcdcont = gcd(p1c,p2c);

    // compute gcd of primitive parts
    p1 /= p1c; p2 /= p2c;

    NT dummy;
    Polynomial<NT> q, r;

    NT g = NT(1), h = NT(1);
    for (;;) { 
        Polynomial<NT>::pseudo_division(p1, p2, q, r, dummy);
        if (r.is_zero()) { break; }
        if (r.degree() == 0) { return Polynomial<NT>(gcdcont); }
        int delta = p1.degree() - p2.degree();
        p1 = p2;
        p2 = r / (g * POLYNOMIAL::ipower(h, delta));
        g = p1.lcoeff();
        // h = h^(1-delta) * g^delta
        POLYNOMIAL::hgdelta_update(h, g, delta);
    }
    

    p2 /= p2.content() * p2.unit_part();

    // combine both parts to proper gcd
    p2 *= gcdcont; 
    p2.simplify_coefficients();
    return p2;
}

#else

template <class NT>
Polynomial<NT> gcd_(
    Polynomial<NT> p1, Polynomial<NT> p2, UFDomain_tag
) {
    // implemented by computing the primitive remainder sequence

    // handle trivial cases
    if (p1.is_zero())
        if (p2.is_zero()) return Polynomial<NT>(NT(1));
        else return p2 / p2.unit_part();
    if (p2.is_zero())
        return p1 / p1.unit_part();
    if (p2.degree() > p1.degree()) {
        Polynomial<NT> p3 = p1; p1 = p2; p2 = p3;
    }

    NT p1c = p1.content(), p2c = p2.content();
    p1 /= p1c; p2 /= p2c;
    NT F = gcd(p1c,p2c);
    Polynomial<NT> q,r;
    NT D;
    while ( ! p2.is_zero() ) { 
        Polynomial<NT>::pseudo_division(p1,p2,q,r,D);
        r /= r.content();
        p1=p2; p2=r;
    }
    p1 /= p1.unit_part();
    return F * p1;
}

#endif // NiX_POLY_USE_PRIMITIVE_GCD

template <class NT>
Polynomial<NT> gcd_(
    Polynomial<NT> p1, Polynomial<NT> p2, Field_tag
) {
    // handle trivial cases
    if (p1.is_zero()){
        if (p2.is_zero()) return Polynomial<NT>(NT(1));
        else return p2 / p2.unit_part();
    }
    if (p2.is_zero())
        return p1 / p1.unit_part();
    if (p2.degree() > p1.degree()) {
        Polynomial<NT> p3 = p1; p1 = p2; p2 = p3;
    }

    Polynomial<NT> q, r;
    while (!p2.is_zero()) { 
        Polynomial<NT>::euclidean_division(p1, p2, q, r);
        p1 = p2; p2 = r;
    }
    p1 /= p1.lcoeff();
    p1.simplify_coefficients();
    return p1;
}

} // namespace POLYNOMIAL

// name gcd() forwarded to the Intern::gcd_() dispatch function
/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
 *  \brief return the greatest common divisor of \c p1 and \c p2
 *
 *  \pre Requires \c NT to be a \c Field or a \c UFDomain.
 *
 *  Internally, computation is performed ``denominator-free'' if
 *  supported by the coefficient type via \c NiX::Fraction_traits.
 *  The gcd is computed from a polynomial remainder sequence (euclidean
 *  or subresultant.) By defining \c NiX_POLY_USE_PRIMITIVE_GCD,
 *  the subresultant PRS can be replaced globally by the primitive PRS.
 *
 *  Mathematically, a gcd in NT[x] is defined only up to multiplication
 *  with a unit (invertible element) of the coefficient ring NT. This
 *  function returns a unit-normal gcd <I>d</I>, i.e.
 *  <I>d</I><TT>.unit_part() == NT(1)</TT>.
 *  If \c NT is a \c Field , all non-zero scalars are units and
 *  unit-normality means that \e d is normalized to have
 *  <I>d</I><TT>.lcoeff() == NT(1)</TT>.
 *  If \c NT is a \c UFDomain, then non-invertible scalar factors 
 *  <B>do</B> matter, and unit-normality typically means something
 *  like <I>d</I><TT>.lcoeff()</TT> being positive.
 */
namespace POLYNOMIAL {
template <class NT> inline
Polynomial<NT> gcd(const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return POLYNOMIAL::gcd_(p1,p2); }
} // namespace POLYNOMIAL

// 2) gcd_utcf computation
//     (gcd up to scalar factors, for non-UFD non-field coefficients)
//     a) first try to decompose the coefficients
//     b) second dispatch depends on the algebra type of NT

namespace POLYNOMIAL {

template <class NT> inline
NT gcd_utcf_(const NT&, const NT&)
{ return NT(1); }

template <class NT> inline
Polynomial<NT> gcd_utcf_(
    const Polynomial<NT>& p1, const Polynomial<NT>& p2
) {
    typedef typename Fraction_traits< Polynomial<NT> >
        ::Is_fraction Is_fraction;

    // dispatch 2a) NT is decomposable or not
    return gcd_utcf_(p1, p2, Is_fraction());
}

template <class NT> inline
Polynomial<NT> gcd_utcf_(
    const Polynomial<NT>& p1, const Polynomial<NT>& p2, ::CGAL::Tag_false
) {
    typedef typename Algebraic_structure_traits<NT>::Algebraic_category Algebraic_category;

    // dispatch 2b) choose the right algorithm for field and ring
    return gcd_utcf_(p1, p2, Algebraic_category());
}

template <class NT>
Polynomial<NT> gcd_utcf_(
    Polynomial<NT> p1, Polynomial<NT> p2, ::CGAL::Tag_true
) {
    typedef Polynomial<NT> POLY;
    typedef typename Fraction_traits<POLY>::Numerator_type INTPOLY;
    typedef typename Fraction_traits<POLY>::Denominator_type DENOM;
    typedef typename INTPOLY::NT INTNT;

    DENOM dummy;
    p1.simplify_coefficients();
    p2.simplify_coefficients();
    INTPOLY p1i = integralize_polynomial(p1, dummy);
    INTPOLY p2i = integralize_polynomial(p2, dummy);

    typedef typename Algebraic_structure_traits<INTNT>::Algebraic_category Algebraic_category;
    INTPOLY d0 = gcd_utcf_(p1i, p2i, Algebraic_category());
    POLY d = fractionalize_polynomial<POLY>(d0, DENOM(1));
    d /= d.unit_part();
    d.simplify_coefficients();
    return d;
}

template <class NT> inline
Polynomial<NT> gcd_utcf_(
    const Polynomial<NT>& p1, const Polynomial<NT>& p2, Field_tag tag
) {
    return gcd_(p1, p2, tag);
}

template <class NT>
NT content_utcf_(const Polynomial<NT>& p)
{
    typename Algebraic_structure_traits<NT>::Integral_division idiv;
    typename Algebraic_structure_traits<NT>::Unit_part upart;
    typedef typename Polynomial<NT>::const_iterator const_iterator;

    const_iterator it = p.begin(), ite = p.end();
    while (*it == NT(0)) it++;
    NT cont = idiv(*it, upart(*it));
    for( ; it != ite; it++) {
        if (cont == NT(1)) break;
        if (*it != NT(0)) cont = gcd_utcf_(cont, *it);
    }
    return cont;
}

template <class NT>
Polynomial<NT> gcd_utcf_(
    Polynomial<NT> p1, Polynomial<NT> p2, Integral_domain_tag
) {
    
    // handle trivial cases
    if (p1.is_zero()){
        if (p2.is_zero()){
            return Polynomial<NT>(NT(1));
        }else{
#if NiX_POLYNOMIAL_GCD_AVOID_CANONICALIZE
            return p2;
#else
            return canonicalize_polynomial(p2);
#endif
        }
    }
    if (p2.is_zero()){
#if NiX_POLYNOMIAL_GCD_AVOID_CANONICALIZE
        return p1;
#else
        return canonicalize_polynomial(p1);
#endif
    }
    if (p2.degree() > p1.degree()) {
        Polynomial<NT> p3 = p1; p1 = p2; p2 = p3;
    }

    // remove redundant scalar factors
    p1=canonicalize_polynomial(p1);
    p2=canonicalize_polynomial(p2); 

    // compute content of p1 and p2
    NT p1c = content_utcf_(p1);
    NT p2c = content_utcf_(p2);

    // compute gcd of content
    NT gcdcont = gcd_utcf_(p1c, p2c);

    // compute gcd of primitive parts
    p1 = div_utcf(p1, p1c, true); 
    p2 = div_utcf(p2, p2c, true); 

 
    Polynomial<NT> q, r;
    
    // TODO measure preformance of both methodes with respect to 
    // univariat polynomials on Integeres
    // univariat polynomials on Sqrt_extension<Integer,Integer>
    // multivariat polynomials
    // May write specializations for different cases 
#if 0
    // implemented using the subresultant algorithm for gcd computation
    // with respect to constant scalar factors
    // see [Cohen, 1993], algorithm 3.3.1
    NT g = NT(1), h = NT(1), dummy;
    for (;;) { 
        Polynomial<NT>::pseudo_division(p1, p2, q, r, dummy);
        if (r.is_zero()) { break; }
        if (r.degree() == 0) { return Polynomial<NT>(gcdcont); }
        int delta = p1.degree() - p2.degree();
        p1 = p2;
        p2 = r / (g * POLYNOMIAL::ipower(h, delta));
        g = p1.lcoeff();
        // h = h^(1-delta) * g^delta
        Intern::hgdelta_update(h, g, delta);
    }
#else
    // implentaion using just the 'naive' methode
    // but performed much better as the one by Cohen
    // (for univariat polynomials with Sqrt_extension coeffs )
    NT dummy;
    for (;;) { 

        Polynomial<NT>::pseudo_division(p1, p2, q, r, dummy);    

        if (r.is_zero()) { break; }
        if (r.degree() == 0) { return Polynomial<NT>(gcdcont); }
        p1 = p2;
        p2 = r ;
        
        p2=canonicalize_polynomial(p2);   
    }
#endif

    Polynomial<NT> cutcf = content_utcf_(p2);
    p2 = div_utcf(p2, cutcf, true);

    // combine both parts to proper gcd
    p2 *= gcdcont;
    
#if NiX_POLYNOMIAL_GCD_AVOID_CANONICALIZE
    return p2;
#else
    // make poly unique
    return canonicalize_polynomial(p2);
#endif
 
}

} // namespace POLYNOMIAL

// name gcd_utcf() forwarded to the Intern::gcd_utcf_() dispatch function
/*! \relates NiX::Polynomial
 *  \brief return a constant multiple of gcd(p1,p2)
 *
 *  Over a non-UFD or non-field NT, the polynomial ring NT[x]
 *  may not possess greatest common divisors. However, concept
 *  \c IntegralDomain requires NT to be an integral domain, such that we
 *  can consider its quotient field Q(NT) over which gcds
 *  of polynomials exist. \c gcd_utcf(p1,p2) is a constant
 *  denominator-free multiple of gcd(p1,p2) in Q(NT)[x].
 *  It may not be a divisor of \c p1 and \c p2 in NT[x].
 *
 *  The result is unit-normal.
 */
template <class NT> inline
Polynomial<NT> gcd_utcf(const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return POLYNOMIAL::gcd_utcf_(p1,p2); }


// 3) extended gcd computation (with cofactors)
//     with dispatch similar to gcd

namespace POLYNOMIAL {

template <class NT> inline
Polynomial<NT> gcdex_(
    Polynomial<NT> x, Polynomial<NT> y,
    Polynomial<NT>& xf, Polynomial<NT>& yf,
    ::CGAL::Tag_false
) {
    typedef typename Algebraic_structure_traits<NT>::Algebraic_category Algebraic_category;
    return gcdex_(x, y, xf, yf, Algebraic_category());
};

template <class NT>
Polynomial<NT> gcdex_(
    Polynomial<NT> x, Polynomial<NT> y,
    Polynomial<NT>& xf, Polynomial<NT>& yf,
    Field_tag
) {
    /* The extended Euclidean algorithm for univariate polynomials.
     * See [Cohen, 1993], algorithm 3.2.2
     */
    typedef Polynomial<NT> POLY;
    typename Algebraic_structure_traits<NT>::Integral_division idiv;

    // handle trivial cases
    if (x.is_zero()) {
        if (y.is_zero()) CGAL_error_msg("gcdex(0,0) is undefined");
        xf = NT(0); yf = idiv(NT(1), y.unit_part());
        return yf * y;
    }
    if (y.is_zero()) {
        yf = NT(0); xf = idiv(NT(1), x.unit_part());
        return xf * x;
    }
    bool swapped = x.degree() < y.degree();
    if (swapped) { POLY t = x; x = y; y = t; }

    // main loop
    POLY u = x, v = y, q, r, m11(1), m21(0), m21old;
    for (;;) {
        /* invariant: (i) There exist m12 and m22 such that
         *   u = m11*x + m12*y
         *   v = m21*x + m22*y
         * (ii) and we have
         *   gcd(u,v) == gcd(x,y)
         */

        // compute next element of remainder sequence
        POLY::euclidean_division(u, v, q, r);  //  u == qv + r
        if (r.is_zero()) break;

        // update u and v while preserving invariant
        u = v; v = r;
        /* Since r = u - qv, this preserves invariant (part ii)
         * and corresponds to the matrix assignment
         *   (u) = (0  1) (u)
         *   (v)   (1 -q) (v)
         */
        m21old = m21; m21 = m11 - q*m21; m11 = m21old;
        /* This simulates the matching matrix assignment
         *   (m11 m12) = (0  1) (m11 m12)
         *   (m21 m22)   (1 -q) (m21 m22)
         * which preserves the invariant (part i)
         */
        if (r.degree() == 0) break;
    }
    /* postcondition:  invariant holds  and  v divides u */

    // make gcd unit-normal
    m21 /= v.unit_part(); v /= v.unit_part();

    // obtain m22 such that  v == m21*x + m22*y
    POLY m22;
    POLY::euclidean_division(v - m21*x, y, m22, r);
    CGAL_assertion(r.is_zero());

    // check computation
    CGAL_assertion(v == m21*x + m22*y);

    // return results
    if (swapped) {
        xf = m22; yf = m21;
    } else {
        xf = m21; yf = m22;
    }
    return v;
}

template <class NT> inline
Polynomial<NT> gcdex_(
    Polynomial<NT> x, Polynomial<NT> y,
    Polynomial<NT>& xf, Polynomial<NT>& yf,
    ::CGAL::Tag_true
) {
    typedef Polynomial<NT> POLY;
    typedef typename Fraction_traits<POLY>::Numerator_type INTPOLY;
    typedef typename Fraction_traits<POLY>::Denominator_type DENOM;
    typedef typename INTPOLY::NT INTNT;

    // rewrite  x as xi/xd  and  y as yi/yd  with integral polynomials xi, yi
    DENOM xd, yd;
    x.simplify_coefficients();
    y.simplify_coefficients();
    INTPOLY xi = integralize_polynomial(x, xd);
    INTPOLY yi = integralize_polynomial(y, yd);

    // compute the integral gcd with cofactors:
    // vi = gcd(xi, yi);  vfi*vi == xfi*xi + yfi*yi
    INTPOLY xfi, yfi; INTNT vfi;
    INTPOLY vi = pseudo_gcdex(xi, yi, xfi, yfi, vfi);

    // proceed to vfi*v == xfi*x + yfi*y  with v = gcd(x,y) (unit-normal)
    POLY v = fractionalize_polynomial<POLY>(vi, vi.lcoeff());
    v.simplify_coefficients();
    CGAL_assertion(v.unit_part() == NT(1));
    vfi *= vi.lcoeff(); xfi *= xd; yfi *= yd;

    // compute xf, yf such that gcd(x,y) == v == xf*x + yf*y
    xf = fractionalize_polynomial<POLY>(xfi, vfi);
    yf = fractionalize_polynomial<POLY>(yfi, vfi);
    xf.simplify_coefficients();
    yf.simplify_coefficients();
    return v;
};

} // namespace POLYNOMIAL

/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
 *  \brief compute gcd with cofactors
 *
 *  This function computes the gcd of polynomials \c p1 and \c p2
 *  along with two other polynomials \c f1 and \c f2 such that
 *  gcd(\e p1, \e p2) = <I>f1*p1 + f2*p2</I>. This is called
 *  <I>extended</I> gcd computation, and <I>f1, f2</I> are called
 *  <I>B&eacute;zout factors</I> or <I>cofactors</I>.
 *
 *  Internally, computation is performed ``denominator-free'' if
 *  supported by the coefficient type via \c NiX::Fraction_traits
 *  (using \c pseudo_gcdex() ), otherwise the euclidean remainder
 *  sequence is used.
 *
 *  \pre \c NT must be a \c Field.
 *
 *  The result <I>d</I> is unit-normal,
 *  i.e. <I>d</I><TT>.lcoeff() == NT(1)</TT>.
 *
 */
template <class NT> inline
Polynomial<NT> gcdex(
    Polynomial<NT> p1, Polynomial<NT> p2,
    Polynomial<NT>& f1, Polynomial<NT>& f2
) {
    typedef typename Fraction_traits< Polynomial<NT> >
        ::Is_fraction Is_fraction;
    return POLYNOMIAL::gcdex_(p1, p2, f1, f2, Is_fraction());
};


/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
 *  \brief compute gcd with ``almost'' cofactors
 *
 *  This is a variant of \c exgcd() for use over non-field \c NT.
 *  It computes the gcd of polynomials \c p1 and \c p2
 *  along with two other polynomials \c f1 and \c f2 and a scalar \c v
 *  such that \e v * gcd(\e p1, \e p2) = <I>f1*p1 + f2*p2</I>,
 *  using the subresultant remainder sequence. That \c NT is not a field
 *  implies that one cannot achieve \e v = 1 for all inputs.
 *
 *  \pre \c NT must be a \c UFDomain of scalars (not polynomials).
 *
 *  The result is unit-normal.
 *
 */
template <class NT>
Polynomial<NT> pseudo_gcdex(
#ifdef DOXYGEN_RUNNING
    Polynomial<NT> p1, Polynomial<NT> p2,
    Polynomial<NT>& f2, Polynomial<NT>& f2, NT& v
#else
    Polynomial<NT> x, Polynomial<NT> y,
    Polynomial<NT>& xf, Polynomial<NT>& yf, NT& vf
#endif // DOXYGEN_RUNNING
) {
    /* implemented using the extended subresultant algorithm
     * for gcd computation with Bezout factors
     *
     * To understand this, you need to understand the computation of
     * cofactors as in the basic extended Euclidean algorithm (see
     * the code above of gcdex_(..., Field_tag)), and the subresultant
     * gcd algorithm, see gcd_(..., UFDomain_tag).
     *
     * The crucial point of the combination of both is the observation
     * that the subresultant factor (called rho here) divided out of the
     * new remainder in each step can also be divided out of the
     * cofactors.
     */

    typedef Polynomial<NT> POLY;
    typename Algebraic_structure_traits<NT>::Integral_division idiv;
    typename Algebraic_structure_traits<NT>::Gcd          gcd;

    // handle trivial cases
    if (x.is_zero()) {
        if (y.is_zero()) CGAL_error_msg("gcdex(0,0) is undefined");
        xf = NT(0); yf = NT(1); vf = y.unit_part();
        return y / vf;
    }
    if (y.is_zero()) {
        xf = NT(1); yf = NT(0); vf = x.unit_part();
        return x / vf;
    }
    bool swapped = x.degree() < y.degree();
    if (swapped) { POLY t = x; x = y; y = t; }

    // compute gcd of content
    NT xcont = x.content(); NT ycont = y.content();
    NT gcdcont = gcd(xcont, ycont);

    // compute gcd of primitive parts
    POLY xprim = x / xcont; POLY yprim = y / ycont;
    POLY u = xprim, v = yprim, q, r;
    POLY m11(1), m21(0), m21old;
    NT g(1), h(1), d, rho;
    for (;;) {
        int delta = u.degree() - v.degree();
        POLY::pseudo_division(u, v, q, r, d);
        CGAL_assertion(d == POLYNOMIAL::ipower(v.lcoeff(), delta+1));
        if (r.is_zero()) break;
        rho = g * POLYNOMIAL::ipower(h, delta);
        u = v; v = r / rho;
        m21old = m21; m21 = (d*m11 - q*m21) / rho; m11 = m21old;
        /* The transition from (u, v) to (v, r/rho) corresponds
         * to multiplication with the matrix
         *   __1__ (0 rho)
         *    rho  (d  -q)
         * The comments and correctness arguments from
         * gcdex(..., Field_tag) apply analogously.
         */
        g = u.lcoeff();
        POLYNOMIAL::hgdelta_update(h, g, delta);
        if (r.degree() == 0) break;
    }

    // obtain v == m21*xprim + m22*yprim
    // the correct m21 was already computed above
    POLY m22;
    POLY::euclidean_division(v - m21*xprim, yprim, m22, r);
    CGAL_assertion(r.is_zero());

    // now obtain gcd(x,y) == gcdcont * v/v.content() == (m21*x + m22*y)/denom
    NT vcont = v.content(), vup = v.unit_part();
    v /= vup * vcont; v *= gcdcont;
    m21 *= ycont; m22 *= xcont;
    vf = idiv(xcont, gcdcont) * ycont * (vup * vcont);
    CGAL_assertion(vf * v == m21*x + m22*y);

    // return results
    if (swapped) {
        xf = m22; yf = m21;
    } else {
        xf = m21; yf = m22;
    }
    return v;
}


// 5) square-free factorization


// ### filtered versions #### 

/*! \brief As NiX::square_free_factorization, but filtered by  
 *  NiX::may_have_multiple_root 
 *  
 *  Use this function if the polynomial might be square free. 
 */  
template <class NT, class OutputIterator1, class OutputIterator2> 
inline
int filtered_square_free_factorization(
                                       Polynomial<NT> p,
                                       OutputIterator1 factors,
                                       OutputIterator2 multiplicities)
{
  if(CGAL::POLYNOMIAL::may_have_multiple_factor(p)){
      return POLYNOMIAL::square_free_factorization(p, factors, multiplicities);
  }else{
#if NiX_POLYNOMIAL_GCD_AVOID_CANONICALIZE
      *factors++ = p;
#else
      *factors++ = canonicalize_polynomial(p);
#endif
      *multiplicities++ = 1;
      return 1;
  }
}

/*! \brief As NiX::square_free_factorization_utcf, but filtered by  
 *  NiX::may_have_multiple_root 
 *  
 *  Use this function if the polynomial might be square free. 
 */  
template <class NT,  class OutputIterator1, class OutputIterator2> 
inline
int filtered_square_free_factorization_utcf( const Polynomial<NT>& p, 
                                         OutputIterator1 factors,
                                         OutputIterator2 multiplicities)
{
    if(CGAL::POLYNOMIAL::may_have_multiple_factor(p)){
        return POLYNOMIAL::square_free_factorization_utcf(p,factors,multiplicities);
        
    }else{
#if NiX_POLYNOMIAL_GCD_AVOID_CANONICALIZE
        *factors++ = p;
#else
        *factors++ = canonicalize_polynomial(p);
#endif   
        *multiplicities++ = 1;
        return 1;

    }
}

CGAL_END_NAMESPACE

#endif // CGAL_POLYNOMIAL_POLYNOMIAL_GCD_H

// EOF
