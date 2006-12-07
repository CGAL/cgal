// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Arno Eigenwillig <arno@mpi-inf.mpg.de>
//                 Michael Seel <seel@mpi-inf.mpg.de>
//                 Michael Hemmer <hemmer@informatik.uni-mainz.de> 
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file NiX/Polynomial.h
 *  \brief Defines class NiX::Polynomial.
 *  
 *  Polynomials in one variable (or more, by recursion)
 */

#ifndef CGAL_POLYNOMIAL_H
#define CGAL_POLYNOMIAL_H

#include <cstdarg>
#include <cctype>
#include <vector>
#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/mpl/if.hpp>
#include <CGAL/Flattening_iterator.h>
//#include <NiX/Modular.h>

#include <CGAL/Exponent_vector.h>

#include <boost/static_assert.hpp>

#ifdef CGAL_USE_LEDA
#include <LEDA/array.h>
#endif // CGAL_USE_LEDA

#include <CGAL/Polynomial/Polynomial_type.h>
#include <CGAL/Polynomial/Algebraic_structure_traits.h>
#include <CGAL/Polynomial/Real_embeddable_traits.h>
#include <CGAL/Polynomial/Fraction_traits.h>
#include <CGAL/Polynomial/Scalar_factor_traits.h>
// #include <CGAL/Polynomial/Modular_traits.h>

CGAL_BEGIN_NAMESPACE

// TODO: copied from EXACUS/NumeriX/include/NiX/number_type_utils.h
template <typename NT>
inline
NT ipower(const NT& base, int expn) {
    // compute base^expn using square-and-multiply
    CGAL_precondition(expn >= 0);

    // handle trivial cases efficiently
    if (expn == 0) return NT(1);
    if (expn == 1) return base;

    // find the most significant non-zero bit of expn
    int e = expn, msb = 0;
    while (e >>= 1) msb++;

    // computing base^expn by square-and-multiply
    NT res = base;
    int b = 1<<msb;
    while (b >>= 1) { // is there another bit right of what we saw so far?
        res *= res;
        if (expn & b) res *= base;
    }
    return res;
}

// TODO: END included from number_type_utils.h








// Internally, Polynomials also need this:

//! used internally for data exchanged between nesting levels of polynomials
// this traits-class provides
//   a) poly_nesting_depth: a counter for identifying a polynomials nesting level
//   b) Innermost_coefficient: a type which defines the numbertype of the
//        innermost nesting level, which is not a polynomial itself
//   c) Innermost_lcoeff: returns the leading coefficient of the polynomial in its
//       common sense
//   d) Innermost_coefficient_to_polynomial: transforms an innermost
//       numbertype into a polynomial like a typecast

 
// fwd Polynomial_traits_d
template <typename Polynomial_d> class Polynomial_traits_d;

namespace INTERN_POLYNOMIAL {

template <class NT>
Polynomial<NT> canonicalize_polynomial_(Polynomial<NT> p, CGAL::Tag_true)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits_d<POLY>::Innermost_leading_coefficient ilcoeff;
    typename Algebraic_extension_traits<IC>::Normalization_factor nfac;
  
    IC tmp = nfac(ilcoeff(p));
    if(tmp != IC(1)){
        p *= POLY(tmp);
    }
    remove_scalar_factor(p);
    p /= p.unit_part();
    p.simplify_coefficients();

    CGAL_postcondition(nfac(ilcoeff(p)) == IC(1));
    return p;
};

template <class NT>
Polynomial<NT> canonicalize_polynomial_(Polynomial<NT> p, CGAL::Tag_false)
{  
    remove_scalar_factor(p);
    p /= p.unit_part();
    p.simplify_coefficients();
    return p;
};

} // namespace INTERN_POLYNOMIAL

/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
 *  
 *  \brief divide a polynomial \c p by its Scalar_factor and Unit_part
 *  
 *  ...making it a canonical representative of all its constant multiples.
 *  Depending on the number type of the innermost coefficient, this
 *  function does
 *    a) dividing \c p by the leading coefficient in fields
 *    b) dividing \c p by the gcd of all coefficients in UFDomains
 *    c) extending the leading coefficient in Sqrt_extensions, so it
 *        becomes integral, and dividing \c p by the gcd of all scalars
 *        \see NiX/Sqrt_extension.h
 *  The result is uniquely determined by setting the leading coefficient
 *  to the minimal integral rational.
 */
template <class NT> inline
  Polynomial<NT> canonicalize_polynomial(const Polynomial<NT>& p)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient IC;
    typedef typename Algebraic_extension_traits<IC>::Is_extended Is_extended;

    if (p.is_zero()) return p;
    return canonicalize_polynomial_(p, Is_extended());
};

namespace INTERN_POLYNOMIAL {

// Polynomial<NT> / Polynomial<NT>  -  coefficient type is extended
template <class NT>
Polynomial<NT> div_utcf_(
    Polynomial<NT> f, const Polynomial<NT>& g, bool, CGAL::Tag_true)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits_d<POLY>::Innermost_leading_coefficient ilcoeff;
    typename Polynomial_traits_d<POLY>::Innermost_coefficient_to_polynomial ictp;
    typename Polynomial_traits_d<POLY>::Innermost_coefficient_begin begin;
    typename Polynomial_traits_d<POLY>::Innermost_coefficient_end end;
    typename Algebraic_extension_traits<IC>::Denominator_for_algebraic_integers dfai;

    IC tmp = ilcoeff(g);
    tmp *= dfai(begin(g), end(g));
    f *= ictp(tmp);
    return canonicalize_polynomial(f / g);
}

// Polynomial<NT> / Polynomial<NT>  -  coefficient type is NOT extended
template <class NT>
Polynomial<NT> div_utcf_(
    Polynomial<NT> f, const Polynomial<NT>& g, bool is_canonicalized, CGAL::Tag_false)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits_d<POLY>::Innermost_leading_coefficient ilcoeff;
    typename Polynomial_traits_d<POLY>::Innermost_coefficient_to_polynomial ictp;

    if (!is_canonicalized) {
        IC lcoeff = ilcoeff(g);
        f *= ictp(lcoeff);
    }
    return canonicalize_polynomial(f / g);
}

// Polynomial<NT> / NT  -  NT is already the coefficient type and is extended
template <class NT>
Polynomial<NT> div_utcf_NT_is_IC(
    Polynomial<NT> f, const NT& g, CGAL::Tag_false)
{
    return canonicalize_polynomial(f);
}

// Polynomial<NT> / NT  -  NT is again a polynomial  -   coefficient type is extended
template <class NT, class Is_nested>
Polynomial<NT> div_utcf_NT_is_IC(
    Polynomial<NT> f, const NT& g, Is_nested)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits_d<POLY>::Innermost_leading_coefficient ilcoeff;
    typename Polynomial_traits_d<POLY>::Innermost_coefficient_to_polynomial ictp;
    typename Polynomial_traits_d<NT>::Innermost_coefficient_begin begin;
    typename Polynomial_traits_d<NT>::Innermost_coefficient_end end;
    typename Algebraic_extension_traits<IC>::Denominator_for_algebraic_integers dfai;

    IC tmp = ilcoeff(g);
    tmp *= dfai(begin(g), end(g));
    f *= ictp(tmp);
    return canonicalize_polynomial(f / g);
}

// Polynomial<NT> / NT  -  coefficient type is extended
template <class NT> inline
Polynomial<NT> div_utcf_(
    const Polynomial<NT>& f, const NT& g, bool, CGAL::Tag_true)
{
    typedef CGAL::Boolean_tag< (Polynomial_traits_d<NT>::d >= 2) > Is_nested;
    return div_utcf_NT_is_IC(f, g, Is_nested() );
}

// Polynomial<NT> / NT  -  coefficient type is NOT extended
template <class NT>
Polynomial<NT> div_utcf_(
    Polynomial<NT> f, const NT& g, bool is_canonicalized, CGAL::Tag_false)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient IC;
    typename Polynomial_traits_d<POLY>::Innermost_leading_coefficient ilcoeff;
   
    if (!is_canonicalized) {
        IC lcoeff = ilcoeff(g);
        f *= POLY(lcoeff);
    }
    return canonicalize_polynomial(f / g);
}

} // namespace INTERN_POLYNOMIAL

//! divide \c f by \c g with respect to constant factors
/*! This function provides a division of two polynomials, which takes
 *  no care of constant factors of the innermost scalar type.
 *  The boolean parameter decides whether the divisor has already been
 *  canonicalized due to running time optimisation.
 *  The result is made unique by canonicalizing it.
 */

template <class NT> inline
Polynomial<NT> div_utcf(
    const Polynomial<NT>& f, const Polynomial<NT>& g, bool is_canonicalized = false)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient IC;
    typedef typename Algebraic_extension_traits<IC>::Is_extended Is_extended;

    return div_utcf_(f, g, is_canonicalized, Is_extended());
}

//! overloaded version for divisors with a by one lower nesting level
template <class NT> inline
Polynomial<NT> div_utcf(
    const Polynomial<NT>& f, const NT& g, bool is_canonicalized = false)
{
    typedef Polynomial<NT> POLY;
    typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient IC;
    typedef typename Algebraic_extension_traits<IC>::Is_extended Is_extended;

    return div_utcf_(f, g, is_canonicalized, Is_extended());
}

//
// Arithmetic Operators, Part III:
// implementation of unary operators and three-address arithmetic
// by friend functions
//

template <class NT> inline
Polynomial<NT> operator + (const Polynomial<NT>& p) {
    CGAL_precondition(p.degree() >= 0);
    return p;
}

template <class NT> inline
Polynomial<NT> operator - (const Polynomial<NT>& p) {
    CGAL_precondition(p.degree()>=0);
    Polynomial<NT> res(p.coeffs().begin(),p.coeffs().end());
    typename Polynomial<NT>::iterator it, ite=res.coeffs().end();
    for(it=res.coeffs().begin(); it!=ite; ++it) *it = -*it;
    return res;
}

template <class NT> inline
Polynomial<NT> operator + (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2)
{ 
    typedef typename Polynomial<NT>::size_type size_type;
    CGAL_precondition(p1.degree()>=0 && p2.degree()>=0);
    bool p1d_smaller_p2d = p1.degree() < p2.degree();
    int min,max,i;
    if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
    else                 { max = p1.degree(); min = p2.degree(); }
    INTERN_POLYNOMIAL::Creation_tag TAG;
    Polynomial<NT> p(TAG, size_type(max + 1));
    for (i = 0; i <= min; ++i ) p.coeff(i) = p1[i]+p2[i];
    if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.coeff(i)=p2[i];
    else /* p1d >= p2d */ for (; i <= max; ++i ) p.coeff(i)=p1[i];
    p.reduce();
    return p;
}

template <class NT> inline
Polynomial<NT> operator - (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2)
{ 
    typedef typename Polynomial<NT>::size_type size_type;
    CGAL_precondition(p1.degree()>=0 && p2.degree()>=0);
    bool p1d_smaller_p2d = p1.degree() < p2.degree();
    int min,max,i;
    if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
    else                 { max = p1.degree(); min = p2.degree(); }
    INTERN_POLYNOMIAL::Creation_tag TAG;
    Polynomial<NT> p(TAG, size_type(max + 1));
    for (i = 0; i <= min; ++i ) p.coeff(i)=p1[i]-p2[i];
    if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.coeff(i)= -p2[i];
    else /* p1d >= p2d */ for (; i <= max; ++i ) p.coeff(i)=  p1[i];
    p.reduce();
    return p;
}

template <class NT> inline
Polynomial<NT> operator * (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2)
{
    typedef typename Polynomial<NT>::size_type size_type;
    CGAL_precondition(p1.degree()>=0 && p2.degree()>=0);
    INTERN_POLYNOMIAL::Creation_tag TAG;
    Polynomial<NT>  p(TAG, size_type(p1.degree()+p2.degree()+1) ); 
    // initialized with zeros
    for (int i=0; i <= p1.degree(); ++i)
        for (int j=0; j <= p2.degree(); ++j)
            p.coeff(i+j) += (p1[i]*p2[j]); 
    p.reduce();
    return p;
}

#ifndef NiX_POLY_USE_SLOW_DIVISION

template <class NT> inline
Polynomial<NT> operator / (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2)
{
    typedef Algebraic_structure_traits< Polynomial<NT> > AST; 
    // Precondition: q with p1 == p2 * q must exist within NT[x].
    // If this holds, we can perform Euclidean division even over a ring NT
    // Proof: The quotients of each division that occurs are precisely
    //   the terms of q and hence in NT.
    CGAL_precondition(!p2.is_zero());
    if (p1.is_zero()) return p1;
    Polynomial<NT> q, r;
    Polynomial<NT>::euclidean_division(p1, p2, q, r);
    // TODO: Replace by correct Makro
    CGAL_postcondition( !AST::Is_exact::value || p2 * q == p1);
    
    return q;
}

#else

template <class NT> inline
Polynomial<NT> operator / (const Polynomial<NT>& p1, 
                           const Polynomial<NT>& p2)
{
    typedef typename Algebraic_structure_traits<NT>::Algebric_structure_tag Algebra_type;
    return division(p1, p2, Algebra_type());
}

template <class NT> inline
Polynomial<NT> division(const Polynomial<NT>& p1, 
                        const Polynomial<NT>& p2,
                        Field_tag)
{
    typedef Algebraic_structure_traits<NT> AST;
    CGAL_precondition(!p2.is_zero());
    if (p1.is_zero()) return p1;
    Polynomial<NT> q,r;
    Polynomial<NT>::euclidean_division(p1,p2,q,r);

    CGAL_postcondition( !AST::Is_exact::value  || p2 * q == p1);
    return q;
}

template <class NT> inline
Polynomial<NT> division(const Polynomial<NT>& p1, 
                        const Polynomial<NT>& p2,
                        Integral_domain_tag)
{
    typedef Algebraic_structure_traits<NT> AST;
    CGAL_precondition(!p2.is_zero());
    if ( p1.is_zero() ) return p1;
    Polynomial<NT> q,r; NT D;
    Polynomial<NT>::pseudo_division(p1,p2,q,r,D);
    q/=D;

    CGAL_postcondition( !AST::Is_exact::value || p2 * q == p1);
    return q;
}

#endif // NiX_POLY_USE_SLOW_DIVISION

//
// Arithmetic Operators, Part IV:
// Mixed-mode three-address arithmetic on top of two-address operators
//

// lefthand side
template <class NT> inline
Polynomial<NT> operator + (const NT& num, Polynomial<NT> p2)
{ p2 += num; return p2; }
template <class NT> inline
Polynomial<NT> operator - (const NT& num, Polynomial<NT> p2)
{ p2 -= num; return -p2; }
template <class NT> inline
Polynomial<NT> operator * (const NT& num, Polynomial<NT> p2)
{ p2 *= num; return p2; }
template <class NT> inline
Polynomial<NT> operator / (const NT& num, const Polynomial<NT>& p2) {
    CGAL_precondition(p2.degree() == 0);
    CGAL_precondition(p2[0] != NT(0));
    typename Algebraic_structure_traits<NT>::Integral_division idiv;
    return Polynomial<NT>(idiv(num, p2[0]));
}

// righthand side
template <class NT> inline
Polynomial<NT> operator + (Polynomial<NT> p1, const NT& num)
{ p1 += num; return p1; }
template <class NT> inline
Polynomial<NT> operator - (Polynomial<NT> p1, const NT& num)
{ p1 -= num; return p1; }
template <class NT> inline
Polynomial<NT> operator * (Polynomial<NT> p1, const NT& num)
{ p1 *= num; return p1; }
template <class NT> inline
Polynomial<NT> operator / (Polynomial<NT> p1, const NT& num)
{ p1 /= num; return p1; }


//
// Comparison Operators
//

// polynomials only
template <class NT> inline
bool operator == (const Polynomial<NT>& p1, const Polynomial<NT>& p2) {
    CGAL_precondition(p1.degree() >= 0);
    CGAL_precondition(p2.degree() >= 0);
    if (p1.is_identical(p2)) return true;
    if (p1.degree() != p2.degree()) return false;
    for (int i = p1.degree(); i >= 0; i--) if (p1[i] != p2[i]) return false;
    return true;
}
template <class NT> inline
bool operator != (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return !(p1 == p2); }    
template <class NT> inline
bool operator < (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return ( p1.compare(p2) < 0 ); }    
template <class NT> inline
bool operator <= (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return ( p1.compare(p2) <= 0 ); }    
template <class NT> inline
bool operator > (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return ( p1.compare(p2) > 0 ); }    
template <class NT> inline
bool operator >= (const Polynomial<NT>& p1, const Polynomial<NT>& p2)
{ return ( p1.compare(p2) >= 0 ); }    

// mixed-mode, lefthand side
template <class NT> inline
bool operator == (const NT& num, const Polynomial<NT>& p)  {
    CGAL_precondition(p.degree() >= 0);
    return p.degree() == 0 && p[0] == num;
}
template <class NT> inline
bool operator != (const NT& num, const Polynomial<NT>& p) 
{ return !(num == p);}
template <class NT> inline
bool operator < (const NT& num, const Polynomial<NT>& p) 
{ return ( p.compare(num) > 0 );}
template <class NT> inline
bool operator <=  (const NT& num, const Polynomial<NT>& p) 
{ return ( p.compare(num) >= 0 );}
template <class NT> inline
bool operator > (const NT& num, const Polynomial<NT>& p) 
{ return ( p.compare(num) < 0 );}
template <class NT> inline
bool operator >= (const NT& num, const Polynomial<NT>& p) 
{ return ( p.compare(num) <= 0 );}

// mixed-mode, righthand side
template <class NT> inline
bool operator == (const Polynomial<NT>& p, const NT& num) 
{ return num == p; }
template <class NT> inline
bool operator != (const Polynomial<NT>& p, const NT& num) 
{ return !(num == p); }
template <class NT> inline
bool operator < (const Polynomial<NT>& p, const NT& num) 
{ return ( p.compare(num) < 0 );}
template <class NT> inline
bool operator <= (const Polynomial<NT>& p, const NT& num) 
{ return ( p.compare(num) <= 0 );}
template <class NT> inline
bool operator > (const Polynomial<NT>& p, const NT& num) 
{ return ( p.compare(num) > 0 );}
template <class NT> inline
bool operator >= (const Polynomial<NT>& p, const NT& num) 
{ return ( p.compare(num) >= 0 );}


//
// I/O Operations
//

/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
 *  \brief output \c p to \c os
 *
 *  Output \c p in a format as specified by
 *  \c LiS::get_mode(os), see \link LiS_io LiS I/O Support \endlink.
 *  Currently, the output for \c LiS::IO::BINARY happens to be
 *  identical to \c LiS::IO::ASCII.
 */
template <class NT>
std::ostream& operator << (std::ostream& os, const Polynomial<NT>& p) {
    switch(CGAL::get_mode(os)) {
    case CGAL::IO::PRETTY:
        p.output_maple(os); break;
    default:
        p.output_ascii(os); break;
    }
    return os;
}

/*! \ingroup NiX_Polynomial
 *  \relates NiX::Polynomial
 *  \brief try to read a polynomial from \c is into \c p
 *
 *  \c is must be in a mode that supports input of polynomials
 *  (\c LiS::IO::ASCII or \c LiS::IO::BINARY) and the input from
 *  \c is must have the format of output to a stream of the same mode.
 */
template <class NT>
std::istream& operator >> (std::istream& is, Polynomial<NT>& p) {
    CGAL_precondition(!CGAL::is_pretty(is));
    p = Polynomial<NT>::input_ascii(is);
    return is;
}


template <class NT> inline
void print_maple_monomial(std::ostream& os, const NT& coeff,
                          const char *var, int expn)
{
    if (expn == 0 || coeff != NT(1)) {
        os << CGAL::oformat(coeff, Parens_as_product_tag());
        if (expn >= 1) os << "*";
    }
    if (expn >= 1) {
        os << var;
        if (expn > 1) os << "^" << CGAL::oformat(expn);
    }
}

template <class NT>
void Polynomial<NT>::output_maple(std::ostream& os) const {
    const Polynomial<NT>& p = *this;
    const char *varname;
    char vnbuf[42];
    
    // use variable names x, y, z, w1, w2, w3, ...
    if (Polynomial_traits_d<NT>::d < 3) {
        static const char *varnames[] = { "x", "y", "z" };
        varname = varnames[Polynomial_traits_d<NT>::d];
    } else {
        sprintf(vnbuf, "w%d", Polynomial_traits_d<NT>::d - 2);
        varname = vnbuf;
    }
    
    int i = p.degree();
    print_maple_monomial(os, p[i], varname, i);
    while (--i >= 0) {
        if (p[i] != NT(0)) {
            os << " + ";
            print_maple_monomial(os, p[i], varname, i);
        }
    }
}

template <class NT>
void Polynomial<NT>::output_ascii(std::ostream &os) const {
    const Polynomial<NT> &p = *this;
    if (p.is_zero()) { os << "P[0 (0," << oformat(NT(0)) << ")]"; return; }

    os << "P[" << oformat(p.degree());
    for (int i = 0; i <= p.degree(); i++) {
        if (p[i] != NT(0)) os << "(" << CGAL::oformat(i) << ","
                              << CGAL::oformat(p[i]) << ")";
    }
    os << "]";
}

template <class NT>
void Polynomial<NT>::output_benchmark(std::ostream &os) const {
    const Polynomial<NT> &p = *this;
    if (p.is_zero()) { 
        os << "Polynomial_1(0)"; 
        return; 
    }
    os << "Polynomial_1(";
    for (int i = 0; i <= p.degree(); i++) {
        os << CGAL::oformat(p[i]);
        if (i != p.degree()) {
            os << ",";
        }
    }
    os << ")";
}

// Moved to internal namespace because of name clashes
// TODO: Is this OK?
namespace INTERN_POLYNOMIAL {

  inline static void swallow(std::istream &is, char d) {
      char c;
      do c = is.get(); while (isspace(c));
      if (c != d) CGAL_assertion_msg( false, "input error: unexpected character in polynomial");
  }
} // namespace INTERN_POLYNOMIAL

template <class NT>
Polynomial<NT> Polynomial<NT>::input_ascii(std::istream &is) {
    char c;
    int degr = -1, i;

    INTERN_POLYNOMIAL::swallow(is, 'P');
    INTERN_POLYNOMIAL::swallow(is, '[');
    is >> CGAL::iformat(degr);
    if (degr < 0) {
        CGAL_assertion_msg( false, "input error: negative degree of polynomial specified");
    }
    INTERN_POLYNOMIAL::Creation_tag TAG;
    Polynomial<NT> p(TAG, degr+1);

    do c = is.get(); while (isspace(c));
    do {
        if (c != '(') CGAL_assertion_msg( false, "input error: ( expected");
        is >> CGAL::iformat(i);
        if (!(i >= 0 && i <= degr && p[i] == NT(0))) {
            CGAL_assertion_msg( false, "input error: invalid exponent in polynomial");
        };
        INTERN_POLYNOMIAL::swallow(is, ',');
        is >> CGAL::iformat(p.coeff(i));
        INTERN_POLYNOMIAL::swallow(is, ')');
        do c = is.get(); while (isspace(c));
    } while (c != ']');

    p.reduce();
    p.simplify_coefficients();
    return p;
}



//
// Non-Member Functions
//



//! return an upper bound on the absolute value of all real roots of \c P.
/*! The upper bound is a power of two. Only works for univariate polynomials.
 *  \pre \c NT must be \c RealComparable.
 *  \relates NiX::Polynomial
 */
template <class NT>
NT weak_upper_root_bound(const Polynomial<NT>& P) { 
    // code comes from Kurt Mehlhorn
    // see [Mignotte, 1992], p.144 for a proof
    CGAL_precondition(Polynomial_traits_d<NT>::d == 0);
    typename Real_embeddable_traits<NT>::Abs abs;
    const int n = P.degree();
    NT x(1); 
    NT val;
    for (;;) {
        val = -abs(P[n]);
        for (int i = n-1; i >= 0; i--) {
            val = val*x + abs(P[i]);
        }
        if (val < NT(0)) return x;
        x *= NT(2);
    }
}

//! return the number of sign variations in the coefficient sequence of \c P.
/*! This is the number of sign changes (+ to - or - to +) in the
 *  coefficient sequence of the polynomial, ignoring zeroes.
 *  Only meaningful for univariate polynomials.
 *  \pre \c NT must be \c RealComparable.
 *  \relates NiX::Polynomial
 */
template <class NT>
int sign_variations(const Polynomial<NT>& P) { 
    typename Real_embeddable_traits<NT>::Sign sign;
    const int n = P.degree();
    int variations = 0;
    int old_sign = sign(P[n]); // never zero unless P is zero
    for (int i = n-1; i >= 0; i--) {
        int s = sign(P[i]);
        if (s == 0) continue;
        if (old_sign != s) {
            old_sign = s;
            variations++;
        }
    }
    return variations;
}


//
// Algebraically non-trivial operations
//

// 1) Euclidean and pseudo-division of polynomials
// (implementation of static member functions)

template <class NT>
void Polynomial<NT>::euclidean_division(
    const Polynomial<NT>& f, const Polynomial<NT>& g,
    Polynomial<NT>& q, Polynomial<NT>& r)
{
    typedef Algebraic_structure_traits<NT> AST;
    typename AST::Integral_division idiv;
    int fd = f.degree(), gd = g.degree();
    if ( fd < gd ) {
        q = Polynomial<NT>(NT(0)); r = f;

        CGAL_postcondition( !AST::Is_exact::value || f == q*g + r); 
        return;
    }
    // now we know fd >= gd 
    int qd = fd-gd, delta = qd+1, rd = fd;

    INTERN_POLYNOMIAL::Creation_tag TAG;    
    q = Polynomial<NT>(TAG, delta ); 
    r = f; r.copy_on_write();
    while ( qd >= 0 ) {
        NT Q = idiv(r[rd], g[gd]);
        q.coeff(qd) += Q;
        r.minus_offsetmult(g,Q,qd);
        r.simplify_coefficients();
        if (r.is_zero()) break;
        rd = r.degree();
        qd = rd - gd;
    }
    q.simplify_coefficients();

    CGAL_postcondition( !AST::Is_exact::value || f == q*g + r);
}

#ifndef NiX_POLY_USE_OLD_PSEUDODIV

template <class NT>
void Polynomial<NT>::pseudo_division(
    const Polynomial<NT>& A, const Polynomial<NT>& B,
    Polynomial<NT>& Q, Polynomial<NT>& R, NT& D)
{
    typedef Algebraic_structure_traits<NT> AST;
    // pseudo-division with incremental multiplication by lcoeff(B)
    // see [Cohen, 1993], algorithm 3.1.2

    CGAL_precondition(!B.is_zero());
    int delta = A.degree() - B.degree();

    if (delta < 0 || A.is_zero()) {
        Q = Polynomial<NT>(NT(0)); R = A; D = NT(1);
       
        CGAL_postcondition( !AST::Is_exact::value || Polynomial<NT>(D)*A == Q*B + R);
        return;
    }
    const NT d = B.lcoeff();
    int e = delta + 1;
    D = ipower(d, e);
    INTERN_POLYNOMIAL::Creation_tag TAG;
    Q = Polynomial<NT>(TAG, e);
    R = A; R.copy_on_write(); R.simplify_coefficients();

    // invariant: d^(deg(A)-deg(B)+1 - e) * A == Q*B + R
    do { // here we have delta == R.degree() - B.degree() >= 0 && R != 0
        NT lR = R.lcoeff();
        for (int i = delta+1; i <= Q.degree(); i++) Q.coeff(i) *= d;
        Q.coeff(delta) = lR;              // Q = d*Q + lR * X^delta
        for (int i = 0; i <= R.degree(); i++) R.coeff(i) *= d;
        R.minus_offsetmult(B, lR, delta); // R = d*R - lR * X^delta * B
        R.simplify_coefficients();
        e--;
        delta = R.degree() - B.degree();
    } while (delta > 0 || delta == 0 && !R.is_zero());
    // funny termination condition because deg(0) = 0, not -\infty

    NT q = ipower(d, e);
    Q *= q; Q.simplify_coefficients();
    R *= q; R.simplify_coefficients();

    CGAL_postcondition( !AST::Is_exact::value || Polynomial<NT>(D)*A == Q*B + R);
}

#else

template <class NT>
void Polynomial<NT>::pseudo_division(
    const Polynomial<NT>& f, const Polynomial<NT>& g, 
    Polynomial<NT>& q, Polynomial<NT>& r, NT& D)
{
    typedef Algebraic_structure_traits<NT> AST;
    // pseudo-division with one big multiplication with lcoeff(g)^{...}
    typename Algebraic_structure_traits<NT>::Integral_division idiv;

    int fd=f.degree(), gd=g.degree();
    if ( fd < gd ) {
        q = Polynomial<NT>(NT(0)); r = f; D = NT(1); 

        CGAL_postcondition( !AST::Is_exact::value  || Polynomial<NT>(D)*f==q*g+r);
        return;
    }
    // now we know rd >= gd 
    int qd = fd-gd, delta = qd+1, rd = fd;
    INTERN_POLYNOMIAL::Creation_tag TAG;
    q = Polynomial<NT>(TAG, delta );
    NT G = g[gd]; // highest order coeff of g
    D = ipower(G, delta);
    Polynomial<NT> res = D*f;
    res.simplify_coefficients();
    while ( qd >= 0 ) {
        NT F = res[rd];    // highest order coeff of res
        NT t = idiv(F, G); // sure to be integral by multiplication of D
        q.coeff(qd) = t;   // store q coeff
        res.minus_offsetmult(g,t,qd); 
        res.simplify_coefficients();
        if (res.is_zero()) break;
        rd = res.degree();
        qd = rd - gd;
    }
    r = res; // already simplified
    q.simplify_coefficients();

    CGAL_postcondition( !AST::Is_exact::value  || Polynomial<NT>(D)*f==q*g+r);
}

#endif // NiX_POLY_USE_OLD_PSEUDODIV

// The following former parts of Polynomial.h have been moved to the new files
// <NiX/polynomial_gcd.h> (items 2,3,5) and <NiX/prs_resultant.h> (item 4):
//
// 2) gcd (basic form without cofactors)
// 3) extended gcd computation (with cofactors)
// 4) resultant computation from polynomial remainder sequences (PRS)
// 5) square-free factorization

// This subroutine has been retained here for use in both new files.
namespace INTERN_POLYNOMIAL {
template <class NT> inline
void hgdelta_update(NT& h, const NT& g, int delta) {
    typename Algebraic_structure_traits<NT>::Integral_division idiv;

    // compute h = h^(1-delta) * g^delta
    switch (delta) {
    case 0:
        // h = h;
        break;
    case 1:
        h = g;
        break;
    default:
        h = idiv(ipower(g, delta), ipower(h, delta-1));
        break;
    }
}
} // namespace INTERN_POLYNOMIAL

// } // namespace NiX
CGAL_END_NAMESPACE

#include <CGAL/polynomial_gcd.h> // used above for NT_traits<Poly...>::Gcd
#include <CGAL/prs_resultant.h>  // for compatibility

CGAL_BEGIN_NAMESPACE

// Cofraction_traits added by Michael Hemmer
/*namespace NiX{
namespace Intern{
    template <class NT, class TAG> class Cofraction_traits_base;

    template <class NT_, class TAG>
    class Cofraction_traits_base<NiX::Polynomial<NT_>, TAG > {
        typedef NT_ NT;
    public:
        typedef Polynomial<NT_>                        Numerator_type;
        typedef ::LiS::False_tag                       Is_composable;
        typedef ::LiS::Null_tag                       Denominator_type;
        typedef ::LiS::Null_tag                       Type;
        typedef ::LiS::Null_tag                       Compose;
    };
    
    template <class NT_>
    class Cofraction_traits_base<NiX::Polynomial<NT_>, LiS::True_tag > {
        typedef NT_ NT;
        typedef NiX::Cofraction_traits<NT_> CFT_NT;
    public:
        typedef Polynomial<NT>                             Numerator_type;
        typedef ::LiS::True_tag                            Is_composable;
        typedef typename CFT_NT::Denominator_type          Denominator_type;
        typedef Polynomial<typename CFT_NT::Type> Type;
        
        class Compose {
        public:
            //! first argument type
            typedef Numerator_type   first_argument_type;
            //! second argument type
            typedef Denominator_type second_argument_type;
            //! result type
            typedef Type    result_type;
            //! Compose fraction
            Type operator() (Numerator_type num, 
                                      Denominator_type den 
                                      = Denominator_type(1)){
                Type tmp1; NiX::convert_to(num,tmp1);
                Type tmp2; NiX::convert_to(den,tmp2);
                return tmp1/tmp2;
            }
        };
    };
} //namespace Intern
*/
/*! \ingroup NiX_Cofraction_traits_specs
 *  \brief Specialization of Cofraction_traits for NiX::Polynomial<NT>.
 */
//template<class NT>
/*class Cofraction_traits<Polynomial<NT> > :
    public Intern::Cofraction_traits_base<
      Polynomial<NT>,
      typename Cofraction_traits<NT>::Is_composable>{
    //nothing new
};*/ 

template <class COEFF>
struct Needs_parens_as_product<Polynomial<COEFF> >{
    typedef Polynomial<COEFF> Poly;
    bool operator()(const Poly& x){ return (x.degree() > 0); }
};


// COERCION_TRAITS BEGIN 

//Coercion_traits_polynomial-----------------------------------
// If there is a Polynomial_traits, valid for more than one Polynomial
// class this part should be adapted, using a Polynomial_traits 
// and the nesting_depth 
template <class A,class B>
class Coercion_traits_for_level<Polynomial<A>, Polynomial<B>, CTL_POLYNOMIAL >{
    typedef Coercion_traits<A,B> CT;            
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;
    typedef Polynomial<typename CT::Type> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const Polynomial<A>& poly) const { 
            typename CT::Cast cast; 
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                    ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
        Type operator()(const Polynomial<B>& poly) const {  
            typename CT::Cast cast;  
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                    ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
    }; 
};
        
template <class A,class B>
class Coercion_traits_for_level<Polynomial<A>,B ,CTL_POLYNOMIAL >{
    typedef Coercion_traits<A,B> CT;
public:
    typedef CGAL::Tag_true  Are_explicit_interoperable;
    typedef CGAL::Tag_false Are_implicit_interoperable;

    typedef Polynomial<typename CT::Type> Type;
    struct Cast{                                      
        typedef Type result_type;                               
        Type operator()(const Polynomial<A>& poly) const {
            typename CT::Cast cast;
            return Type(::boost::make_transform_iterator(poly.begin(),cast),
                       ::boost::make_transform_iterator(poly.end()  ,cast));
        } 
        Type operator()(const B& x) const {
            typename CT::Cast cast;
            return Type(cast(x));
        } 
    };                                                        
}; 
template <class A,class B> 
class Coercion_traits_for_level<B,Polynomial<A>,CTL_POLYNOMIAL  >
    :public Coercion_traits_for_level<Polynomial<A>,B,CTL_POLYNOMIAL >
{};

// COERCION_TRAITS END



CGAL_END_NAMESPACE

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/polynomial_utils.h>

//
// trailing documentation
//

// Literature reference
//
// [Akritas, 1989]
//      Alkiviadis G. Akritas
//      Elements of Computer Algebra With Applications
//      Wiley, New York, 1989.
//
// [Cohen, 1993]
//      Cohen, Henri
//      A Course in Computational Algebraic Number Theory
//      Springer GTM 138, 1993
//
// [Cox et al, 1997]
//      David Cox; John Little; Donal O'Shea
//      Ideals, Varieties, and Algorithms
//      2nd ed., Springer UTM, 1997
//
// [Geddes et al, 1992]
//      Geddes, Keith O. and Czapor, Stephen R. and Labahn, George
//      Algorithms for Computer Algebra
//      Kluwer, 1992
//
// [Mignotte, 1992]
//      Mignotte, Maurice
//      Mathematics for Computer Algebra
//      Springer, 1992
//
// [PARI]
//      (a computer algebra system by Henri Cohen and collaborators)
//      http://www.parigp-home.de/

#endif  // NiX_POLYNOMIAL_H

// EOF
