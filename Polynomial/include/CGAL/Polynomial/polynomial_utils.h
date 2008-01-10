
#include <CGAL/basic.h>


#ifndef CGAL_POLYNOMIAL_POLYNOMIAL_UTILS_H
#define CGAL_POLYNOMIAL_POLYNOMIAL_UTILS_H

CGAL_BEGIN_NAMESPACE

//! returns the total degree of the polynomial
template <class Polynomial> 
int total_degree(const Polynomial& polynomial){
    typedef Polynomial_traits_d<Polynomial> PT;
    typename PT::Total_degree total_degree;
    return total_degree(polynomial);
}


// sign() forwarded to the sign() member function
template <class NT> inline 
CGAL::Sign sign(const Polynomial<NT>& p) { return p.sign(); }

// the non-member variants of diff() etc.
template <class NT> inline
Polynomial<NT> diff(const Polynomial<NT>& p)
{ Polynomial<NT> q(p); q.diff(); return q; }

template<class NT> inline
Polynomial<NT> scale_up(const Polynomial<NT>& p, const NT& a)
{ Polynomial<NT> q(p); q.scale_up(a); return q; }

template<class NT> inline
Polynomial<NT> scale_down(const Polynomial<NT>& p, const NT& b)
{ Polynomial<NT> q(p); q.scale_down(b); return q; }

template<class NT> inline
Polynomial<NT> scale(const Polynomial<NT>& p, const NT& a, const NT& b)
{ Polynomial<NT> q(p); q.scale(a, b); return q; }

template<class NT> inline
Polynomial<NT> translate_by_one(const Polynomial<NT>& p)
{ Polynomial<NT> q(p); q.translate_by_one(); return q; }

template<class NT> inline
Polynomial<NT> translate(const Polynomial<NT>& p, const NT& c)
{ Polynomial<NT> q(p); q.translate(c); return q; }

template<class NT> inline
Polynomial<NT> translate(const Polynomial<NT>& p, const NT& a, const NT& b)
{ Polynomial<NT> q(p); q.translate(a, b); return q; }

template<class NT> inline
Polynomial<NT> reversal(const Polynomial<NT>& p)
{ Polynomial<NT> q(p); q.reversal(); return q; }



// Canonicalize_polynomial functions ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


namespace POLYNOMIAL {

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

} // namespace POLYNOMIAL

// div_utfc functions //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace POLYNOMIAL {

    // Polynomial<NT> / Polynomial<NT>  -  coefficient type is extended
    template <class NT>
    Polynomial<NT> div_utcf_(
        Polynomial<NT> f, const Polynomial<NT>& g, bool, CGAL::Tag_true)
    {
        typedef Polynomial<NT> POLY;
        typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient IC;
        typename Polynomial_traits_d<POLY>::Innermost_leading_coefficient ilcoeff;
        typename Polynomial_traits_d<POLY>::Innermost_coefficient_begin begin;
        typename Polynomial_traits_d<POLY>::Innermost_coefficient_end end;
        typename Algebraic_extension_traits<IC>::Denominator_for_algebraic_integers dfai;
    
        IC tmp = ilcoeff(g);
        tmp *= dfai(begin(g), end(g));
        f *= POLY(tmp);
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
    
        if (!is_canonicalized) {
            IC lcoeff = ilcoeff(g);
            f *= POLY(lcoeff);
        }
        return canonicalize_polynomial(f / g);
    }
    
    // Polynomial<NT> / NT  -  NT is already the coefficient type and is extended
    template <class NT>
    Polynomial<NT> div_utcf_NT_is_IC(
        Polynomial<NT> f, const NT&, CGAL::Tag_false)
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
        typename Polynomial_traits_d<NT>::Innermost_coefficient_begin begin;
        typename Polynomial_traits_d<NT>::Innermost_coefficient_end end;
        typename Algebraic_extension_traits<IC>::Denominator_for_algebraic_integers dfai;
    
        IC tmp = ilcoeff(g);
        tmp *= dfai(begin(g), end(g));
        f *= POLY(tmp);
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

} // namespace POLYNOMIAL

CGAL_END_NAMESPACE
#endif // CGAL_POLYNOMIAL_POLYNOMIAL_UTILS_H
