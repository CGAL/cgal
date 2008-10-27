
#include <CGAL/basic.h>


#ifndef CGAL_POLYNOMIAL_CGALi_UTILS_H
#define CGAL_POLYNOMIAL_CGALi_UTILS_H

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

namespace CGALi {

    template <class NT>
    Polynomial<NT> canonicalize_polynomial_(Polynomial<NT> p, CGAL::Tag_true)
    {
        typedef Polynomial<NT> POLY;
        typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient_type IC;
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


    /*! \ingroup CGAL_Polynomial
     *  \relates CGAL::Polynomial
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
     *        \see CGAL/Sqrt_extension.h
     *  The result is uniquely determined by setting the leading coefficient
     *  to the minimal integral rational.
     */
    template <class NT> inline
      Polynomial<NT> canonicalize_polynomial(const Polynomial<NT>& p)
    {
        typedef Polynomial<NT> POLY;
        typedef typename Polynomial_traits_d<POLY>::Innermost_coefficient_type IC;
        typedef typename Algebraic_extension_traits<IC>::Is_extended Is_extended;
    
        if (p.is_zero()) return p;
        return canonicalize_polynomial_(p, Is_extended());
    };

} // namespace CGALi

// div_utfc functions //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


//! divide \c f by \c g with respect to constant factors
/*! This function provides a division of two polynomials, which takes
 *  no care of constant factors of the innermost scalar type.
 *  The result is made unique by canonicalizing it.
 */
    
template <class NT> inline
Polynomial<NT> integral_division_up_to_constant_factor(
    const Polynomial<NT>& f, const Polynomial<NT>& g)
{
  typedef Polynomial<NT> POLY;
  typedef Polynomial_traits_d<POLY> PT;
  typename PT::Integral_division_up_to_constant_factor idiv_utcf; 
  return idiv_utcf(f, g);
}


CGAL_END_NAMESPACE
#endif // CGAL_POLYNOMIAL_CGALi_UTILS_H
