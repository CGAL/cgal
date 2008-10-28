
#include <CGAL/basic.h>


#ifndef CGAL_POLYNOMIAL_CGALi_UTILS_H
#define CGAL_POLYNOMIAL_CGALi_UTILS_H

CGAL_BEGIN_NAMESPACE

//! returns the total degree of the polynomial
template <class Polynomial_d> inline
typename Polynomial_traits_d<Polynomial_d>::Total_degree::result_type
total_degree(const Polynomial_d& polynomial){
    typedef Polynomial_traits_d<Polynomial_d> PT;
    typename PT::Total_degree total_degree;
    return total_degree(polynomial);
}

//! returns the total degree of the polynomial
template <class Polynomial_d> inline
typename Polynomial_traits_d<Polynomial_d>::Canonicalize::result_type
canonicalize(const Polynomial_d& polynomial){
    typedef Polynomial_traits_d<Polynomial_d> PT;
    typename PT::Canonicalize canonicalize;
    return canonicalize(polynomial);
}

template <class Polynomial_d> inline
typename Polynomial_traits_d<Polynomial_d>
::Integral_division_up_to_constant_factor::result_type 
integral_division_up_to_constant_factor(
    const Polynomial_d& f, const Polynomial_d& g)
{
  typedef Polynomial_traits_d<Polynomial_d> PT;
  typename PT::Integral_division_up_to_constant_factor idiv_utcf; 
  return idiv_utcf(f, g);
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


CGAL_END_NAMESPACE
#endif // CGAL_POLYNOMIAL_CGALi_UTILS_H
