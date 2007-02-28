
#include <CGAL/basic.h>


#ifndef CGAL_POLYNOMIAL_UTILS_H
#define CGAL_POLYNOMIAL_UTILS_H

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

CGAL_END_NAMESPACE
#endif // CGAL_POLYNOMIAL_UTILS_H
