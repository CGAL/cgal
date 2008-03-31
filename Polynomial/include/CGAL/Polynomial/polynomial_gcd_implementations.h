#ifndef CGAL_POLYNOMIAL_GCD_IMPLEMENTATIONS_H
#define CGAL_POLYNOMIAL_GCD_IMPLEMENTATIONS_H


#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Polynomial/polynomial_utils.h>
#include <CGAL/Polynomial/hgdelta_update.h>
#include <CGAL/Polynomial/polynomial_gcd.h>


namespace CGAL {  
namespace CGALi {

template <class NT>
inline
Polynomial<NT> gcd_utcf_UFD(
        Polynomial<NT> p1, Polynomial<NT> p2
) {
    typedef CGAL::Polynomial_traits_d<Polynomial<NT> > PT;
    // implemented using the subresultant algorithm for gcd computation
    // see [Cohen, 1993], algorithm 3.3.1
    // handle trivial cases
    if (p1.is_zero())
        if (p2.is_zero()) return Polynomial<NT>(NT(1));
        else {
            return CGAL::CGALi::canonicalize_polynomial(p2);
        }
    if (p2.is_zero()){
        return CGAL::CGALi::canonicalize_polynomial(p1);
    }
    if (p2.degree() > p1.degree()) {
        Polynomial<NT> p3 = p1; p1 = p2; p2 = p3;
    }

    // compute gcd of content
    NT p1c = p1.content(), p2c = p2.content();
    NT gcdcont = CGALi::gcd(p1c,p2c);

    // compute gcd of primitive parts
    p1 /= p1c; p2 /= p2c;

    NT dummy;
    Polynomial<NT> q, r;

    NT g = NT(1), h = NT(1);
    for (;;) { 
        Polynomial<NT>::pseudo_division(p1, p2, q, r, dummy);
        if (r.is_zero()) { break; }
        
        if (r.degree() == 0) { 
            return CGAL::CGALi::canonicalize_polynomial(Polynomial<NT>(gcdcont));
        }
        int delta = p1.degree() - p2.degree();
        p1 = p2;
        p2 = r / (g * ipower(h, delta));
        g = p1.lcoeff();
        // h = h^(1-delta) * g^delta
        CGAL::CGALi::hgdelta_update(h, g, delta);
    }

    p2 /= p2.content() * p2.unit_part();

    // combine both parts to proper gcd
    p2 *= gcdcont; 

    return CGAL::CGALi::canonicalize_polynomial(p2);
}

template <class NT>
inline
Polynomial<NT> gcd_Euclidean_ring(
        Polynomial<NT> p1, Polynomial<NT> p2
) {
//    std::cout<<" gcd_Field"<<std::endl;
    // handle trivial cases
    if (p1.is_zero())
        if (p2.is_zero()) return Polynomial<NT>(NT(1));
        else return p2 / p2.unit_part();
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


template <class NT>
inline
Polynomial<NT> gcd_for_decomposible_fields(
        Polynomial<NT> p1, Polynomial<NT> p2
) {
    // std::cout<<" gcd_decompos_field"<<std::endl;
    typedef Polynomial<NT> POLY;
    typedef typename CGAL::Fraction_traits<POLY>::Numerator_type INTPOLY;
    typedef typename CGAL::Fraction_traits<POLY>::Denominator_type DENOM;
    typedef typename INTPOLY::NT INTNT;
    typename CGAL::Fraction_traits<POLY>::Decompose decompose;
    typename CGAL::Fraction_traits<POLY>::Compose   compose;
    
    DENOM dummy;
    p1.simplify_coefficients();
    p2.simplify_coefficients();
    INTPOLY p1i, p2i; 
    decompose(p1,p1i, dummy);
    decompose(p2,p2i, dummy);

    INTPOLY d0 = CGALi::gcd_utcf(p1i, p2i);
    POLY d = compose(d0, DENOM(1));
    d /= d.unit_part();
    d.simplify_coefficients();
    return d;
}


template <class NT>
inline
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
        if (*it != NT(0)) cont = CGALi::gcd_utcf_(cont, *it);
    }
    
    return cont;
}


template <class NT>
inline 
Polynomial<NT> gcd_utcf_Integral_domain( Polynomial<NT> p1, Polynomial<NT> p2){
    // std::cout<<" gcd_utcf_Integral_domain"<<std::endl;
        
// handle trivial cases
    if (p1.is_zero())
        if (p2.is_zero()){
            return Polynomial<NT>(NT(1));
        }else{
            return CGAL::CGALi::canonicalize_polynomial(p2);
        }
    if (p2.is_zero()){
        return CGAL::CGALi::canonicalize_polynomial(p1);
    }

    if (p2.degree() > p1.degree()) {
        Polynomial<NT> p3 = p1; p1 = p2; p2 = p3;
    }
    
    // remove redundant scalar factors
    p1=canonicalize_polynomial(p1);
    p2=canonicalize_polynomial(p2); 

    // compute content of p1 and p2
    NT p1c = CGALi::content_utcf_(p1);
    NT p2c = CGALi::content_utcf_(p2);
    
    
    
    
    // compute gcd of content
    NT gcdcont = CGALi::gcd_utcf_(p1c, p2c);

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
        p2 = r / (g * ipower(h, delta));
        g = p1.lcoeff();
        // h = h^(1-delta) * g^delta
        CGAL::CGALi::hgdelta_update(h, g, delta);
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

    p2 = CGALi::div_utcf(p2, content_utcf_(p2), true);

    // combine both parts to proper gcd
    p2 *= gcdcont;

    Polynomial<NT> result; 
    
    // make poly unique
    result = CGAL::CGALi::canonicalize_polynomial(p2);
    return result; 
}


}  // namespace CGALi

}  // namespace CGAL

#endif //CGAL_POLYNOMIAL_GCD_IMPLEMENTATIONS_H
