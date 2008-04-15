#ifndef CGAL_NEW_RESULTANT_H
#define CGAL_NEW_RESULTANT_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Interpolator.h>
#include <CGAL/Polynomial/resultant.h>

CGAL_BEGIN_NAMESPACE

template <class Coeff> 
inline Coeff new_resultant( CGAL::Polynomial<Coeff>, CGAL::Polynomial<Coeff>, int);
namespace CGALi{

template <class Coeff> 
inline Coeff new_resultant_interpolate( CGAL::Polynomial<Coeff>, CGAL::Polynomial<Coeff> );
template <class Coeff> 
inline Coeff new_resultant_modularize( CGAL::Polynomial<Coeff>, CGAL::Polynomial<Coeff>, CGAL::Tag_true);
template <class Coeff> 
inline Coeff new_resultant_modularize( CGAL::Polynomial<Coeff>, CGAL::Polynomial<Coeff>, CGAL::Tag_false);
template <class Coeff> 
inline Coeff new_resultant_decompose( CGAL::Polynomial<Coeff>, CGAL::Polynomial<Coeff>, CGAL::Tag_true);
template <class Coeff> 
inline Coeff new_resultant_decompose( CGAL::Polynomial<Coeff>, CGAL::Polynomial<Coeff>, CGAL::Tag_false);
template <class Coeff> 
inline Coeff new_resultant_( CGAL::Polynomial<Coeff>, CGAL::Polynomial<Coeff>);

} // namespace CGALi

namespace CGALi{

template <class IC> 
inline IC 
new_resultant_interpolate( CGAL::Polynomial<IC> F, CGAL::Polynomial<IC> G){
    CGAL_precondition(CGAL::Polynomial_traits_d<CGAL::Polynomial<IC> >::d == 1);
    return CGALi::resultant(F,G); 
}

template <class Coeff_2> 
inline
CGAL::Polynomial<Coeff_2>  new_resultant_interpolate(
        CGAL::Polynomial<CGAL::Polynomial<Coeff_2> > F, 
        CGAL::Polynomial<CGAL::Polynomial<Coeff_2> > G, 
        int index = 0 ){
    
    typedef CGAL::Polynomial<Coeff_2> Coeff_1;
    typedef CGAL::Polynomial<Coeff_1> POLY;
    typedef CGAL::Polynomial_traits_d<POLY> PT;
    typedef typename PT::Innermost_coefficient IC; 

    CGAL_precondition(PT::d >= 2);

    CGAL_precondition(index >= 0);
    CGAL_precondition(index <  PT::d);
        
    POLY F_ = typename PT::Move()(F, index, 0);
    POLY G_ = typename PT::Move()(G, index, 0);
    CGAL_postcondition(index != 0 || F_ == F);
    
    typename PT::Degree degree; 
    int maxdegree = degree(F_,0)*degree(G_,PT::d-1) + degree(F_,PT::d-1)*degree(G_,0); 

    typedef std::pair<IC,Coeff_2> Point; 
    std::vector<Point> points; // interpolation points  

    typename CGAL::Polynomial_traits_d<Coeff_1>::Degree_vector degree_vector; 

    int i(0);
    CGAL::Exponent_vector ev_f(degree_vector(Coeff_1()));
    CGAL::Exponent_vector ev_g(degree_vector(Coeff_1()));
    while((int) points.size() <= maxdegree + 1){
        i++;
        Coeff_1 F_at_i = typename PT::Evaluate()(F_,IC(i));
        Coeff_1 G_at_i = typename PT::Evaluate()(G_,IC(i));
        if(degree_vector(F_at_i) >  ev_f || degree_vector(G_at_i) >  ev_g){
            points.clear();
            ev_f  = degree_vector(F_at_i);
            ev_g  = degree_vector(G_at_i);
            CGAL_postcondition(points.size() == 0);
        }
        if(degree_vector(F_at_i) ==  ev_f && degree_vector(G_at_i) ==  ev_g){
            Coeff_2 res_at_i = new_resultant(F_at_i, G_at_i, 0);
            points.push_back(Point(IC(i),res_at_i));
        }
    }
    Coeff_1 result = CGAL::Interpolator<Coeff_1>(points.begin(),points.end()).get_interpolant();
    
#ifndef CGAL_NDEBUG
    while((int) points.size() <= maxdegree + 3){
        i++;
        Coeff_1 F_at_i = typename PT::Evaluate()(F_,IC(i));
        Coeff_1 G_at_i = typename PT::Evaluate()(G_,IC(i));
        
        assert(degree_vector(F_at_i) <= ev_f);
        assert(degree_vector(G_at_i) <= ev_g);
        
        if(degree_vector(F_at_i) ==  ev_f && degree_vector(G_at_i) ==  ev_g){
            Coeff_2 res_at_i = new_resultant(F_at_i, G_at_i, 0);
            points.push_back(Point(IC(i), res_at_i));
        }
    }
    
    Coeff_1 result_ = CGAL::Interpolator<Coeff_1>(points.begin(), points.end()).get_interpolant();
    // the interpolate polynomial has to be stable !
    assert(result_ == result); 
#endif 
    return result; 
}


template <class Coeff> 
inline
Coeff new_resultant_modularize( 
        CGAL::Polynomial<Coeff> F, CGAL::Polynomial<Coeff> G, CGAL::Tag_false){
    return new_resultant_interpolate(F,G);
}

template <class Coeff> 
inline
Coeff new_resultant_modularize( 
        CGAL::Polynomial<Coeff> F, CGAL::Polynomial<Coeff> G, CGAL::Tag_true){
    return new_resultant_interpolate(F,G);
}


template <class Coeff> 
inline
Coeff new_resultant_decompose( 
        CGAL::Polynomial<Coeff> F, CGAL::Polynomial<Coeff> G, CGAL::Tag_false){
    typedef CGAL::Polynomial<Coeff> Polynomial; 
    typedef typename Modular_traits<Polynomial>::Is_modularizable Is_modularizable; 
    return new_resultant_modularize(F,G,Is_modularizable());
}

template <class Coeff> 
inline
Coeff new_resultant_decompose( 
        CGAL::Polynomial<Coeff> F, CGAL::Polynomial<Coeff> G, CGAL::Tag_true){  
    
    typedef Polynomial<Coeff> POLY;
    typedef typename Fraction_traits<POLY>::Numerator_type Numerator;
    typedef typename Fraction_traits<POLY>::Denominator_type Denominator;
    typename Fraction_traits<POLY>::Decompose decompose;
    typedef typename Numerator::NT RES;
    
    Denominator a, b;
    F.simplify_coefficients();
    G.simplify_coefficients();
    Numerator F0; decompose(F,F0,a);
    Numerator G0; decompose(G,G0,b);
    Denominator c = CGAL::ipower(a, G.degree()) * CGAL::ipower(b, F.degree());
    typedef typename Algebraic_structure_traits<RES>::Algebraic_category Algebraic_category;
    RES res0 =  CGAL::CGALi::new_resultant_(F0, G0);
    typename Fraction_traits<Coeff>::Compose comp_frac;
    Coeff res = comp_frac(res0, c);
    typename Algebraic_structure_traits<Coeff>::Simplify simplify;
    simplify(res);
    return res;
}


template <class Coeff> 
inline
Coeff new_resultant_( 
        CGAL::Polynomial<Coeff> F, CGAL::Polynomial<Coeff> G){
    typedef CGAL::Fraction_traits<Polynomial<Coeff > > FT;
    typedef typename FT::Is_fraction Is_fraction; 
    return new_resultant_decompose(F,G,Is_fraction());
}


} // namespace CGALi

template <class Coeff> 
inline
Coeff  new_resultant( 
        CGAL::Polynomial<Coeff> F, 
        CGAL::Polynomial<Coeff> G, 
        int index = CGAL::Polynomial_traits_d< CGAL::Polynomial<Coeff> >::d-1){
    typedef CGAL::Polynomial_traits_d<CGAL::Polynomial<Coeff> > PT;
    CGAL::Polynomial<Coeff> F_ = typename PT::Move()(F, index, 0);
    CGAL::Polynomial<Coeff> G_ = typename PT::Move()(G, index, 0);
    return CGALi::new_resultant_(F_,G_);
}
    
CGAL_END_NAMESPACE



#endif // CGAL_NEW_RESULTANT_H

