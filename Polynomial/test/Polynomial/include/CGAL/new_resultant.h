#ifndef CGAL_NEW_RESULTANT_H
#define CGAL_NEW_RESULTANT_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Interpolator.h>
#include <CGAL/Polynomial/resultant.h>

CGAL_BEGIN_NAMESPACE

template <class Coeff> 
inline Coeff new_resultant( const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>&, int);
namespace CGALi{

template <class Coeff> 
inline Coeff new_resultant_interpolate( const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>& );
template <class Coeff> 
inline Coeff new_resultant_modularize( const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>&, CGAL::Tag_true);
template <class Coeff> 
inline Coeff new_resultant_modularize( const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>&, CGAL::Tag_false);
template <class Coeff> 
inline Coeff new_resultant_decompose( const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>&, CGAL::Tag_true);
template <class Coeff> 
inline Coeff new_resultant_decompose( const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>&, CGAL::Tag_false);
template <class Coeff> 
inline Coeff new_resultant_( const CGAL::Polynomial<Coeff>&, const CGAL::Polynomial<Coeff>&);

} // namespace CGALi

namespace CGALi{

template <class IC> 
inline IC 
new_resultant_interpolate( const CGAL::Polynomial<IC>& F, const CGAL::Polynomial<IC>& G){
    CGAL_precondition(CGAL::Polynomial_traits_d<CGAL::Polynomial<IC> >::d == 1);
    return CGALi::resultant(F,G); 
}

template <class Coeff_2> 
inline
CGAL::Polynomial<Coeff_2>  new_resultant_interpolate(
        const CGAL::Polynomial<CGAL::Polynomial<Coeff_2> >& F, 
        const CGAL::Polynomial<CGAL::Polynomial<Coeff_2> >& G){
    
    typedef CGAL::Polynomial<Coeff_2> Coeff_1;
    typedef CGAL::Polynomial<Coeff_1> POLY;
    typedef CGAL::Polynomial_traits_d<POLY> PT;
    typedef typename PT::Innermost_coefficient IC; 

    CGAL_precondition(PT::d >= 2);
    
    typename PT::Degree degree; 
    typename CGAL::Polynomial_traits_d<Coeff_1>::Degree_vector degree_vector; 

    int maxdegree = degree(F,0)*degree(G,PT::d-1) + degree(F,PT::d-1)*degree(G,0); 

    typedef std::pair<IC,Coeff_2> Point; 
    std::vector<Point> points; // interpolation points  
    
   
    int i(0);
    CGAL::Exponent_vector ev_f(degree_vector(Coeff_1()));
    CGAL::Exponent_vector ev_g(degree_vector(Coeff_1()));
    
   
    while((int) points.size() <= maxdegree + 1){
        i++;
        // timer1.start();
        Coeff_1 Fat_i = F.evaluate(Coeff_1(i));
        Coeff_1 Gat_i = G.evaluate(Coeff_1(i));
        // timer1.stop();
        
        if(degree_vector(Fat_i) >  ev_f || degree_vector(Gat_i) >  ev_g){
            points.clear();
            ev_f  = degree_vector(Fat_i);
            ev_g  = degree_vector(Gat_i);
            CGAL_postcondition(points.size() == 0);
        }
        if(degree_vector(Fat_i) ==  ev_f && degree_vector(Gat_i) ==  ev_g){
            // timer2.start();
            Coeff_2 res_at_i = new_resultant_interpolate(Fat_i, Gat_i);
            // timer2.stop();
            points.push_back(Point(IC(i),res_at_i));
        }      
    }
   
    // timer3.start();
    CGAL::Interpolator<Coeff_1> interpolator(points.begin(),points.end());
    Coeff_1 result = interpolator.get_interpolant();
    // timer3.stop();

#ifndef CGAL_NDEBUG
    while((int) points.size() <= maxdegree + 3){
        i++;
        Coeff_1 Fat_i = typename PT::Evaluate()(F,IC(i));
        Coeff_1 Gat_i = typename PT::Evaluate()(G,IC(i));
        
        assert(degree_vector(Fat_i) <= ev_f);
        assert(degree_vector(Gat_i) <= ev_g);
        
        if(degree_vector(Fat_i) ==  ev_f && degree_vector(Gat_i) ==  ev_g){
            Coeff_2 res_at_i = new_resultant(Fat_i, Gat_i, 0);
            points.push_back(Point(IC(i), res_at_i));
        }
    }
    CGAL::Interpolator<Coeff_1> interpolator_(points.begin(),points.end());
    Coeff_1 result_= interpolator_.get_interpolant();
    
     // the interpolate polynomial has to be stable !
    assert(result_ == result); 
#endif 
    return result; 
}


template <class Coeff> 
inline
Coeff new_resultant_modularize( 
        const CGAL::Polynomial<Coeff>& F, const CGAL::Polynomial<Coeff>& G, CGAL::Tag_false){
    return new_resultant_interpolate(F,G);
};

template <class Coeff> 
inline
Coeff new_resultant_modularize( 
        const CGAL::Polynomial<Coeff>& F, const CGAL::Polynomial<Coeff>& G, CGAL::Tag_true){
    
    typedef Polynomial_traits_d<CGAL::Polynomial<Coeff> > PT;
    typedef typename PT::Polynomial_d Polynomial;
    typedef typename PT::Innermost_coefficient IC;
    
    
    typedef Chinese_remainder_traits<Coeff> CRT;
    typedef typename CRT::Scalar_type Scalar;


    typedef typename CGAL::Modular_traits<Polynomial>::Modular_NT MPolynomial; 
    typedef typename CGAL::Modular_traits<Coeff>::Modular_NT      MCoeff; 
    typedef typename CGAL::Modular_traits<Scalar>::Modular_NT     MScalar;  
        
    typename CRT::Chinese_remainder chinese_remainder; 
    typename CGAL::Modular_traits<Coeff>::Modular_image_inv inv_map;


    typename PT::Degree_vector                                     degree_vector; 
    typename CGAL::Polynomial_traits_d<MPolynomial>::Degree_vector mdegree_vector;

    bool solved = false; 
    int prime_index = 0; 
    int n = 0;
    Scalar p,q,pq,s,t; 
    Coeff R, R_old; 
    
    CGAL::Timer timer_evaluate, timer_resultant, timer_cr; 
    
    do{
        MPolynomial mF, mG;
        MCoeff mR;
        timer_evaluate.start();
        do{
            // select a prime number
            int current_prime = -1;
            prime_index++;
            if(prime_index >= 2000){
                std::cerr<<"primes in the array exhausted"<<std::endl;
                assert(false);
                current_prime = CGALi::get_next_lower_prime(current_prime);
            } else{
                current_prime = CGALi::primes[prime_index];
            }
            CGAL::Modular::set_current_prime(current_prime);
            
            mF = CGAL::modular_image(F);
            mG = CGAL::modular_image(G);
            
        }while( degree_vector(F) != mdegree_vector(mF) || 
                degree_vector(G) != mdegree_vector(mG));
        timer_evaluate.stop();
        
        timer_resultant.start();
        n++;
        mR = new_resultant_interpolate(mF,mG);
        timer_resultant.stop();
        timer_cr.start();
        if(n == 1){ 
            // init chinese remainder
            q =  CGAL::Modular::get_current_prime(); // implicit ! 
            R = inv_map(mR);
        }else{
            // continue chinese remainder
            p = CGAL::Modular::get_current_prime(); // implicit!  
            R_old  = R ;
//            chinese_remainder(q,Gs ,p,inv_map(mG_),pq,Gs);             
//            cached_extended_euclidean_algorithm(q,p,s,t);
            CGALi_algorithm_M::Cached_extended_euclidean_algorithm
                <typename CRT::Scalar_type> ceea;
            ceea(q,p,s,t);
            pq =p*q;
            chinese_remainder(q,p,pq,s,t,R_old,inv_map(mR),R);
            q=pq;
        }
        solved = (R==R_old);
        timer_cr.stop();       
    } while(!solved);
        
    std::cout << "Time Evaluate   : " << timer_evaluate.time() << std::endl; 
    std::cout << "Time Resultant  : " << timer_resultant.time() << std::endl; 
    std::cout << "Time Chinese R  : " << timer_cr.time() << std::endl; 
    // CGAL_postcondition(R == new_resultant_interpolate(F,G));
    return R;
    // return new_resultant_interpolate(F,G);
}


template <class Coeff> 
inline
Coeff new_resultant_decompose( 
        const CGAL::Polynomial<Coeff>& F, const CGAL::Polynomial<Coeff>& G, CGAL::Tag_false){
#ifdef CGAL_RESULTANT_NUSE_MODULAR_ARITHMETIC
    return  new_resultant_modularize(F,G,CGAL::Tag_false());
#else
    typedef CGAL::Polynomial<Coeff> Polynomial; 
    typedef typename Modular_traits<Polynomial>::Is_modularizable Is_modularizable; 
    return new_resultant_modularize(F,G,Is_modularizable());
#endif
}

template <class Coeff> 
inline
Coeff new_resultant_decompose( 
        const CGAL::Polynomial<Coeff>& F, const CGAL::Polynomial<Coeff>& G, CGAL::Tag_true){  
    
    typedef Polynomial<Coeff> POLY;
    typedef typename Fraction_traits<POLY>::Numerator_type Numerator;
    typedef typename Fraction_traits<POLY>::Denominator_type Denominator;
    typename Fraction_traits<POLY>::Decompose decompose;
    typedef typename Numerator::NT RES;
    
    Denominator a, b;
    // F.simplify_coefficients(); not const 
    // G.simplify_coefficients(); not const 
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
        const CGAL::Polynomial<Coeff>& F, const CGAL::Polynomial<Coeff>& G){
    typedef CGAL::Fraction_traits<Polynomial<Coeff > > FT;
    typedef typename FT::Is_fraction Is_fraction; 
    return new_resultant_decompose(F,G,Is_fraction());
}


} // namespace CGALi

template <class Coeff> 
inline
Coeff  new_resultant( 
        const CGAL::Polynomial<Coeff>& F_, 
        const CGAL::Polynomial<Coeff>& G_, 
        int index = CGAL::Polynomial_traits_d< CGAL::Polynomial<Coeff> >::d-1){
    typedef CGAL::Polynomial_traits_d<CGAL::Polynomial<Coeff> > PT;
    CGAL::Polynomial<Coeff> F = typename PT::Move()(F_, index, 0);
    CGAL::Polynomial<Coeff> G = typename PT::Move()(G_, index, 0);
    return CGALi::new_resultant_(F,G);
}
    
CGAL_END_NAMESPACE



#endif // CGAL_NEW_RESULTANT_H

