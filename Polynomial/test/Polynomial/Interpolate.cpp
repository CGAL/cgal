#include <iostream>
#include <CGAL/basic.h>
#include <cassert>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Modular.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

static CGAL::Random my_rnd(346); // some seed 

#define CGAL_SNAP_CGALi_TRAITS_D(T)                             \
    typedef T PT;                                               \
    typedef typename PT::Polynomial_d          Polynomial_d;    \
    typedef typename PT::Coefficient           Coeff;           \
    typedef typename PT::Innermost_coefficient ICoeff;          \
    typedef CGAL::Polynomial_traits_d<Coeff> PTC;               \
    typedef CGAL::Exponent_vector Exponent_vector;              \
    typedef std::pair< CGAL::Exponent_vector , ICoeff > Monom;  \
    typedef std::vector< Monom > Monom_rep;      

template <class Polynomial_d_>
Polynomial_d_
generate_sparse_random_polynomial(int max_degree = 6){
    typedef CGAL::Polynomial_traits_d<Polynomial_d_> PT;
    CGAL_SNAP_CGALi_TRAITS_D(PT);
    typename PT::Construct_polynomial construct; 


    typedef CGAL::Exponent_vector Exponent_vector;
    typedef std::pair< CGAL::Exponent_vector , ICoeff > Monom;
    typedef std::vector< Monom > Monom_rep;

    int range = 20;
    int number_of_variables = PT::d;
    int number_of_coeffs = 
        CGAL::min(number_of_variables * (int)ceil(log(max_degree+1))+1,4);
    
    Polynomial_d result; 
    for(int i = 0; i < number_of_coeffs; i++){
        CGAL::Exponent_vector exps(PT::d);
        for(int j = 0; j < PT::d; j++){
            exps[j]=my_rnd.get_int(0,max_degree);
        }
        ICoeff c = ICoeff(my_rnd.get_int(-range,range));
        Monom_rep monom_rep;
        monom_rep.push_back(Monom(exps,c));
        result += construct(monom_rep.begin(), monom_rep.end());
    }
    
    return result;
}





template<class IC>
inline IC interpolate_resultant(CGAL::Polynomial<IC> F1, CGAL::Polynomial<IC> F2){
    typedef CGAL::Polynomial_traits_d<CGAL::Polynomial<IC> > PT;

    //std::cout << "F1 " << F1 << std::endl;
    //std::cout << "F2 " << F2 << std::endl;
    //std::cout  << typename PT::Resultant()(F1,F2) << std::endl;
    return typename PT::Resultant()(F1,F2);
}


template<class T>
inline
CGAL::Polynomial<T> interpolate_resultant(
        const CGAL::Polynomial<CGAL::Polynomial<T> >& F1, 
        const CGAL::Polynomial<CGAL::Polynomial<T> >& F2){

    typedef CGAL::Polynomial<CGAL::Polynomial<T> > POLY_d;
    typedef CGAL::Polynomial_traits_d<CGAL::Polynomial<CGAL::Polynomial<T> > > PT_d;
    typedef typename PT_d::Innermost_coefficient IC;
    
    typedef CGAL::Polynomial<IC> POLY_1; 
    
    CGAL::Polynomial<T> result_pq(0); 
    POLY_1 pq(1);
    int i = 0;
    CGAL::Timer evaluate_timer; 
    CGAL::Timer interpolate_resultant_timer; 
    CGAL::Timer interpolate_timer; 
    
    while(true){
        i++;
        POLY_1 p(-IC(i),IC(1));
        evaluate_timer.start();
        CGAL::Polynomial<T> F1_at_i = typename PT_d::Evaluate()(F1,i);
        CGAL::Polynomial<T> F2_at_i = typename PT_d::Evaluate()(F2,i);
        evaluate_timer.stop();

        interpolate_resultant_timer.start();
        CGAL::Polynomial<T> result_p= interpolate_resultant(F1_at_i, F2_at_i);
        interpolate_resultant_timer.stop();
        POLY_1 q(pq); 
        CGAL::Polynomial<T> result_q= result_pq;
        interpolate_timer.start();
        typename CGAL::Polynomial_traits_d<CGAL::Polynomial<T> >::Interpolate()(
                p, result_p,
                q, result_q,
                pq, result_pq);
        interpolate_timer.stop();
        if(result_pq == result_q) break; 
        else{
            q = pq; 
            result_q = result_pq;
        }
    }   
    return result_pq; 
}


template <class Coeff> 
inline
Coeff my_resultant(CGAL::Polynomial<Coeff> F1, CGAL::Polynomial<Coeff> F2){
    typedef CGAL::Polynomial<Coeff> POLY;
    typedef CGAL::Polynomial_traits_d<POLY> PT;
    POLY F1_ = typename PT::Move()(F1,PT::d-1, 0);
    POLY F2_ = typename PT::Move()(F2,PT::d-1, 0);
    return interpolate_resultant(F1_,F2_);
}



int main(){

    CGAL::set_pretty_mode(std::cout);

    typedef CGAL::CORE_arithmetic_kernel AK;
    typedef AK::Integer Integer; 
    typedef CGAL::Polynomial<Integer>      Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;
    Polynomial_3 F = generate_sparse_random_polynomial<Polynomial_3>();
    
    

    typedef CGAL::Polynomial<CGAL::Modular>  MPolynomial_1;
    typedef CGAL::Polynomial<MPolynomial_1>  MPolynomial_2;
    typedef CGAL::Polynomial<MPolynomial_2>  MPolynomial_3;
    
    {
        MPolynomial_2 F1 = CGAL::modular_image(generate_sparse_random_polynomial<Polynomial_2>(2));
        MPolynomial_2 F2 = CGAL::modular_image(generate_sparse_random_polynomial<Polynomial_2>(2));
        std::cout << "start gcd "<< std::endl;
        std::cout << F1 << std::endl;
        std::cout << F2 << std::endl;
        MPolynomial_2 G = CGAL::gcd(F1,F2);
        std::cout << "end gcd  "<< std::endl;
        F1  = CGAL::integral_division(F1,G);
        F2  = CGAL::integral_division(F2,G);
    
        //assert(CGAL::is_one(CGAL::gcd(F1,F2)));
//    std::cout << CGAL::Modular_traits<Polynomial_2>::Modular_image_inv()(my_resultant(F1,F2)) << std::endl;
        //    assert(CGAL::Polynomial_traits_d<MPolynomial_2>::Resultant()(F1,F2) == my_resultant(F1,F2));
    }
    {
        MPolynomial_3 F1 = CGAL::modular_image(generate_sparse_random_polynomial<Polynomial_3>(2));
        MPolynomial_3 F2 = CGAL::modular_image(generate_sparse_random_polynomial<Polynomial_3>(2));
        std::cout << "start gcd "<< std::endl;
        //std::cout << F1 << std::endl;
        //std::cout << F2 << std::endl;
        MPolynomial_3 G = CGAL::gcd(F1,F2);
        std::cout << "end gcd  "<< std::endl;
        F1  = CGAL::integral_division(F1,G);
        F2  = CGAL::integral_division(F2,G);
    
        assert(CGAL::is_one(CGAL::gcd(F1,F2)));
//    std::cout << CGAL::Modular_traits<Polynomial_3>::Modular_image_inv()(my_resultant(F1,F2)) << std::endl;
        assert(CGAL::Polynomial_traits_d<MPolynomial_3>::Resultant()(F1,F2) == my_resultant(F1,F2));
    }
    std::cout << " DONE " << std::endl;
}

