
// #define CGAL_RESULTANT_NUSE_MODULAR_ARITHMETIC 1

#include <CGAL/basic.h>

#include <CGAL/Timer.h>
CGAL::Timer timer1, timer2, timer3; 

#include <cassert>
#include <iostream>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Modular.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Random.h>

#include <CGAL/gen_sparse_polynomial.h>
#include <CGAL/new_resultant.h>


static CGAL::Random my_rnd(346); // some seed 

template<class Polynomial_d> 
void test_resultant(){
    
    typedef typename CGAL::Polynomial_traits_d<Polynomial_d> PT;
    typedef typename PT::Coefficient Coeff; 
    for (int k = 0; k < 1; k++){
        Polynomial_d F2 = 
            CGAL::generate_sparse_random_polynomial<Polynomial_d>(my_rnd,12); 
        Polynomial_d F1 = 
            CGAL::generate_sparse_random_polynomial<Polynomial_d>(my_rnd,12);
        Polynomial_d G1 = typename PT::Move()(F1,PT::d-1,0);
        Polynomial_d G2 = typename PT::Move()(F2,PT::d-1,0);
        
        CGAL::Timer timer_new;
        timer_new.start();
        Coeff F = CGAL::new_resultant(F1,F2,PT::d-1);
        // Coeff G = CGAL::new_resultant(G1,G2,0);
        timer_new.stop();

        assert(F==F);
        std::cout <<"timer1: "<< timer1.time() << std::endl;
        std::cout <<"timer2: "<< timer2.time() << std::endl;
        std::cout <<"timer3: "<< timer3.time() << std::endl;
        timer1.reset();
        timer2.reset();
        timer3.reset();
        std::cout <<"new time: "<< timer_new.time() << std::endl;
    }
    std::cout << "end resultant: " << PT::d << std::endl;
}
   

template<class Coeff> 
void test_multi_resultant(){
    
    typedef CGAL::Polynomial<Coeff>        Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;

    Polynomial_3 Q  =
        CGAL::generate_polynomial_degree_each<Polynomial_3>(my_rnd,4,100); 
    Polynomial_2 B1 =
        CGAL::generate_polynomial_degree_each<Polynomial_2>(my_rnd,6,80); 
    Polynomial_2 B2 = 
        CGAL::generate_polynomial_degree_each<Polynomial_2>(my_rnd,6,80); 
    std::cout << "Q:  "<< Q   << std::endl;
    std::cout << "B1: "<< B1 << std::endl;
    std::cout << "B2: "<< B2 <<std::endl;
    
    CGAL::Timer timer_new;
    timer_new.start();
    Polynomial_2 R1 = CGAL::new_resultant(Q,Polynomial_3(B1),0);
    Polynomial_1 R2 = CGAL::new_resultant(R1,B2,0);
    std::cout << std::endl;
    std::cout << R1.degree() << std::endl;
    std::cout << R2.degree() << std::endl;
    timer_new.stop();
    
    std::cout <<"timer1: "<< timer1.time() << std::endl;
    std::cout <<"timer2: "<< timer2.time() << std::endl;
    std::cout <<"timer3: "<< timer3.time() << std::endl;
    timer1.reset();
    timer2.reset();
    timer3.reset();

    std::cout <<"new time: "<< timer_new.time() << std::endl;
    std::cout << "end multi resultant: " << std::endl;
}
   

int main(){

    //CGAL::force_ieee_double_precision();
    CGAL::set_pretty_mode(std::cout);

    typedef CGAL::Arithmetic_kernel AK;
    typedef AK::Integer Integer; 
    typedef CGAL::Polynomial<Integer>      Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;

    test_multi_resultant<Integer>();
    
    // test_resultant<Polynomial_1>();
    // test_resultant<Polynomial_2>();
    //test_resultant<Polynomial_3>();
    

    typedef CGAL::Polynomial<CGAL::Modular>  MPolynomial_1;
    typedef CGAL::Polynomial<MPolynomial_1>  MPolynomial_2;
    typedef CGAL::Polynomial<MPolynomial_2>  MPolynomial_3;

    //test_resultant<MPolynomial_1>();
    //test_resultant<MPolynomial_2>();
    //test_resultant<MPolynomial_3>();

    typedef CGAL::Polynomial<AK::Rational>   RPolynomial_1;
    typedef CGAL::Polynomial<RPolynomial_1>  RPolynomial_2;
    typedef CGAL::Polynomial<RPolynomial_2>  RPolynomial_3;

    //test_resultant<RPolynomial_1>();
    //test_resultant<RPolynomial_2>();
    //test_resultant<RPolynomial_3>();
}

 

