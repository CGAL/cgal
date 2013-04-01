
#ifndef CGAL_GEN_SPARSE_POLYNOMIAL_H
#define CGAL_GEN_SPARSE_POLYNOMIAL_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Random.h>
#include <cmath>

namespace CGAL {

template <class Polynomial_d_>
Polynomial_d_
generate_sparse_random_polynomial(CGAL::Random random, int max_degree = 6){

    typedef CGAL::Polynomial_traits_d<Polynomial_d_> PT;
    typedef typename PT::Innermost_coefficient_type IC;
    typedef typename PT::Polynomial_d Polynomial_d;
    typename PT::Construct_polynomial construct; 

    typedef std::pair< CGAL::Exponent_vector , IC > Monom;
    typedef std::vector< Monom > Monom_rep;

    int range = 2000000;
    int number_of_variables = PT::d;
    double mdp = max_degree+1;
    int number_of_coeffs = 
        (CGAL::max)(number_of_variables * (int)std::ceil(std::log(mdp))+1,100);
    
    Polynomial_d result; 
    for(int i = 0; i < number_of_coeffs; i++){
        CGAL::Exponent_vector exps;
        for(int j = 0; j < PT::d; j++){
          exps.push_back(random.get_int(0,max_degree));
        }
        IC c = IC(random.get_int(-range,range));
        Monom_rep monom_rep;
        monom_rep.push_back(Monom(exps,c));
        result += construct(monom_rep.begin(), monom_rep.end());
        CGAL_postcondition(result.degree() >= 0);
    }
    // std::cout << result << std::endl;
    return result;
}



template <class Coeff> Coeff 
gen_randome_coeff_with_bits(CGAL::Random random, int bits = 5){
    Coeff c = 0; 
    int range = CGAL::ipower(2,30);
    for(int i = 0; i < (bits/30) ;i++){
        c *= range; 
        c += Coeff(random.get_int(0,range));
    }
    for(int i = 0; i < bits%30 ;i++){
        c += Coeff(random.get_int(0,2));
    }
    c *= CGAL::ipower(Coeff(-1),random.get_int(0,2));
    return c;
}

template <class Polynomial_d>
Polynomial_d
generate_polynomial_degree_each(CGAL::Random random, int max_degree = 6, int bits = 5){

    typedef CGAL::Polynomial_traits_d<Polynomial_d> PT;
    typedef typename PT::Innermost_coefficient_type IC;
    typename PT::Construct_polynomial construct; 

    typedef CGAL::Exponent_vector Exponent_vector;
    typedef std::pair< CGAL::Exponent_vector , IC > Monom;
    typedef std::vector< Monom > Monom_rep;


    
    int number_of_coeffs = CGAL::ipower(max_degree+1,PT::d);
    Monom_rep monom_rep;
   
    for(int i = 0; i < number_of_coeffs; i++){
        std::vector<int> v; 
        IC c = gen_randome_coeff_with_bits<IC>(random, bits);
        for(int j = 0; j < PT::d; j++){
            v.push_back((i/CGAL::ipower(max_degree+1,j))%(max_degree+1));
        }
        Exponent_vector ev(v.begin(),v.end());
        monom_rep.push_back(Monom(ev,c));
    }
    
    return construct(monom_rep.begin(),monom_rep.end());
}



} //namespace CGAL

#endif //  CGAL_GEN_SPARSE_POLYNOMIAL_H
