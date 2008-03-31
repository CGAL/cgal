#include <iostream>
#include <CGAL/basic.h>
#include <cassert>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial/ipower.h>
#include <CGAL/Random.h>

static CGAL::Random my_rnd(346); // some seed 

#define CGAL_SNAP_CGALi_TRAITS_D(T)                        \
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
generate_sparse_random_polynomial(int max_degree = 10){
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


int main(){

    typedef CGAL::Arithmetic_kernel AK;
    typedef AK::Integer Integer; 
    typedef CGAL::Polynomial<Integer>      Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;

    typedef CGAL::Polynomial<CGAL::Modular>  MPolynomial_1;
    typedef CGAL::Polynomial<MPolynomial_1>  MPolynomial_2;
    typedef CGAL::Polynomial<MPolynomial_2>  MPolynomial_3;

    Polynomial_3 F = generate_sparse_random_polynomial<Polynomial_3>();
    MPolynomial_3 mF = CGAL::Modular_traits<Polynomial_3>::Modular_image()(F);
    
    CGAL::Polynomial_traits_d<MPolynomial_3>::Degree mdegree_3; 
    
    for(int i = 0; i <= mdegree_3(mF);i++){}


    
    std::cout << F << std::endl;
    std::cout << mF << std::endl;
    
    
}

