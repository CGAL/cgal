#include <iostream>
#include <CGAL/basic.h>
#include <cassert>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Modular.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

#include <CGAL/gen_sparse_polynomial.h>
#include <CGAL/Interpolator.h>

static CGAL::Random my_rnd(346); // some seed 

template<class Polynomial_d> 
void test_interpolate(){
    typedef typename CGAL::Polynomial_traits_d<Polynomial_d> PT;
    typedef typename PT::Innermost_coefficient IC; 
    typedef typename PT::Coefficient Coeff; 
   
    for(int k = 0 ; k < 5; k++){
        Polynomial_d F = 
            CGAL::generate_sparse_random_polynomial<Polynomial_d>
            (my_rnd, 20/PT::d);
        
        typedef CGAL::Polynomial_traits_d<Polynomial_d> PT;
        typedef std::pair<IC,Coeff> Point; 
        std::vector<Point> points;
        for(int i = 0; i <= F.degree(); i++){
            points.push_back(Point(IC(i),typename PT::Evaluate()(F,IC(i))));
        }
        CGAL::Interpolator<Polynomial_d> 
            interpolator(points.begin(),points.end());
        Polynomial_d F_new = interpolator.get_interpolant();
        assert(F_new == F);
    }
    std::cout << "end interpolate: " << PT::d << std::endl;
}

int main(){

    CGAL::set_pretty_mode(std::cout);

    typedef CGAL::Arithmetic_kernel AK;
    typedef AK::Integer Integer; 
    typedef CGAL::Polynomial<Integer>      Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;
    

    typedef CGAL::Polynomial<CGAL::Modular>  MPolynomial_1;
    typedef CGAL::Polynomial<MPolynomial_1>  MPolynomial_2;
    typedef CGAL::Polynomial<MPolynomial_2>  MPolynomial_3;


    // test_resultant<Polynomial_3>();

    test_interpolate<Polynomial_1>();
    test_interpolate<Polynomial_2>();
    test_interpolate<Polynomial_3>();
    test_interpolate<MPolynomial_1>();
    test_interpolate<MPolynomial_2>();
    test_interpolate<MPolynomial_3>();
}

