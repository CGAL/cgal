#include <iostream>
#include <CGAL/basic.h>
#include <cassert>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Residue.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

#include <CGAL/gen_sparse_polynomial.h>
#include <CGAL/Polynomial/Interpolator.h>

static CGAL::Random my_rnd(346); // some seed 

template<class Polynomial_d> 
void test_interpolator(){
    typedef typename CGAL::Polynomial_traits_d<Polynomial_d> PT;
    typedef typename PT::Innermost_coefficient_type IC; 
    typedef typename PT::Coefficient_type Coeff; 
    
    for(int k = 0 ; k < 5; k++){
      // gen some polynomial
      Polynomial_d F = 
            CGAL::generate_sparse_random_polynomial<Polynomial_d>
            (my_rnd, 20/PT::d);
        
        typedef CGAL::Polynomial_traits_d<Polynomial_d> PT;
        typedef std::pair<IC,Coeff> Point; 
        std::vector<Point> points;
        for(int i = 0; i <= F.degree(); i++){
            points.push_back(Point(IC(i),typename PT::Evaluate()(F,IC(i))));
        }

        
        typedef CGAL::internal::Interpolator<Polynomial_d> Interpolator_d;
        Interpolator_d interpolator_1(points.begin(),points.end());
        assert(F == interpolator_1.get_interpolant());
        
        Interpolator_d interpolator_2;
        assert(Polynomial_d(0) == interpolator_2.get_interpolant());
        for(unsigned int i = 0; i < points.size(); i++){
          interpolator_2.add_interpolation_point(points[i]);
          if(i%10 == 0 ) 
            assert(F != interpolator_2.get_interpolant());
        }
        assert(F == interpolator_2.get_interpolant());
    }
    //std::cout << "end interpolate: " << PT::d << std::endl;
}

#if CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

int main(){

    CGAL::set_pretty_mode(std::cout);

    typedef CGAL::Arithmetic_kernel AK;
    typedef AK::Integer Integer; 
    typedef CGAL::Polynomial<Integer>      Polynomial_1;
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
    typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;
     

    typedef CGAL::Sqrt_extension<Integer,Integer> EXT;
    typedef CGAL::Polynomial<EXT>            EPolynomial_1;
    typedef CGAL::Polynomial<EPolynomial_1>  EPolynomial_2;
    typedef CGAL::Polynomial<EPolynomial_2>  EPolynomial_3;

    test_interpolator<Polynomial_1>();
    test_interpolator<Polynomial_2>();
    test_interpolator<Polynomial_3>();
    test_interpolator<EPolynomial_1>();
    test_interpolator<EPolynomial_2>();
    test_interpolator<EPolynomial_3>();
    
    // Enforce IEEE double precision and to nearest for modular arithmetic
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);
    
    typedef CGAL::Polynomial<CGAL::Residue>  MPolynomial_1;
    typedef CGAL::Polynomial<MPolynomial_1>  MPolynomial_2;
    typedef CGAL::Polynomial<MPolynomial_2>  MPolynomial_3;
    test_interpolator<MPolynomial_1>();
    test_interpolator<MPolynomial_2>();
    test_interpolator<MPolynomial_3>();
}

#else

int main(){
  std::cout << " Test needs a default arithmetic kernel " << std::endl;  
  return 0; 
}

#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
