 // TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: Algebraic_curve_kernel_2.C,v 1.2 2007/10/12 11:43:30 emeliyan Exp $
// 
//
// Author(s)     : Michael Kerber <mkerber@mpi-sb.mpg.de>
//
// ============================================================================

#include <CGAL/basic.h>

#include <CGAL/Benchmark/Benchmark.hpp>
#include <CGAL/Benchmark/Option_parser.hpp>


#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial/resultant.h>

template<class NT>
std::vector<NT> randvector(int degree,long bitsize) {
  std::vector<NT> coeffs(degree+1);
  for(int i=0;i<=degree;i++) {
    // Creating the coefficients VERY elementary...
    NT coeff=0;
    for(int j=0;j<bitsize-1;j++) {
      coeff = 2*coeff + (lrand48()%2);
    }
    // The last bit determines the sign
    if(lrand48()%2==0) {
      coeff=-coeff;
    }    
    coeffs[i]=coeff;
  }
  return coeffs;
}


int main( int argc, char** argv ) {

    if(argc <4) {
        std::cerr << argv[0] << " degree coeff_size number_of_surfaces" 
                  << std::endl;
        std::exit(1);
    }
    int deg = std::atoi(argv[1]);
    int cs  = std::atoi(argv[2]);
    int no  = std::atoi(argv[3]);

    typedef CGAL::CORE_arithmetic_kernel Arithmetic_kernel;
    typedef Arithmetic_kernel::Integer Integer;
    typedef CGAL::Polynomial<Integer> Poly_int1;
    typedef CGAL::Polynomial<Poly_int1> Poly_int2;
    typedef CGAL::Polynomial<Poly_int2> Poly_int3;

    srand48(time(NULL));
    // Create random polynomial of given degree
    
    std::vector<Poly_int3> surfaces;

    for( int j=0; j<no; j++ ) {
        
        std::vector<Poly_int2> bi_coeffs;
        for(int i=0;i<=deg;i++) {
            Poly_int2 curr_coeff;
            std::vector<Poly_int1> uni_coeffs;
            for(int k=0;k<=deg-i;k++) {
                std::vector<Integer> curr_inner_coeffs 
                    = randvector<Integer>(deg-i-k,cs);
                uni_coeffs.push_back
                    (Poly_int1(curr_inner_coeffs.begin(),
                               curr_inner_coeffs.end()));
            }
            bi_coeffs.push_back(Poly_int2(uni_coeffs.begin(),
                                          uni_coeffs.end()));
        }
        Poly_int3 f(bi_coeffs.begin(),bi_coeffs.end());
        surfaces.push_back(f);
    }
    
    CGAL_assertion(static_cast<int>(surfaces.size()) == no);

    CGAL::set_ascii_mode(std::cout);

    std::cout << no + (no*(no-1)/2) << std::endl;

    for(int i = 0 ; i < no; i++) {
        std::cout << CGAL::CGALi::resultant(surfaces[i],
                                            CGAL::diff(surfaces[i])) 
                  << std::endl;
    }

    for(int i = 0 ; i < no; i++) {
        for(int j = i+1; j < no; j++) {
            std::cout << CGAL::CGALi::resultant(surfaces[i],surfaces[j]) 
                      << std::endl;
        }
    }
            
    return 0;    
}
