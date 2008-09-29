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

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>

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
        std::cerr << argv[0] << " degree coeff_size number_of_curves" 
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

    srand48(time(NULL));
    // Create random polynomial of given degree
    
    std::vector<Poly_int2> curves;

    for( int j=0; j<no; j++ ) {
        
        std::vector<Poly_int1> coeffs;
        for(int i=0;i<=deg;i++) {
            std::vector<Integer> curr_coeffs = randvector<Integer>(deg-i,cs);
            coeffs.push_back(Poly_int1(curr_coeffs.begin(),curr_coeffs.end()));
        }
        Poly_int2 f(coeffs.begin(),coeffs.end());
        curves.push_back(f);
    }
    
    CGAL_assertion(static_cast<int>(curves.size()) == no);

    CGAL::set_ascii_mode(std::cout);

    std::cout << no << std::endl;

    for(int j = 0 ; j < no; j++) {
        std::cout << curves[j] << std::endl;
    }
            
    return 0;    
}
