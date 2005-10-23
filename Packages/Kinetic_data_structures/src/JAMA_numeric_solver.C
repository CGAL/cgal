#include <CGAL/Polynomial/internal/numeric_solvers_support.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>


#include <TNT/tnt_array2d.h>
#include <TNT/tnt_array1d.h>
#include <TNT/jama_eig.h>

#include <algorithm>
#include <functional>
#include <iterator>

//#include <iomanip>


CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

static const double max_error_value =0.00005;





template <bool CLEAN, class NT> 
static void jama_compute_roots(const NT *begin, const NT *end,  NT lb, 
			       NT ub, std::vector<NT> &roots){
  int degree= end-begin-1;
  TNT::Array2D<NT> arr(degree, degree, 0.0);
  for (int i=0; i< degree; ++i){
    arr[0][i]=-begin[degree-i-1]/begin[degree];
  }
  for (int i=0; i+1< degree; ++i){
    arr[i+1][i]=1;
  }
    
  /*std::cout << begin << std::endl;
    std::cout << "[";
    for (int i=0; i< arr.dim1(); ++i){
    std::cout << "[";
    for (int j=0; j< arr.dim2(); ++j){
    std::cout << arr[i][j];
    if (j != arr.dim2()-1) std::cout << ',';
    }
    std::cout <<"]," << std::endl;
    }
    std::cout << "]" << std::endl;*/
  //std::cout << arr << std::endl;

  JAMA::Eigenvalue<NT> ev(arr);
  TNT::Array1D<NT> real, imag;
  ev.getImagEigenvalues(imag);
  ev.getRealEigenvalues(real);
  CGAL_Polynomial_assertion(imag.dim1()== real.dim1());

  NT tol;
  if (CLEAN) tol=.00005;
  else tol=0;

  for (int i=0; i< real.dim1(); ++i){
    //std::cout << "Testing " << real[i] << "+" << imag[i] << "i" << std::endl;
    if (root_is_good(real[i], imag[i], lb-tol, ub)){
      roots.push_back(real[i]/*polish_root(begin, end, real[i])*/);
      //std::cout << "Good was " <<  real[i] << "+" <<  imag[i] << "i\n";
    } else {
      //std::cout << "Rejected " << real[i] << "+" << imag[i] << "i\n";
    }
  }
  std::sort(roots.begin(), roots.end(), std::greater<NT>());
  if (CLEAN) check_first_root(begin, end, lb, roots);
}
 
void jama_polynomial_compute_roots(const double *begin, const double *end,
				   double lb, double ub,
				   std::vector<double> &roots){
 int degree= end-begin-1;
  switch( degree) {
  case -1:
  case 0:
    break;
  case 1:
    compute_linear_roots(begin,end, lb, ub, roots);
    break;
  case 2:
    compute_quadratic_roots(begin, end, lb, ub, roots);
    break;
  default:
    jama_compute_roots<false>(begin, end, lb, ub, roots);
  }
}

void jama_polynomial_compute_cleaned_roots(const double *begin, const double *end, 
					   double lb, double ub, 
					   std::vector<double> &roots){
 int degree= end-begin-1;
  switch( degree) {
  case -1:
  case 0:
    break;
  case 1:
    compute_linear_cleaned_roots(begin,end, lb, ub, roots);
    break;
  case 2:
    compute_quadratic_cleaned_roots(begin, end, lb, ub, roots);
    break;
  default:
    jama_compute_roots<true>(begin, end, lb, ub, roots);
  }
}


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
