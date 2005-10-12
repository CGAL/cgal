#include <CGAL/Polynomial/internal/numeric_solvers.h>


#include "tnt_array2d.h"
#include "tnt_array1d.h"
#include "jama_eig.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <vector>

#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/internal/Rational/Derivative.h>

#ifdef CGAL_POLYNOMIAL_USE_GSL
#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>
#endif

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

static const double max_error_value =0.0000005;

/*double polish_root(const double *begin, const double *end,
		   double rt) {
  std::cout << "Polished from " << rt << " to ";
  typedef POLYNOMIAL_NS::Polynomial<double> Fn;
  typedef POLYNOMIAL_NS::internal::Derivative<Fn> Diff;
  Diff dx;
  Fn f(begin, end);
  Fn d= dx(f);

  //int numit=0;
  for (int i=0; i< 10; ++i) {
    double delta= -f(rt)/d(rt);
    rt= rt+delta;
  }
  std::cout << rt << std::endl;
  return rt;
  }*/


template <class NT> 
static inline void check_first_root(const NT *begin, const NT *end, 
				    NT lb, std::vector<NT> &roots){
  // if we are not close to the current time, then we are fine
  if (roots.empty()) return;
  if (roots.back() > lb+ .0005) return;
  
  typedef CGAL_POLYNOMIAL_NS::Polynomial<NT> Fn;
  typedef CGAL_POLYNOMIAL_NS::internal::Derivative<Fn> Diff;
  Diff dx;
  Fn f(begin, end);
  Fn d= dx(f);

  // switch
  //while (sign(d(roots.back().representation()))== ZERO) d= dx_(d);
  while (d(roots.back())==0) d= dx(d);
    
  // switch
  //if (sign(d(roots.back().representation()))==POSITIVE){
  if (d(roots.back()) > 0){
    roots.pop_back();
  }
}




template <class NT> 
static inline bool is_real(NT r, NT i) {
  return abs(i)< .00005*abs(r);
}




template <class NT>
struct root_compare {
  bool operator()(NT a, NT b) const {
    return a > b;
  }
};





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
  assert(imag.dim1()== real.dim1());

  NT tol;
  if (CLEAN) tol=.00005;
  else tol=0;

  for (int i=0; i< real.dim1(); ++i){
    //std::cout << "Testing " << real[i] << "+" << imag[i] << "i" << std::endl;
    if (is_real(real[i], imag[i]) && real[i] < ub && real[i] > lb-tol){
      roots.push_back(real[i]/*polish_root(begin, end, real[i])*/);
    }
  }
  std::sort(roots.begin(), roots.end(), root_compare<NT>());
  if (CLEAN) check_first_root(begin, end, lb, roots);
}






template <bool CLEAN, class NT>
static inline void compute_quadratic_roots(const NT *begin, const NT *end,  NT lb, NT ub,
			     std::vector<NT> &roots){
  NT max_error=0;
  if (CLEAN) max_error=max_error_value;
  CGAL_Polynomial_assertion(begin[2] != 0);
    
  NT desc= begin[1]*begin[1]-4*begin[0]*begin[2];
  if (desc < 0) return;
    
  NT ur= (-begin[1]-sqrt(desc))/(2*begin[2]);
  NT lr= (-begin[1]+sqrt(desc))/(2*begin[2]);

  // drop even roots
  if (ur >lb-max_error && ur < ub){
    if (!CLEAN || sign(begin[2]) != POSITIVE){
      roots.push_back(ur);
      if (lr >lb-max_error && lr < ub){
	roots.push_back(lr);
      }
    }
  } else {
    if (lr > lb-max_error && lr <ub){
      if (!CLEAN || sign(begin[2]) != NEGATIVE){
	roots.push_back(lr);
      }
    }
  }
}





template <bool CLEAN, class NT>
static inline void compute_linear_roots(const NT *begin, const NT *, 
				 NT lb, NT ub, 
				 std::vector<NT> &roots){
  if (CLEAN &&  begin[1]>0 ) return;
  NT max_error=0;
  if (CLEAN) max_error=max_error_value;
  NT r= -to_double(begin[0]/begin[1]);
  if (r > lb-max_error && r < ub) {
    roots.push_back(r);
  }
}
  




#ifdef CGAL_POLYNOMIAL_USE_GSL

template <bool CLEAN>
static inline void gsl_compute_roots(const double *begin, const double *end, 
				     double lb, double ub, 
				     std::vector<double> &roots_){

  double max_error=0;
  if (CLEAN) max_error=max_error_value;

  gsl_poly_complex_workspace *workspace;
  int degree= end-begin-1;
  workspace= gsl_poly_complex_workspace_alloc(degree+1);
    
  //! I can't use auto_ptr because it is an array
  double *roots= new double[2*degree];
    
  /*for ( int i=0; i< 2*fn.degree(); ++i){
    roots[i]=0;
    }*/
    
  int ret = gsl_poly_complex_solve(begin, degree+1, workspace,
				   roots);
  roots_.reserve(ret);

  if (ret!=GSL_SUCCESS){
    std::cerr << "GSL solver did not converge on ";
    std::copy(begin, end, std::ostream_iterator<double>(std::cerr, " "));
    std::cerr << std::endl;
  }
    
  for ( int i=degree-1; i>=0; --i){
    double r= roots[2*i];
    double c= roots[2*i+1];
    if (std::abs(c) < .000005 && r > lb -max_error && r < ub) {
      roots_.push_back(r);
    }
  }

  std::sort(roots_.begin(), roots_.end(), std::greater<double>() );
  delete roots;

  gsl_poly_complex_workspace_free(workspace);

  if (CLEAN) check_first_root(begin,end, lb, roots_);
  return;
}





template <bool CLEAN>
static inline void gsl_compute_cubic_roots(const double *begin, const double *end, 
				    double lb, double ub, 
				    std::vector<double> &roots_){
  Polynomial_precondition(begin[3] != 0);
  double max_error=0;
  if (CLEAN) max_error=max_error_value;
  double r[3];
  int num_roots= gsl_poly_solve_cubic(begin[2]/begin[3], 
				      begin[1]/begin[3],
				      begin[0]/begin[3], &r[0],&r[1],&r[2]);
  roots_.reserve(num_roots);
  // I want reverse sorted roots
  for (int i=num_roots-1; i>=0; --i){
    if (r[i]>  lb- max_error && r[i] < ub) roots_.push_back(r[i]);
  }
  if (CLEAN) check_first_root(begin, end, lb, roots_);
}





template <bool CLEAN>
static inline void gsl_polynomial_compute_roots(const double *begin, const double *end, 
					 double lb, double ub, 
					 std::vector<double> &roots){
  int degree= end-begin-1;
  switch( degree) {
  case -1:
  case 0:
    break;
  case 1:
    compute_linear_roots<CLEAN>(begin,end, lb, ub, roots);
    break;
  case 2:
    compute_quadratic_roots<CLEAN>(begin, end, lb, ub, roots);
    break;
  case 3:
    gsl_compute_cubic_roots<CLEAN>(begin, end, lb, ub, roots);
    break;
  default:
    gsl_compute_roots<CLEAN>(begin, end, lb, ub, roots);
  }
}




static inline double gsl_evaluate_polynomial(const double *b, const double *e, double t) {
  if (b==e) return 0;
  /*double *d= new double[coefs.size()];
  for (unsigned int i=0; i< coefs.size(); ++i){
    d[i]= coefs[i];
    }*/
  double v= gsl_poly_eval(b, e-b, t);
  return v;
}


#endif






template <bool CLEAN, class NT>
static inline void jama_polynomial_compute_roots_template(const NT *begin, const NT *end,  NT lb, 
						   NT ub, std::vector<NT> &roots){
 int degree= end-begin-1;
  switch( degree) {
  case -1:
  case 0:
    break;
  case 1:
    compute_linear_roots<CLEAN>(begin,end, lb, ub, roots);
    break;
  case 2:
    compute_quadratic_roots<CLEAN>(begin, end, lb, ub, roots);
    break;
  default:
    jama_compute_roots<CLEAN>(begin, end, lb, ub, roots);
  }; 
}





/*template <bool CLEAN>
inline void jama_polynomial_compute_roots_bool(const double *begin, const double *end,  double lb, 
					       double ub, std::vector<double> &roots){
  jama_polynomial_compute_roots_template<CLEAN>(begin, end, lb, ub, roots);

  if (false) {
    POLYNOMIAL_NS::Interval_arithmetic_guard gd;
    POLYNOMIAL_NS::To_interval<double> ti;
    std::vector<POLYNOMIAL_NS::Interval_nt> iroots;
    try {
      std::vector<POLYNOMIAL_NS::Interval_nt> c(end-begin);
      for (int i=0; i<= end-begin; ++i){
	c[i]= ti(begin[i]);
      }
      jama_polynomial_compute_roots_template<CLEAN>(&c[0], &c[0]+(end-begin), 
						    POLYNOMIAL_NS::Interval_nt(ti(lb)), 
						    POLYNOMIAL_NS::Interval_nt(ti(ub)), 
						    iroots);
    } catch (...) {
      std::cerr << "Computation failed" << std::endl;
    }
    std::copy(iroots.begin(), iroots.end(), 
	      std::ostream_iterator<POLYNOMIAL_NS::Interval_nt>(std::cout, " "));
    std::cout << std::endl;
  }
  
  }*/
 
void jama_polynomial_compute_roots(const double *begin, const double *end,
				   double lb, double ub,
				   std::vector<double> &roots){
    jama_polynomial_compute_roots_template<false>(begin, end, lb, ub, roots);
}

void jama_polynomial_compute_cleaned_roots(const double *begin, const double *end, 
					   double lb, double ub, 
					   std::vector<double> &roots){
    jama_polynomial_compute_roots_template<true>(begin, end, lb, ub, roots);
}

void polynomial_compute_roots(const double *begin, const double *end,  double lb, 
			      double ub, std::vector<double> &roots){
#ifdef POLYNOMIAL_USE_GSL
    gsl_polynomial_compute_roots<false>(begin, end, lb, ub, roots);
#else
    jama_polynomial_compute_roots(begin, end, lb, ub, roots);
#endif
}

void polynomial_compute_cleaned_roots(const double *begin, const double *end,  double lb, 
				      double ub, std::vector<double> &roots){
#ifdef POLYNOMIAL_USE_GSL
    gsl_polynomial_compute_roots<true>(begin, end, lb, ub, roots);
#else
    jama_polynomial_compute_cleaned_roots(begin, end, lb, ub, roots);
#endif
}


double evaluate_polynomial(const double *b, const double *e, double t){
#ifdef POLYNOMIAL_USE_GSL
  return gsl_evaluate_polynomial(b, e, t);
#else
  if (b==e) return 0.0;

  const double *rit=e-1;
  double result = *rit;
  --rit;
  for (; rit != b-1; --rit){
    result *= t;
    result += (*rit);
  }
  return result;
#endif
}


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
