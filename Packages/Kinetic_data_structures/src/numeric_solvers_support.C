#include <CGAL/Polynomial/internal/numeric_solvers.h>

#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/internal/Rational/Derivative.h>



CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

static double max_error_value=.00005;

template <bool CLEAN, class NT>
static inline void compute_quadratic_roots_t(const NT *begin, const NT *end,  NT lb, NT ub,
					     std::vector<NT> &roots){
  NT max_error=0;
  if (CLEAN) max_error=max_error_value;
  CGAL_Polynomial_assertion(begin[2] != 0);
    
  NT desc= begin[1]*begin[1]-4*begin[0]*begin[2];
  if (desc <= 0) return;
    
  NT ur= (-begin[1]+sqrt(desc))/(2*begin[2]);
  NT lr= (-begin[1]-sqrt(desc))/(2*begin[2]);
  if (begin[2]< 0) std::swap(lr, ur);
  if (lr > lb-max_error && lr < ub) {
    roots.push_back(ur);
    if (lr > lb-max_error && lr < ub && (!CLEAN || /*lr > lb+max_error ||*/ begin[2] >0)){
      roots.push_back(lr);
    }
  } else {
    // only upper
    if (ur > lb-max_error && ur < ub && (!CLEAN || /*ur > lb+max_error ||*/ begin[2] <0)){
      roots.push_back(ur);
    }
  }
  
  // drop even roots
  /*if (ur >lb-max_error && ur < ub){
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
    }*/
}


void compute_quadratic_roots(const double *begin, const double *end, 
			     double lb, double ub, std::vector<double> &roots){
  return compute_quadratic_roots_t<false>(begin, end, lb, ub, roots);
}

void compute_quadratic_cleaned_roots(const double *begin, const double *end, 
				     double lb, double ub, std::vector<double> &roots){
  return compute_quadratic_roots_t<true>(begin, end, lb, ub, roots);
}


template <bool CLEAN, class NT>
static inline void compute_linear_roots_t(const NT *begin, const NT *, 
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


void compute_linear_roots(const double *begin, const double *end, 
			  double lb, double ub, std::vector<double> &roots){
  return compute_linear_roots_t<false>(begin, end, lb, ub, roots);
}

void compute_linear_cleaned_roots(const double *begin, const double *end, 
				  double lb, double ub, std::vector<double> &roots){
  return compute_linear_roots_t<true>(begin, end, lb, ub, roots);
}


template <class NT> 
static inline void check_first_root_t(const NT *begin, const NT *end, 
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

void check_first_root(const double *begin, const double *end,
		      double lb, std::vector<double> &roots){
  check_first_root_t(begin, end, lb, roots);
}

/*void polynomial_compute_roots(const double *begin, const double *end,  double lb, 
			      double ub, std::vector<double> &roots){
#ifdef POLYNOMIAL_USE_GSL
    gsl_polynomial_compute_roots(begin, end, lb, ub, roots);
#else
    jama_polynomial_compute_roots(begin, end, lb, ub, roots);
#endif
}

void polynomial_compute_cleaned_roots(const double *begin, const double *end,  double lb, 
				      double ub, std::vector<double> &roots){
#ifdef POLYNOMIAL_USE_GSL
    gsl_polynomial_compute_cleaned_roots(begin, end, lb, ub, roots);
#else
    jama_polynomial_compute_cleaned_roots(begin, end, lb, ub, roots);
#endif
}*/


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
