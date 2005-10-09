#ifndef CGAL_POLYNOMIAL_INTERNAL_NUMERIC_SOLVERS_H
#define CGAL_POLYNOMIAL_INTERNAL_NUMERIC_SOLVERS_H
#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

void polynomial_compute_roots(const double *begin, const double *end,  
			      double lb, double ub,
			      std::vector<double> &roots);
void polynomial_compute_cleaned_roots(const double *begin, const double *end,  
				      double lb, double ub, 
				      std::vector<double> &roots);
double evaluate_polynomial(const double *b, const double *e, double t);


/*void compute_linear_roots(const double *begin, const double *, double lb, double ub, bool CLEAN, std::vector<double> &roots);
void compute_quadratic_roots(const double *begin, const double *,  double lb, double ub, bool CLEAN, 
			     std::vector<double> &roots);
			     void jama_compute_roots(const double *begin, const double *end,  double lb, double ub, bool CLEAN, std::vector<double> &roots);*/
void jama_polynomial_compute_roots(const double *begin, const double *end, 
				   double lb, double ub, std::vector<double> &roots);

void jama_polynomial_compute_cleaned_roots(const double *begin, const double *end, 
					   double lb, double ub, std::vector<double> &roots);
CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
