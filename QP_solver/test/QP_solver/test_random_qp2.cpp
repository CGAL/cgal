#include <cstdlib>
#include <cassert>
#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Random.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_options.h>
#include <CGAL/QP_functions.h>

// choose exact floating point type
#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#else
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#endif

// program and solution types
typedef CGAL::Quadratic_program<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

// random number generator
CGAL::Random rd;

CGAL::Comparison_result random_rel()
{
  int z = rd.get_int(-1,2);
  return CGAL::Comparison_result(z);
}

void statistics (const Solution& s, 
		 unsigned int& o, unsigned int& i, unsigned int& u)
{
    switch (s.status()) {
    case CGAL::QP_OPTIMAL:
      o++;
      break;
    case CGAL::QP_INFEASIBLE:
      i++;
      break;
    case CGAL::QP_UNBOUNDED:
      u++;
      break;
    default:
      assert(false);
    }
}

unsigned int qp_optimal = 0;
unsigned int qp_infeasible = 0;
unsigned int qp_unbounded = 0;
unsigned int nqp_optimal = 0;
unsigned int nqp_infeasible = 0;
unsigned int nqp_unbounded = 0;
unsigned int lp_optimal = 0;
unsigned int lp_infeasible = 0;
unsigned int lp_unbounded = 0;
unsigned int nlp_optimal = 0;
unsigned int nlp_infeasible = 0;
unsigned int nlp_unbounded = 0;

// parameters
int tries = 5000;
int max_dim = 11; // must be >1
int max_entry = 11; // must be >0

int main() {  
  // print seed
  std::cout << "Random seed: " << rd.get_seed() << std::endl;

  // options
  CGAL::Quadratic_program_options options;
  options.set_auto_validation(true);

  // generate a set of small random qp's
  for (int i=0; i<tries; ++i) {
    // first choose dimensions
    int n = rd.get_int(1,max_dim);
    int m = rd.get_int(1,max_dim);

    // construct matrix D as C^T C, for C randomly chosen with n columns
    int k = rd.get_int (1, 2*n); // number of rows of C
    std::vector<std::vector<int> > C (k, std::vector<int>(n, 0));
    for (int j=0; j<k+n; ++j)  // sparse C
      C[rd.get_int(0, k)][rd.get_int(0,n)] = 
	rd.get_int(-max_entry, max_entry);

    // now fill the program 
    Program p;
    // A
    for (int j=0; j<n+m; ++j)
      p.set_a (rd.get_int(0,n), rd.get_int(0,m), rd.get_double());
    // b, r
    for (int i=0; i<m/2; ++i) {
      p.set_b (rd.get_int(0,m), rd.get_double());
      p.set_r (rd.get_int(0,m), random_rel());
    }
    // fl, l, fu, u
    for (int j=0; j<n; ++j) {
      double l = rd.get_double();
      double u = rd.get_double();
      if (l > u) std::swap (l, u); 
      p.set_l(j, rd.get_bool(), l);
      p.set_u(j, rd.get_bool(), u);
    }
    // D
    for (int i=0; i<n; ++i)
      for (int j=0; j<=i; ++j) {
	double entry = 0;
	for (int l=0; l<k; ++l) 
	  entry += C[l][i] * C[l][j];
	p.set_d(i, j, entry);
      }
    // c
    for (int j=0; j<n/2; ++j)
      p.set_c (rd.get_int(0, n), rd.get_double());
    // c0
    p.set_c0(rd.get_double());
    
    // solve it
    Solution s = CGAL::solve_quadratic_program (p, ET(), options);
    assert(s.is_valid());
    statistics (s, qp_optimal, qp_infeasible, qp_unbounded);

    // also solve it as nqp, lp, nlp
    s = CGAL::solve_nonnegative_quadratic_program (p, ET(), options); 
    assert(s.is_valid());
    statistics (s, nqp_optimal, nqp_infeasible, nqp_unbounded);
    s = CGAL::solve_linear_program (p, ET(), options);    
    assert(s.is_valid()); 
    statistics (s, lp_optimal, lp_infeasible, lp_unbounded);
    s = CGAL::solve_nonnegative_linear_program (p, ET(), options);   
    assert(s.is_valid());  
    statistics (s, nlp_optimal, nlp_infeasible, nlp_unbounded);
  }
  
  // output statistics
  std::cout << "Solved " << tries 
	    << " random QP / NQP  / LP / NLP .\n"
	    << " Optimal:    " 
	    << qp_optimal << " / " 
	    << nqp_optimal << " / " 
	    << lp_optimal << " / " 
	    << nlp_optimal << "\n"
	    << " Infeasible: "
	    << qp_infeasible << " / " 
	    << nqp_infeasible << " / " 
	    << lp_infeasible << " / " 
	    << nlp_infeasible << "\n"
	    << " Unbounded:  "
	    << qp_unbounded << " / " 
	    << nqp_unbounded << " / " 
	    << lp_unbounded << " / " 
	    << nlp_unbounded << std::endl;

  return 0;
}
