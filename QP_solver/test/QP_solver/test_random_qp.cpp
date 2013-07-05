#include <cstdlib>
#include <cassert>
#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Random.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_options.h>
#include <CGAL/QP_functions.h>

// choose exact integral type
#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#else
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz ET;
#endif

// program and solution types
typedef CGAL::Quadratic_program_from_iterators
<int**,                                                // for A
 int*,                                                 // for b
 CGAL::Comparison_result*,                             // for r
 bool*,                                                // for fl
 int*,                                                 // for l
 bool*,                                                // for fu
 int*,                                                 // for u
 int**,                                                // for D
 int*>                                                 // for c 
Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

// random number generator
CGAL::Random rd;

// random entries
int random_unsigned() 
{
  bool z = rd.get_bool();
  if (z) 
    return 0;
  else
    return rd.get_int (0, 21);
}
int random_signed() 
{
  bool z = rd.get_bool();
  if (z) 
    return 0;
  else
    return rd.get_int (-10, 11);
}

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

int tries = 50000;

int main() {
  // print seed
  std::cout << "Random seed: " << rd.get_seed() << std::endl;

  // options
  CGAL::Quadratic_program_options options;
  options.set_auto_validation(true);

  // generate a set of small random qp's
  for (int i=0; i<tries; ++i) {
    int  Ax[] = {random_signed(), random_signed()};         
    int  Ay[] = {random_signed(), random_signed()};         
    int*  A[] = {Ax, Ay};                       
    int   b[] = {random_signed(), random_signed()};         
    CGAL::Comparison_result
      r[] = {random_rel(), random_rel()};
    bool fl[] = {rd.get_bool(), rd.get_bool()};                   
    int   l[] = {random_signed(),random_signed()};
    bool fu[] = {rd.get_bool(), rd.get_bool()};               
    int   u[] = {random_signed(),random_signed()};   
    // make sure that l<=u
    if (l[0] > u[0]) {int h = l[0]; l[0] = u[0]; u[0] = h;}
    if (l[1] > u[1]) {int h = l[1]; l[1] = u[1]; u[1] = h;}
    int  D1[] = {random_unsigned()};                    
    int  D2[] = {0, random_unsigned()};
    // can still change D_21 as long as D remains positive-semidefinite
    if (D1[0] < D2[1]) 
      D2[0] = rd.get_int(-D1[0], D1[0]+1);
    else
      D2[0] = rd.get_int(-D2[1], D2[1]+1);
    assert(D1[0] * D2[1] >= D2[0] * D2[0]);  
    int*  D[] = {D1, D2};                      
    int   c[] = {random_signed(), random_signed()};
    int  c0   = random_signed();                     

    // now construct the quadratic program; the first two parameters are
    // the number of variables and the number of constraints (rows of A)
    Program qp (2, 2, A, b, r, fl, l, fu, u, D, c, c0);

    // write/read it and check equality
    std::stringstream inout;
    CGAL::print_quadratic_program (inout, qp);
    CGAL::Quadratic_program_from_mps<int> qp2 (inout);
    assert(CGAL::QP_functions_detail::are_equal_qp (qp, qp2));
  
    // solve it
    Solution s = CGAL::solve_quadratic_program (qp, ET(), options);
    assert(s.is_valid());
    statistics (s, qp_optimal, qp_infeasible, qp_unbounded);

    // also solve it as nqp, lp, nlp
    s = CGAL::solve_nonnegative_quadratic_program (qp, ET(), options); 
    assert(s.is_valid());
    statistics (s, nqp_optimal, nqp_infeasible, nqp_unbounded);
    s = CGAL::solve_linear_program (qp, ET(), options);    
    assert(s.is_valid()); 
    statistics (s, lp_optimal, lp_infeasible, lp_unbounded);
    s = CGAL::solve_nonnegative_linear_program (qp, ET(), options);   
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
