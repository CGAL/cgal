// This program shows how the QP-solver can be called from a C++
// program. PLEASE NOTE: the API of the QP-solver has been submitted 
// to the CGAL editorial board and is currently under revision. IT IS
// NOT GUARANTEED (ACTUALLY UNLIKELY) THAT THE CODE AS BELOW WILL STILL
// WORK IN FUTURE RELEASES. It will have to be adpated according to the
// final API of the package. Please contact Bernd Gaertner 
// (gaertner at inf.ethz.ch) in case of problems or feedback. 


// This example program solves the quadratic program 
//    minimize    x^T D x + c^T x
//    subject to  A x <rel> b
//                L <= x <= u
// where <rel> stands for a vector of relations from {<=, ==, >=}
//
// In the concrete example below, we have the following data:
//     D      c     A      b    L     U
// ---------------------------------------------------------------------
//     1 0    0 3   1 2    1    0 0   infinity 1
//     0 0
//
// In other words, we minimize x^2 + 3y subject to x + 2y = 1
// and x >= 0, 0 <= y <= 1. Let's see what this gives: if we
// substitute x with 1 - 2y in the objective function, we get
// 4y^2 - y + 1. The derivative is 8y-1, so this is minimzed
// for y = 1/8. This value (as well as x = 1 - 2(1/8) = 3/4)
// are within the bounds, so the solution should be x = 3/4
// and y = 1/8. 

#include <CGAL/Gmpz.h>
#include <CGAL/QP_solver.h>

// here we declare the types of the various QP entries
struct Traits {
  enum Row_type {LESS_EQUAL, EQUAL, GREATER_EQUAL}; 
  typedef Row_type *Row_type_iterator; // iterator for the constraint type

  // the exact internal type (input type must be convertible to this)
  typedef CGAL::Gmpz ET;

  // the input type is int, and the problem comes via arrays 
  typedef int **A_iterator;  // for the columns of A
  typedef int *B_iterator;   // for the right-hand side b
  typedef int *C_iterator;   // for the objective function c
  typedef int **D_iterator;  // for the rows of D
  typedef bool *FL_iterator;    // for finiteness of lower bound L
  typedef bool *FU_iterator;    // for finiteness of upper bound U
  typedef int *L_iterator;   // for lower bounds L
  typedef int *U_iterator;   // for upper bounds U

  typedef CGAL::Tag_false Is_linear;     // we have a proper QP
  typedef CGAL::Tag_true Is_symmetric;   // the matrix D is symmetric
  typedef CGAL::Tag_true Has_equalities_only_and_full_rank;
     // A has indeed full rank 1, and there are only equalities
  typedef CGAL::Tag_false Is_in_standard_form; // not only x >=0 bounds
};

typedef CGAL::QP_solver<Traits> Solver;  // the solver type


int main() {
  int c[]       = {0, 3};
  int D_row_0[] = {1, 0};
  int D_row_1[] = {0, 0};
  int A_col_0[] = {1};
  int A_col_1[] = {2};
  int b[]       = {1};

  Traits::Row_type rt[] = {Traits::EQUAL};
  bool fl[]        = {true, true};
  bool fu[]        = {false, true};  // x has upper bound infinity...
  int l[]       = {0, 0};
  int u[]       = {0, 1};         // ... and it's u-entry is ignored

  int *rows_of_D[] = {D_row_0, D_row_1};
  int *cols_of_A[] = {A_col_0, A_col_1};

  // now call the solver; the first two arguments are
  // the number of variables and the number of constraints (rows of A) 
  Solver solver(2, 1, cols_of_A, b, c, rows_of_D, rt, fl, l, fu, u);

  if (solver.status() == Solver::OPTIMAL) { // we know that, don't we?
    std::cout << "Optimal feasible solution x: ";
    for (Solver::Variable_value_iterator it = solver.variables_value_begin(); 
	 it != solver.variables_value_end(); ++it)
      std::cout << *it << " "; // variable values
    std::cout << "f(x): " << solver.solution() << std::endl;
  }

  return 0;
} 
