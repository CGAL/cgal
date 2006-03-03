#include <CGAL/QP_solver/gmp_double.h>
#include <CGAL/Gmpz.h>
#include <CGAL/QP_solver/Double.h>
#include <CGAL/QP_solver.h>
#include <iostream>


struct Traits {
  enum Row_type {LESS_EQUAL, EQUAL, GREATER_EQUAL};
  typedef CGAL::Double ET;
  
  typedef double **A_iterator;
  typedef double *B_iterator;
  typedef double *C_iterator;
  typedef double **D_iterator;
  typedef bool *FL_iterator;
  typedef bool *FU_iterator;
  typedef double *L_iterator;
  typedef double *U_iterator;
  typedef Row_type *Row_type_iterator;

  typedef CGAL::Tag_false Is_linear;
  typedef CGAL::Tag_true Is_symmetric;
  typedef CGAL::Tag_false Has_equalities_only_and_full_rank;
  typedef CGAL::Tag_false Is_in_standard_form;
};

typedef CGAL::QP_solver<Traits> Solver;


int main() {

  double c[5] = {1, 1, 1, 1, 1};
  double D_row_0[5] = {0, 0, 0, 0, 0};
  double D_row_1[5] = {0, 0, 0, 0, 0};
  double D_row_2[5] = {0, 0, 0, 0, 0};
  double D_row_3[5] = {0, 0, 0, 0, 0};
  double D_row_4[5] = {0, 0, 0, 0, 0};

  double A_col_0[10] = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
  double A_col_1[10] = {1, 0, 0, 0, 1, 1, 1, 0, 0, 0};
  double A_col_2[10] = {0, 1, 0, 0, 1, 0, 0, 1, 1, 0};
  double A_col_3[10] = {0, 0, 1, 0, 0, 1, 0, 1, 0, 1};
  double A_col_4[10] = {0, 0, 0, 1, 0, 0, 1, 0, 1, 1};
  double b[10] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

  Traits::Row_type rt[10] = {Traits::GREATER_EQUAL, Traits::GREATER_EQUAL, Traits::GREATER_EQUAL, Traits::GREATER_EQUAL, Traits::GREATER_EQUAL,
			     Traits::GREATER_EQUAL, Traits::GREATER_EQUAL, Traits::GREATER_EQUAL, Traits::GREATER_EQUAL, Traits::GREATER_EQUAL};
  bool fl[5] = {false, true, false, false, false};
  bool fu[5] = {false, true, false, false, false};
  double l[5] = {0, 2, 0, 0, 0};
  double u[5] = {10000, 4, 10000, 10000, 10000};
  
  double* rows_of_D[5] = {D_row_0, D_row_1, D_row_2, D_row_3, D_row_4};
  double* cols_of_A[5] = {A_col_0, A_col_1, A_col_2, A_col_3, A_col_4};
  
  Solver solver(5, 10, cols_of_A, b, c, rows_of_D, rt, fl, l, fu, u);
  
  if (solver.status() != Solver::INFEASIBLE) {
    std::cout << "Basic variables x: ";
    for (Solver::Basic_variable_value_iterator it = solver.basic_original_variables_value_begin(); it != solver.basic_original_variables_value_end(); ++it)
      std::cout << to_double(*it) << " ";

    std::cout << std::endl << "Full variables x: ";
    for (Solver::Variable_value_iterator it = solver.variables_value_begin(); it != solver.variables_value_end(); ++it)
      std::cout << to_double(*it) << " ";

    std::cout << std::endl << "f(x): " << to_double(solver.solution()) << std::endl;
  }  
  else
    std::cout << "No feasible solution found" << std::endl;
  
return 0;

}
