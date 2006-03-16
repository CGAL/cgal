#include <CGAL/QP_solver/gmp_double.h>
#include <CGAL/Gmpz.h>
#include <CGAL/QP_solver/Double.h>
#include <CGAL/QP_solver.h>

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
 typedef CGAL::Tag_true Has_equalities_only_and_full_rank;
 typedef CGAL::Tag_false Is_in_standard_form;
};

typedef CGAL::QP_solver<Traits> Solver;


int main() {
 double c[] = {0, 3};
 double D_row_0[] = {1, 0};
 double D_row_1[] = {0, 0};
 double A_col_0[] = {1};
 double A_col_1[] = {2};
 double b[] = {1};

 Traits::Row_type rt[] = {Traits::EQUAL};
 bool fl[] = {true, true};
 bool fu[] = {false, true};
 double l[] = {0, 0};
 double u[] = {0, 1};

 double *rows_of_D[] = {D_row_0, D_row_1};
 double *cols_of_A[] = {A_col_0, A_col_1};

 Solver solver(2, 1, cols_of_A, b, c, rows_of_D, rt, fl, l, fu, u);

 if (solver.status() != Solver::INFEASIBLE) {
   std::cout << "Optimal feasible solution x: ";
   for (Solver::Variable_value_iterator it = solver.variables_value_begin(); it != solver.variables_value_end(); ++it)
     std::cout << *it << " ";
   std::cout << "f(x): " << solver.solution() << std::endl;
 }

 return 0;
} 
