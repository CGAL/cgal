#include <CGAL/Accelerate_sparse_matrix.h>
#include <CGAL/Accelerate_vector.h>
#include <CGAL/Accelerate_solver_traits.h>

int main()
{
  {
    CGAL::Accelerate_sparse_matrix<double> A(3);

    A.set_coef(0,0,1);
    A.set_coef(0,1,2);
    A.set_coef(0,2,3);
    A.set_coef(1,1,5);
    A.set_coef(2,2,9);

    A.assemble_matrix();

    assert(A.get_coef(1,0) == 0);

    CGAL::Accelerate_vector<double> B(3), X(3);
    B[0] = 1.0;
    B[1] = 2.0;
    B[2] = 3.0;

    CGAL::Accelerate_solver_traits<double> ast;
    double D;
    ast.linear_solver(A, B, X, D);
    assert(D == 1.0);

  }
  std::cout << std::endl;
  {
    CGAL::Accelerate_sparse_matrix<double> A(3, true);// symmetric only set lower left

    A.set_coef(0,0,1);
    A.set_coef(1,0,2);
    A.set_coef(2,0,3);
    A.set_coef(1,1,5);
    A.set_coef(2,2,9);

    A.assemble_matrix();

    assert(A.get_coef(1,0) == A.get_coef(0, 1));
  }

  return 0;
}
