#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_matrix.h>

typedef CGAL::Eigen_solver_traits<> Eigen_solver;
typedef Eigen_solver::NT            FT;
typedef Eigen_solver::Matrix        Eigen_matrix;
typedef Eigen_solver::Vector        Eigen_vector;

int main(void)
{
  srand(static_cast<unsigned int>(time (NULL)));
  std::size_t degree = 3000;
  std::size_t nb_nonzero_coef = 100;

  Eigen_vector B(degree); // Zero vector
  Eigen_matrix A(degree);

  // Randomly make some coefficients of the matrix non-zero
  for(std::size_t i=0; i<nb_nonzero_coef; ++i)
  {
    int x = rand() % degree;
    int y = rand() % degree;

    FT value = rand() / (FT)RAND_MAX;

    A.add_coef(x, y, value);
    A.add_coef(y, x, value);
  }

  Eigen_vector X(degree);
  FT d;

  Eigen_solver solver;
  if(!(solver.linear_solver(A, B, X, d)))
  {
    std::cerr << "Error: linear solver failed" << std::endl;
    return -1;
  }

  std::cerr << "Linear solve succeeded" << std::endl;
  return 0;
}
