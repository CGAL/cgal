#include <set>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_matrix.h>


typedef CGAL::Eigen_solver_traits<Eigen::ConjugateGradient<CGAL::Eigen_sparse_matrix<double>::EigenType> > Eigen_solver;
typedef Eigen_solver::NT FT;
typedef Eigen_solver::Matrix Eigen_matrix;
typedef Eigen_solver::Vector Eigen_vector;


int main(void)
{
  srand (time (NULL));
  std::size_t degree = 3000;
  std::size_t nb_nonzero_coef = 100;
  
  Eigen_vector B (degree);
  Eigen_matrix A (degree);

  Eigen_vector diag (degree);
  for (std::size_t i = 0; i < nb_nonzero_coef; ++ i)
    {
      std::size_t x = rand () % degree;
      std::size_t y = rand () % degree;
      if (x == y)
	continue;

      FT value = rand () / static_cast<FT> (RAND_MAX);
	
      A.add_coef (x, y, value);
      diag.set (x, diag.vector()[x] - value);
      
      B.set (x, 1.);
    }

  for (std::size_t i = 0; i < degree; ++ i)
    A.add_coef (i, i, diag.vector()[i]);
  
  A.assemble_matrix();
  
  Eigen_vector X (degree);
  FT d;

  Eigen_solver solver;
  if (!(solver.linear_solver (A, B, X, d)))
    {
      std::cerr << "Error: linear solver failed" << std::endl;
      return -1;
    }
  
  // Print extract of result
  std::cout << "Vector X (non-zero coefficients) = [ ";
  for (std::size_t i = 0; i < degree; ++ i)
    if (X.vector()[i] != 0.)
      std::cout << X.vector()[i] << " ";
  std::cout << "]" << std::endl;
  
  return 0;
}
