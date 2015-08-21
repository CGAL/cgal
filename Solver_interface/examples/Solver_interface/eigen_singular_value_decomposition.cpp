#include <CGAL/Simple_cartesian.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <CGAL/Eigen_svd.h>

typedef CGAL::Eigen_svd::FT     FT;
typedef CGAL::Eigen_svd::Vector Eigen_vector;
typedef CGAL::Eigen_svd::Matrix Eigen_matrix;

int main(void)
{
  std::size_t degree = 3;
  
  Eigen_vector B (degree);
  Eigen_matrix M (degree, degree);


  // Fill B and M with random numbers
  for (std::size_t i = 0; i < degree; ++ i)
    {
      B.set (i, rand());
      for (std::size_t j = 0; j < degree; ++ j)
	M.set (i, j, rand());
    }

  // Solve MX=B
  std::cout << CGAL::Eigen_svd::solve(M, B) << std::endl;

  // Print result
  std::cout << "Solution of SVD = [ ";
  for (std::size_t i = 0; i < degree; ++ i)
    std::cout << B.vector()[i] << " ";
  std::cout << "]" << std::endl;
  
  return 0;
}
