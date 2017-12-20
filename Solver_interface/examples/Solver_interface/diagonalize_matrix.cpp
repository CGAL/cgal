#include <iostream>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_diagonalize_traits.h>
#else
#include <CGAL/Diagonalize_traits.h>
#endif

typedef double                                          FT;
typedef CGAL::cpp11::array<FT, 6>                       Eigen_matrix;
typedef CGAL::cpp11::array<FT, 3>                       Eigen_vector;
typedef CGAL::cpp11::array<FT, 9>                       Eigen_three_vectors;

// If Eigen is enabled, use it, otherwise fallback to the internal model
#ifdef CGAL_EIGEN3_ENABLED
typedef CGAL::Eigen_diagonalize_traits<FT, 3>           Diagonalize_traits;
#else
typedef CGAL::Diagonalize_traits<FT, 3>                 Diagonalize_traits;
#endif

int main(void)
{
  Eigen_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};

  // Fill matrix with random numbers
  for(std::size_t i=0; i<6; ++i)
    covariance[i] = rand();

  Eigen_vector eigenvalues;
  Eigen_three_vectors eigenvectors;

  if(!(Diagonalize_traits::diagonalize_selfadjoint_covariance_matrix(covariance,
                                                                     eigenvalues,
                                                                     eigenvectors)))
  {
    std::cerr << "Error: cannot diagonalize matrix" << std::endl;
    return -1;
  }

  // Print result
  for(std::size_t i=0; i<3; ++i)
  {
    std::cout << "Eigenvalue " << i+1 << " = " << eigenvalues[i] << std::endl
              << "  with eigenvector [ ";
    for(std::size_t j=0; j<3; ++j)
      std::cout << eigenvectors[3*i + j] << " ";
    std::cout << "]" << std::endl;
  }

  return 0;
}
