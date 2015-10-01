#ifndef CGAL_DEFAULT_DIAGONALIZE_TRAITS_H
#define CGAL_DEFAULT_DIAGONALIZE_TRAITS_H

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_diagonalize_traits.h>
#else
#include <CGAL/Diagonalize_traits.h>
#endif


namespace CGAL {


  // Wrapper class designed to automatically use
  // Eigen_diagonalize_traits if Eigen is available and otherwise use
  // the fallback Diagonalize_traits class.
  
template <typename FT, unsigned int dim = 3>
class Default_diagonalize_traits{

#ifdef CGAL_EIGEN3_ENABLED
  typedef Eigen_diagonalize_traits<FT, dim> Base;
#else
  typedef Diagonalize_traits<FT, dim> Base;
#endif

public:

  typedef cpp11::array<FT, dim> Vector;
  typedef cpp11::array<FT, dim*dim> Matrix;
  typedef cpp11::array<FT, (dim * (dim+1) / 2)> Covariance_matrix;
  
  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const Covariance_matrix& cov,
    Vector& eigenvalues)
  {
    return Base::diagonalize_selfadjoint_covariance_matrix (cov, eigenvalues);
  }

  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const Covariance_matrix& cov,
    Vector& eigenvalues,
    Matrix& eigenvectors)
  {
    return Base::diagonalize_selfadjoint_covariance_matrix (cov, eigenvalues, eigenvectors);
  }

  // Extract the eigenvector associated to the largest eigenvalue
  static bool
  extract_largest_eigenvector_of_covariance_matrix (
    const Covariance_matrix& cov,
    Vector& normal)
  {
    return Base::extract_largest_eigenvector_of_covariance_matrix (cov, normal);
  }
};


} // namespace CGAL

#endif // CGAL_DEFAULT_DIAGONALIZE_TRAITS_H
