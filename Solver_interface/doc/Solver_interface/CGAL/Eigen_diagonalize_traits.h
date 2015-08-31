namespace CGAL {

/*!
\ingroup PkgSolver

The class `Eigen_diagonalize_traits` provides an interface to the
diagonalization of covariance matrices of \ref thirdpartyEigen
"Eigen".

The version 3.1 (or greater) of \ref thirdpartyEigen "Eigen" must be
available on the system.

\cgalModels `DiagonalizeTraits`


`FT`: floating type

`dim`: dimension of the matrices and vectors

\sa http://eigen.tuxfamily.org

Example 
-------------- 


*/


  
template <typename FT, unsigned int dim = 3>
class Eigen_diagonalize_traits{

public:
  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const cpp11::array<FT, (dim * (dim+1) / 2)>& cov,
    cpp11::array<FT, dim>& eigenvalues);

  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const cpp11::array<FT, (dim * (dim+1) / 2)>& cov,
    cpp11::array<FT, dim>& eigenvalues,
    cpp11::array<FT, dim * dim>& eigenvectors);

  static bool
  extract_largest_eigenvector_of_covariance_matrix (
    const cpp11::array<FT, (dim * (dim+1) / 2)>& cov,
    cpp11::array<FT,dim> &normal);


};

} // namespace CGAL

