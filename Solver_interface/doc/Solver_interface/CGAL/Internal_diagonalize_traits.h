namespace CGAL {

/*!
\ingroup PkgSolver

The class `Internal_diagonalize_traits` provides an internal
implementation for the diagonalization of Variance-Covariance
Matrices.

\cgalModels `DiagonalizeTraits`

\tparam FT Number type
\tparam dim Dimension of the matrices and vectors

Example 
-------------- 

*/


  
template <typename FT, unsigned int dim = 3>
class Internal_diagonalize_traits{

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
  extract_smallest_eigenvector_of_covariance_matrix (
    const cpp11::array<FT, (dim * (dim+1) / 2)>& cov,
    cpp11::array<FT,dim> &normal);


};

} // namespace CGAL

