namespace CGAL {

/*!
\ingroup PkgSolver

The class `Internal_vcm_traits` provides an internal implementation for the diagonalization of Variance-Covariance Matrices.

\cgalModels `VCMTraits`


`FT`: floating type

`dim`: dimension of the matrices and vectors

Example 
-------------- 

*/


  
template <typename FT, unsigned int dim = 3>
class Internal_vcm_traits{

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

