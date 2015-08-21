namespace CGAL {

/*!
\ingroup PkgSolver

The class `Eigen_vcm_traits` provides an interface to the diagonalization of Variance-Covariance Matrices of \ref thirdpartyEigen "Eigen".
The version 3.1 (or greater) of \ref thirdpartyEigen "Eigen" must be available on the system. 

\cgalModels `VCMTraits`


`T`: floating type

\sa http://eigen.tuxfamily.org

Example 
-------------- 


*/


  
template <typename FT>
class Eigen_vcm_traits{

public:
  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const cpp11::array<FT,6>& cov,
    cpp11::array<FT, 3>& eigenvalues);
  
  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const cpp11::array<FT,6>& cov,
    cpp11::array<FT, 3>& eigenvalues,
    cpp11::array<FT, 9>& eigenvectors);
  
  static bool
  extract_largest_eigenvector_of_covariance_matrix (
    const cpp11::array<FT,6>& cov,
    cpp11::array<FT,3> &normal);

};

} // namespace CGAL

