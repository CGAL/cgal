/*!
\ingroup PkgSolverConcepts
\cgalConcept
Concept providing functions to extract eigenvectors and eigenvalue from
covariance matrices represented by an array `a` as follows:
<center>
\f$ \begin{bmatrix}
  a[0] & a[1] & a[2] \\
  a[1] & a[3] & a[4] \\
  a[2] & a[4] & a[5] \\
 \end{bmatrix}\f$
</center>
\cgalHasModel `CGAL::Eigen_vcm_traits`
*/
class VCMTraits
{
public:
  /// fill `eigenvalues` with the eigenvalues of the covariance matrix represented by `cov`.
  /// Eigenvalues are sorted by increasing order.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const cpp11::array<double,6>& cov,
    cpp11::array<double, 3>& eigenvalues);

  /// Extract the eigenvector associated to the largest eigenvalue
  /// of the covariance matrix represented by `cov`.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool
  extract_largest_eigenvector_of_covariance_matrix (
    const cpp11::array<double,6>& cov,
    cpp11::array<double,3> &normal);
};

