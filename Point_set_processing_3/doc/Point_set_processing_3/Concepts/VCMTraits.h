/*!
\ingroup PkgPointSetProcessingConcepts
\cgalConcept
Concept providing functions to extract eigen vectors and eigen value from
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
  /// fill `eigen_values` with the eigen values of the covariance matrix represented by `cov`.
  /// Eigen values are sorted by increasing order.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const cpp11::array<double,6>& cov,
    cpp11::array<double, 3>& eigen_values);

  /// Extract the eigenvector associated to the greatest eigen value
  /// of the covariance matrix represented by `cov`.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool
  extract_greatest_eigenvector_of_covariance_matrix (
    const cpp11::array<double,6>& cov,
    cpp11::array<double,3> &normal);
};

