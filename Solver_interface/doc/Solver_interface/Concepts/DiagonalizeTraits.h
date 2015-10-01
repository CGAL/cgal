/*!
\ingroup PkgSolverConcepts
\cgalConcept

Concept providing functions to extract eigenvectors and eigenvalues
from covariance matrices represented by an array `a`, using symmetric
diagonalization. For example, a matrix of dimension 3 is defined as
follows:
<center>
\f$ \begin{bmatrix}
  a[0] & a[1] & a[2] \\
  a[1] & a[3] & a[4] \\
  a[2] & a[4] & a[5] \\
 \end{bmatrix}\f$
</center>

\tparam FT Number type
\tparam dim Dimension of the matrices and vectors

\cgalHasModel `CGAL::Eigen_diagonalize_traits`
\cgalHasModel `CGAL::Diagonalize_traits`

*/

template <typename FT, unsigned int dim = 3>
class DiagonalizeTraits
{
public:

  typedef cpp11::array<FT, dim> Vector;
  typedef cpp11::array<FT, dim*dim> Matrix;
  typedef cpp11::array<FT, (dim * (dim+1) / 2)> Covariance_matrix;

  /// fill `eigenvalues` with the eigenvalues of the covariance matrix represented by `cov`.
  /// Eigenvalues are sorted by increasing order.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const Covariance_matrix& cov,
    Vector& eigenvalues);

  /// fill `eigenvalues` with the eigenvalues and `eigenvectors` with
  /// the eigenvectors of the covariance matrix represented by `cov`.
  /// Eigenvalues are sorted by increasing order.
  /// \return `true` if the operation was successful and `false`
  /// otherwise.
  static bool
  diagonalize_selfadjoint_covariance_matrix(
    const Covariance_matrix& cov,
    Vector& eigenvalues,
    Matrix& eigenvectors);


  /// Extract the eigenvector associated to the largest eigenvalue
  /// of the covariance matrix represented by `cov`.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool
  extract_largest_eigenvector_of_covariance_matrix (
    const Covariance_matrix& cov,
    Vector& normal);
};

