/*!
\ingroup PkgSolverInterfaceConcepts
\cgalConcept

Concept describing the set of requirements for solving the normal equation \f$ A^t A X = A^t B \f$,
\f$ A \f$ being a matrix, \f$ At \f$ its transpose matrix, \f$ B \f$ and \f$ X \f$ being two vectors.

\sa `SparseLinearAlgebraTraits_d`

\cgalHasModelsBegin
\cgalHasModels{CGAL::Eigen_solver_traits<T>}
\cgalHasModelsEnd
*/
class NormalEquationSparseLinearAlgebraTraits_d
{
public:
/// \name Types
/// @{

/*!
Matrix type model of `SparseLinearAlgebraTraits_d::Matrix`
*/
typedef unspecified_type Matrix;

/*!
Vector type model of `SparseLinearAlgebraTraits_d::Vector`
*/
typedef unspecified_type Vector;

/*!
Number type
*/
typedef unspecified_type NT;

/// @}

/// \name Creation
/// @{

/*!
Default constructor.
*/
NormalEquationSparseLinearAlgebraTraits_d();

/// @}

/// \name Operations
/// @{

/*!
Factorize the sparse matrix `At * A`.
This factorization is used in `normal_equation_solver()` to solve the system for different right-hand side vectors.
@return `true` if the factorization is successful and `false` otherwise.
*/
bool normal_equation_factor(const Matrix& A);

/*!
Solve the sparse linear system `At * A * X = At * B`, with `A` being the matrix
provided in `#normal_equation_factor()`, and `At` its transpose matrix.
@return `true` if the solver is successful and `false` otherwise.
*/
bool normal_equation_solver(const Vector& B, Vector& X);

/*!
Equivalent to a call to \link normal_equation_factor() `normal_equation_factor(A)` \endlink
followed by a call to \link normal_equation_solver `normal_equation_solver(B,X)` \endlink .
*/
bool normal_equation_solver(const Matrix& A, const Vector& B, Vector& X);

/// @}

}; /* end NormalEquationSparseLinearAlgebraTraits_d */
