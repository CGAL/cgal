/*!
\ingroup PkgSolverInterfaceConcepts
\cgalConcept

@brief Concept describing the set of requirements for a direct sparse linear system solver with factorization.
A model of this concept stores the left-hand matrix (denoted \f$ A \f$) and provides an additional factorization
method to solve the system for different right-hand vectors.

\cgalRefines{SparseLinearAlgebraTraits_d}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Eigen_solver_traits<T>}
\cgalHasModelsEnd
*/
class SparseLinearAlgebraWithFactorTraits_d {
public:

/// \name Creation
/// @{

/*!
Default constructor.
*/
SparseLinearAlgebraWithFactorTraits_d();

/// @}

/// \name Operations
/// @{

/// Factorize the sparse matrix `A`.
/// This factorization is used in `SparseLinearAlgebraWithFactorTraits_d::linear_solver()`
/// to solve the system for different right-hand side vectors.
/// See `::SparseLinearAlgebraTraits_d::linear_solver()` for the description of `D`.
/// \return `true` if the factorization is successful and `false` otherwise.
bool factor(const Matrix& A, NT& D);

/// Solve the sparse linear system \f$ A \times X = B\f$, with \f$ A \f$ being the matrix
/// provided in `SparseLinearAlgebraWithFactorTraits_d::factor()`.
/// \return `true` if the solver is successful and `false` otherwise.
bool linear_solver(const Vector& B, Vector& X);

/// Solve the sparse linear system \f$ A \times X = B\f$, with \f$ A \f$ being the matrix
/// provided in `SparseLinearAlgebraWithFactorTraits_d::factor()`.
/// \return `true` if the solver is successful and `false` otherwise.
bool linear_solver(const Matrix& B, Vector& X);

/// @}

}; /* end SparseLinearAlgebraWithFactorTraits_d */
