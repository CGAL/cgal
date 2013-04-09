
/*!
\ingroup PkgSurfaceModelingConcepts
\cgalConcept

@brief Concept describing the set of requirements for a direct sparse linear system solver with factorization.
A model of this concept provides an additional factorization method to solve the system for different right-hand vectors.

\cgalRefines `SparseLinearAlgebraTraits_d`

\cgalHasModel `CGAL::Eigen_solver_traits<T>`

*/

class SparseLinearAlgebraTraitsWithPreFactor_d {
public:

/// \name Creation
/// @{

/*!

Default constructor.

*/
SparseLinearAlgebraTraitsWithPreFactor_d();

/// @}

/// \name Operations
/// @{

/*!
Factorize the sparse matrix A.
This factorization is used in SparseLinearAlgebraTraitsWithPreFactor_d::linear_solver to solve the system for different right-hand side vectors.
@return true if the factorization is successful
*/
bool pre_factor(const Matrix& A, NT& D);

/*!
Solve the sparse linear system <I>A\f$ \times \f$ X = B</I>, with the matrix provided in SparseLinearAlgebraTraitsWithPreFactor_d::pre_factor.
@return true if the solver is successful
*/
bool linear_solver(const Vector& B, Vector& X);

/// @}

}; /* end SparseLinearAlgebraTraitsWithPreFactor_d */

