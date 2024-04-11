/*!
\ingroup PkgSolverInterfaceConcepts
\cgalConcept

The concept `SparseLinearAlgebraTraits_d` is used to solve sparse linear systems <I>A\f$ \times \f$ X = B</I>.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Eigen_solver_traits<T>}
\cgalHasModelsEnd
*/
class SparseLinearAlgebraTraits_d
{
public:
/// \name Types
/// @{

/*!

*/
typedef unspecified_type Matrix;

/*!

*/
typedef unspecified_type Vector;

/*!

*/
typedef unspecified_type NT;

/// @}

/// \name Creation
/// @{

/*!
Default constructor.
*/
SparseLinearAlgebraTraits_d();

/// @}

/// \name Operations
/// @{

/*!
Solve the sparse linear system \f$A \times X = B\f$. Return `true` on success.
The solution is then \f$(1/D) \times X \f$.

\pre `A.row_dimension()` == `B.dimension()`
\pre `A.column_dimension()` == `X.dimension()`
*/
bool linear_solver(const Matrix& A, const Vector& B, Vector& X, NT& D);

/// @}
}; /* end SparseLinearAlgebraTraits_d */

/*!
\cgalConcept

`SparseLinearAlgebraTraits_d::Vector` is a concept of a vector that can be multiplied
by a sparse matrix.

\cgalRefines{DefaultConstructible}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Eigen_vector<T>}
\cgalHasModelsEnd

\sa `SparseLinearAlgebraTraits_d`
\sa `SparseLinearAlgebraTraits_d::Matrix`
*/
class SparseLinearAlgebraTraits_d::Vector
{
public:
/// \name Types
/// @{

/*!

*/
typedef unspecified_type NT;

/*!
Index type
*/
typedef unspecified_type Index;

/// @}

/// \name Creation
/// @{

/*!
Create a vector initialized with zeros.
*/
Vector(Index rows);

/*!
Copy constructor.
*/
Vector(const Vector& toCopy);
/// @}

/// \name Operations
/// @{

/*!
Return the vector's number of coefficients.
*/
Index dimension() const;

/*!
Read/write access to a vector coefficient.
\pre `0 <= row < dimension()`.
*/
NT operator[](Index row) const;

/*!

*/
NT& operator[](Index row);

/// @}

}; /* end Vector */

/*!
\cgalConcept

`SparseLinearAlgebraTraits_d::Matrix` is a concept of a sparse matrix class.

\cgalRefines{Assignable,DefaultConstructible}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Eigen_sparse_matrix<T>}
\cgalHasModels{CGAL::Eigen_sparse_symmetric_matrix<T>}
\cgalHasModelsEnd

\sa `SparseLinearAlgebraTraits_d`
\sa `SparseLinearAlgebraTraits_d::Vector`
*/
class SparseLinearAlgebraTraits_d::Matrix
{
public:
/// \name Types
/// @{

/*!
Index type
*/
typedef unspecified_type Index;

/*!

*/
typedef unspecified_type NT;

/// @}

/// \name Creation
/// @{

/*!
Create a square matrix initialized with zeros.
*/
Matrix(Index dimension);

/*!
Create a rectangular matrix initialized with zeros.
*/
Matrix(Index rows, Index columns);

/// @}

/// \name Operations
/// @{

/*!
Return the matrix number of rows.
*/
Index row_dimension() const;

/*!
Return the matrix number of columns.
*/
Index column_dimension() const;

/*!
Read access to a matrix coefficient.

\pre `0 <= row < row_dimension()`
\pre `0 <= column < column_dimension()`
*/
NT get_coef(Index row, Index column) const;

/*!
Write access to a matrix coefficient: `a_ij = a_ij + val`.

\pre `0 <= row < row_dimension()`
\pre `0 <= column < column_dimension()`
*/
void add_coef(Index row, Index column, NT value);

/*!
Write access to a matrix coefficient: `a_ij = val`.

Optimization: Users can indicate that the coefficient does not already exist
in the matrix by setting `new_coef` to `true`.

\pre `0 <= i < row_dimension()`
\pre `0 <= j < column_dimension()`
*/
void set_coef(Index row, Index column, NT value, bool new_coef = false);

/*!
Swaps the content of `*this` and `m`.
*/
void swap(Matrix& m);

/// Multiplication with a scalar.
friend Matrix
operator*(const NT& c, const Matrix& M);

/// Sum of two matrices.
friend Matrix
operator+(const Matrix& M0, const Matrix& M1);

/// @}

}; /* end Matrix */
