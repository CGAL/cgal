
/*!
\ingroup PkgSolverConcepts

\cgalConcept

The concept `SparseLinearAlgebraTraits_d` is used to solve sparse linear systems <I>A\f$ \times \f$ X = B</I>. 


\cgalHasModel `CGAL::Eigen_solver_traits<T>`


*/

class SparseLinearAlgebraTraits_d {
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

Solve the sparse linear system <I>A\f$ \times \f$ X = B</I>. Return true on success. The solution is then (1/D) \f$ \times \f$ X. 

\pre `A.row_dimension()` == `B.dimension()`
\pre `A.column_dimension()` == `X.dimension()`

*/ 
bool linear_solver(const Matrix& A, const Vector& B, Vector& X, NT& D); 

/// @}

}; /* end SparseLinearAlgebraTraits_d */


/*!
\cgalConcept

`SparseLinearAlgebraTraits_d::Vector` is a concept of a vector that can be multiplied by a sparse matrix. 

\cgalHasModel `CGAL::Eigen_vector<T>`

\sa `SparseLinearAlgebraTraits_d`
\sa `SparseLinearAlgebraTraits_d::Matrix`

*/


class SparseLinearAlgebraTraits_d::Vector {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef unspecified_type NT; 

/// @} 

/// \name Creation 
/// @{

/*!

Create a vector initialized with zeros. 

*/ 
Vector(int rows); 

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
int dimension() const; 

/*!

Read/write access to a vector coefficient. 
\pre `0 <= row < dimension()`. 

*/ 
NT operator[](int row) const; 

/*!

*/ 
NT& operator[](int row); 

/// @}

}; /* end Vector */

/*!

\cgalConcept

`SparseLinearAlgebraTraits_d::Matrix` is a concept of a sparse matrix class. 

\cgalHasModel `CGAL::Eigen_sparse_matrix<T>`
\cgalHasModel `CGAL::Eigen_sparse_symmetric_matrix<T>`

\sa `SparseLinearAlgebraTraits_d`
\sa `SparseLinearAlgebraTraits_d::Vector`

*/

  class SparseLinearAlgebraTraits_d::Matrix {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef unspecified_type NT; 

/// @} 

/// \name Creation 
/// @{

/*!

Create a square matrix initialized with zeros. 

*/ 
Matrix(int dimension); 

/*!

Create a rectangular matrix initialized with zeros. 

*/ 
Matrix(int rows, int columns); 

/// @} 

/// \name Operations 
/// @{

/*!

Return the matrix number of rows. 

*/ 
int row_dimension() const; 

/*!

Return the matrix number of columns. 

*/ 
int column_dimension() const; 

/*!

Read access to a matrix coefficient. 

\pre `0 <= row < row_dimension()`
\pre `0 <= column < column_dimension()`

*/ 
NT get_coef(int row, int column) const; 

/*!

Write access to a matrix coefficient: `a_ij = a_ij + val`. 

\pre `0 <= row < row_dimension()`
\pre `0 <= column < column_dimension()`

*/ 
void add_coef(int row, int column, NT value); 

/*!

Write access to a matrix coefficient: `a_ij = val`. 

Optimization: Caller can optimize this call by setting `new_coef` to true if the coefficient does not already exist in the matrix. 

\pre `0 <= i < row_dimension()`
\pre `0 <= j < column_dimension()`

*/ 
void set_coef(int row, int column, NT value, bool new_coef = false); 

/// @}

}; /* end Matrix */

