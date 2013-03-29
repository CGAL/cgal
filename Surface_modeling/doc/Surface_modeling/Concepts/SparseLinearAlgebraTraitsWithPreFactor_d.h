
/*!
\ingroup PkgSurfaceModeling
\cgalConcept

@brief Concept describing the set of requirements for direct sparse linear system solver with factorization. 
The aim is using the same factorization to solve the system for different right-hand sides.

\cgalRefines `LinearAlgebraTraits_d` 

\cgalHasModel `CGAL::Eigen_solver_traits<T>`

\sa `SparseLinearAlgebraTraitsWithPreFactor_d::Matrix`
\sa `SparseLinearAlgebraTraitsWithPreFactor_d::Vector`

*/

class SparseLinearAlgebraTraitsWithPreFactor_d {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef Hidden_type Matrix; 

/*! 

*/ 
typedef Hidden_type Vector; 

/*! 

*/ 
typedef Hidden_type NT; 

/// @} 

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
Factorize the solver with sparse matrix A. 
This factorization is used in SparseLinearAlgebraTraitsWithPreFactor_d::linear_solver to solve the system for different right-hand sides.
@return true if the factorization is successful
*/ 
bool pre_factor(const Matrix& A, NT& D);

/*! 
Solve the sparse linear system <I>A\f$ \times \f$ X = B</I>, with factorization computed in SparseLinearAlgebraTraitsWithPreFactor_d::pre_factor.
@return true if the solver is successful
*/ 
bool linear_solver(const Vector& B, Vector& X);

/// @}

/*!
\ingroup PkgSurfaceModeling
\cgalConcept

Concept describing the set of requirements for vector that can be multiplied by a sparse matrix. 

\cgalRefines `LinearAlgebraTraits_d::Vector` 

\cgalHasModel `CGAL::Eigen_vector<T>`

\sa `SparseLinearAlgebraTraitsWithPreFactor_d`
\sa `SparseLinearAlgebraTraitsWithPreFactor_d::Matrix`

*/

class Vector {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef Hidden_type NT; 

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
\ingroup PkgSurfaceModeling
\cgalConcept

Concept describing the set of requirements for a sparse matrix class. 

\cgalRefines `LinearAlgebraTraits_d::Matrix` 

\cgalHasModel `CGAL::Eigen_sparse_matrix<T>`
\cgalHasModel `CGAL::Eigen_sparse_symmetric_matrix<T>`

\sa `SparseLinearAlgebraTraitsWithPreFactor_d`
\sa `SparseLinearAlgebraTraitsWithPreFactor_d::Vector`

*/

class Matrix {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef Hidden_type NT; 

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


}; /* end SparseLinearAlgebraTraitsWithPreFactor_d */

