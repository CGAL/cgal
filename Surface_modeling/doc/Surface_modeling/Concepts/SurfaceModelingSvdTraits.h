
/*!
\ingroup PkgSurfaceModelingConcepts
\cgalConcept

@brief Concept describing the set of requirements for SVD factorization.

\sa `SurfaceModelingSvdTraits::Matrix`
\sa `SurfaceModelingSvdTraits::Vector`

*/

class SurfaceModelingSvdTraits {
public:

/// \name Types 
/// @{

/*! 
<I>3x3</I> Matrix Type
*/ 
typedef Hidden_type Matrix; 

/*! 
<I>3x1</I>  Vector Type
*/ 
typedef Hidden_type Vector; 

/// @}

/// \name Creation
/// @{

/*!
Default constructor.
*/
SurfaceModelingSvdTraits();

/// @}

/// \name Operations
/// @{

/*!
Compute SVD factorization of matrix A. Singular values have to be sorted in decreasing order.
*/
bool compute(const Matrix& A);

/*!
Get matrix U 
*/
const Matrix& matrixU() const;

/*!
Get matrix V
*/
const Matrix& matrixV() const;

/// @}

}; /* end SparseLinearAlgebraTraitsWithPreFactor_d */


/*!
\ingroup PkgSurfaceModelingConcepts
\cgalConcept

`SurfaceModelingSvdTraits::Vector` is a concept of a vector that can be multiplied by a sparse matrix. 

\sa `SurfaceModelingSvdTraits`
\sa `SurfaceModelingSvdTraits::Matrix`
*/
class Vector {
public:
/// \name Creation 
/// @{

/*! 
Create a vector with three parameters. 
*/ 
Vector(double x, double y, double z); 

/// @} 

/// \name Operations 
/// @{

/*! 
Read/write access to a vector coefficient. 
*/ 
double operator[](int row) const; 

/*! 
Read/write access to a vector coefficient. 
*/ 
double& operator[](int row); 

/*! 
Multiply current vector with \f$ v^T \f$.
*/ 
Matrix mult(Vector v);

/// @}

}; /* end Vector */


/*!
\ingroup PkgSurfaceModelingConcepts
\cgalConcept

`SurfaceModelingSvdTraits::Matrix` is a concept of a 3x3 matrix. 

\sa `SurfaceModelingSvdTraits`
\sa `SurfaceModelingSvdTraits::Vector`
*/
class Matrix {
public:
/// \name Creation 
/// @{

/*! 
Default constructor.
*/ 
Matrix(); 

/// @} 

/// \name Operations 
/// @{

/*!
Assignment operator.
*/
Matrix& operator= (const Matrix & other);

/*!
Multiply all coefficients with `scalar`.
*/
Matrix operator* (double scalar) const;

/*!
Apply a coefficient based summation.
*/
Matrix& operator+= (const Matrix& other);

/*!
Multiply this with `other`.
*/
Matrix operator* (const Matrix& other) const;

/*! 
Read/write access to a matrix coefficient. 
*/ 
double get_coef(int row, int column) const; 

/*! 
Read/write access to a matrix coefficient. 
*/ 
double set_coef(int row, int column, double value); 

/*!
Get transpose.
*/
Matrix transpose() const;

/*!
Get determinant.
*/
double determinant() const;

/*! 
Set all coefficients to zero.
*/ 
void set_zero();
/// @}

}; /* end Matrix */




