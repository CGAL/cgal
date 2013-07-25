
/*!
\ingroup PkgKernelDLinAlgConcepts
\cgalConcept

The data type `LinearAlgebraTraits_d` encapsulates two classes 
`Matrix`, `Vector` and many functions of basic linear algebra. 
An instance of data type `Matrix` is a matrix of variables of type 
`NT`. Accordingly, `Vector` implements vectors of variables of 
type `NT`. Most functions of linear algebra are <I>checkable</I>, 
i.e., the programs can be asked for a proof that their output is 
correct. For example, if the linear system solver declares a linear 
system \f$ A x = b\f$ unsolvable it also returns a vector \f$ c\f$ such that 
\f$ c^T A = 0\f$ and \f$ c^T b \neq 0\f$. 

\cgalHasModel `CGAL::Linear_algebraHd<RT>`
\cgalHasModel `CGAL::Linear_algebraCd<FT>`

*/

class LinearAlgebraTraits_d {
public:

/// \name Types 
/// @{

/*!
the number type of the components. 
*/ 
typedef unspecified_type NT; 

/*!
the vector type. 
*/ 
typedef unspecified_type Vector; 

/*!
the matrix type. 
*/ 
typedef unspecified_type Matrix; 

/// @} 

/// \name Operations 
/// @{

/*!
returns \f$ M^T\f$ 
(a `M.column_dimension()` \f$ \times\f$ `M.column_dimension()` - 
matrix). 
*/ 
static Matrix transpose(const Matrix& M); 

/*!
determines whether `M` has an inverse. It also 
computes either the inverse as \f$ (1/D) \cdot I\f$ or when no 
inverse exists, a vector \f$ c\f$ such that \f$ c^T \cdot M = 0 \f$. 
\pre \f$ M\f$ is square. 
*/ 
static bool inverse(const Matrix& M, Matrix& I, NT& D, 
Vector& c); 

/*!
returns the 
inverse matrix of `M`. More precisely, \f$ 1/D\f$ times the matrix 
returned is the inverse of `M`. 
\pre `determinant(M) != 0`.

\pre \f$ M\f$ is square. 
*/ 
static Matrix inverse(const Matrix& M, NT& D) ; 

/*!
returns the determinant \f$ D\f$ of 
`M` and sufficient information to verify that the value of the 
determinant is correct. If the determinant is zero then \f$ c\f$ is a 
vector such that \f$ c^T \cdot M = 0\f$. If the determinant is non-zero 
then \f$ L\f$ and \f$ U\f$ are lower and upper diagonal matrices respectively 
and \f$ q\f$ encodes a permutation matrix \f$ Q\f$ with \f$ Q(i,j) = 1\f$ iff \f$ i = 
q(j)\f$ such that \f$ L \cdot M \cdot Q = U\f$, \f$ L(0,0) = 1\f$, \f$ L(i,i) = U(i 
- 1,i - 1)\f$ for all \f$ i\f$, \f$ 1 \le i < n\f$, and \f$ D = s \cdot U(n - 1,n - 
1)\f$ where \f$ s\f$ is the determinant of \f$ Q\f$.

\pre `M` is square. 
*/ 
static NT determinant (const Matrix& M, Matrix& L, Matrix& 
U, std::vector<int>& q, Vector& c); 

/*!
verifies the conditions stated above. 
*/ 
static bool verify_determinant (const Matrix& M, NT D, 
Matrix& L, Matrix& U, const std::vector<int>& q, Vector& 
c); 

/*!
returns the 
determinant of `M`.

\pre `M` is square. 
*/ 
static NT determinant (const Matrix& M); 

/*!
returns 
the sign of the determinant of `M`.

\pre `M` is square. 
*/ 
static int sign_of_determinant (const Matrix& M); 

/*!
determines 
the complete solution space of the linear system \f$ M\cdot x = b\f$. If 
the system is unsolvable then \f$ c^T \cdot M = 0\f$ and \f$ c^T \cdot b 
\not= 0\f$. If the system is solvable then \f$ (1/D) x\f$ is a solution, 
and the columns of `spanning_vectors` are a maximal set of 
linearly independent solutions to the corresponding homogeneous 
system.

\pre `M.row_dimension() = b.dimension()`. 
*/ 
static bool linear_solver(const Matrix& M, const Vector& b, 
Vector& x, NT& D, Matrix& spanning_vectors, Vector& c); 

/*!
determines whether the linear 
system \f$ M\cdot x = b\f$ is solvable. If yes, then \f$ (1/D) x\f$ is a 
solution, if not then \f$ c^T \cdot M = 0\f$ and \f$ c^T \cdot b \not= 0\f$. 
\pre `M.row_dimension() = b.dimension()`. 
*/ 
static bool linear_solver(const Matrix& M, const Vector& b, 
Vector& x, NT& D, Vector& c) ; 

/*!
as above, but without the witness \f$ c\f$ 
\pre `M.row_dimension() = b.dimension()`. 
*/ 
static bool linear_solver(const Matrix& M, const Vector& b, 
Vector& x, NT& D) ; 

/*!
determines whether the system \f$ M \cdot x = b\f$ is solvable 

\pre `M.row_dimension() = b.dimension()`. 
*/ 
static bool is_solvable(const Matrix& M, const Vector& b) 
; 

/*!
determines whether the homogeneous linear system 
\f$ M\cdot x = 0\f$ has a non - trivial solution. If yes, then \f$ x\f$ is 
such a solution. 
*/ 
static bool homogeneous_linear_solver (const Matrix& M, 
Vector& x); 

/*!
determines the solution space of the 
homogeneous linear system \f$ M\cdot x = 0\f$. It returns the dimension 
of the solution space. Moreover the columns of `spanning_vecs` 
span the solution space. 
*/ 
static int homogeneous_linear_solver (const Matrix& M, 
Matrix& spanning_vecs); 

/*!
returns the indices of a maximal subset 
of independent columns of `M`. 
*/ 
static int independent_columns (const Matrix& M, 
std::vector<int>& columns); 

/*!
returns the rank of 
matrix `M` 
*/ 
static int rank (const Matrix & M); 

/// @}

}; /* end LinearAlgebraTraits_d */

