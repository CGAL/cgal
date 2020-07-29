
/*!
\ingroup PkgKernelDLinAlgConcepts
\cgalConcept

An instance of data type `Matrix` is a matrix of
variables of number type `NT`. The types `Matrix` and `Vector`
together realize many functions of basic linear algebra.

*/

class Matrix {
public:

/// \name Types
/// @{

/*!
the ring type of the components.
*/
typedef unspecified_type NT;

/*!
bidirectional iterator for accessing
all components row-wise.
*/
typedef unspecified_type iterator;

/*!
bidirectional iterator for accessing
all components row-wise.
*/
typedef unspecified_type const_iterator;


/*!
random access iterator for accessing row
entries.
*/
typedef unspecified_type row_iterator;

/*!
random access iterator for accessing row
entries.
*/
typedef unspecified_type const_row_iterator;

/*!
random access iterator for accessing
column entries.
*/
typedef unspecified_type column_iterator;

/*!
random access iterator for accessing
column entries.
*/
typedef unspecified_type const_column_iterator;

/*!
a tag class for identity initialization
*/
typedef unspecified_type Identity;

/*!
the vector type used.
*/
typedef unspecified_type Vector;

/// @}

/// \name Creation
/// @{

/*!
creates an instance `M` of type
`Matrix`.
*/
Matrix();

/*!
creates an instance `M` of type
`Matrix` of dimension \f$ n \times n\f$ initialized to the zero matrix.

*/
Matrix(int n);

/*!
creates an instance `M` of
type `Matrix` of dimension \f$ m \times n\f$ initialized to the zero
matrix.
*/
Matrix(int m, int n);

/*!
creates an instance
`M` of type `Matrix` of dimension
`p.first`\f$ \times\f$`p.second` initialized to the zero matrix.
*/
Matrix(std::pair<int,int> p);

/*!
creates an
instance `M` of type `Matrix` of dimension \f$ n \times n\f$
initialized to the identity matrix (times `x`).
*/
Matrix(int n, Identity, NT x = NT(1));

/*!
creates an instance `M`
of type `Matrix` of dimension \f$ m \times n\f$ initialized to the
matrix with `x` entries.
*/
Matrix(int m, int n, NT x);

/*!
creates an
instance `M` of type `Matrix`. Let \f$ S\f$ be the ordered set of
\f$ n\f$ column-vectors of common dimension \f$ m\f$ as given by the iterator
range `[first,last)`. `M` is initialized to an \f$ m \times n\f$
matrix with the columns as specified by \f$ S\f$.

\pre `Forward_iterator` has a value type `V` from which we require to provide a iterator type `V::const_iterator`, to have `V::value_type == NT`.

Note that `Vector` or `std::vector<NT>` fulfill these requirements.
*/
template <class Forward_iterator>
Matrix(Forward_iterator first, Forward_iterator last);

/*!
creates an instance
`M` of type `Matrix`. Let \f$ A\f$ be an array of \f$ n\f$
column-vectors of common dimension \f$ m\f$. `M` is initialized to an
\f$ m \times n\f$ matrix with the columns as specified by \f$ A\f$.
*/
Matrix(std::vector< Vector > A);

/// @}

/// \name Operations
/// @{

/*!
returns \f$ n\f$, the number of rows of
`M`.
*/
int row_dimension() ;

/*!
returns \f$ m\f$, the number of columns
of `M`.
*/
int column_dimension() ;

/*!
returns \f$ (m,n)\f$, the
dimension pair of `M`.
*/
std::pair<int,int> dimension() ;

/*!
returns the \f$ i\f$-th row of `M` (an
\f$ m\f$ - vector).

\pre \f$ 0 \le i \le m - 1\f$.
*/
Vector row(int i) ;

/*!
returns the \f$ i\f$-th column of `M`
(an \f$ n\f$ - vector).

\pre \f$ 0 \le i \le n - 1\f$.
*/
Vector column(int i) ;

/*!
returns \f$ M_{i,j}\f$.
\pre \f$ 0\le i\le m-1\f$ and \f$ 0\le j\le n-1\f$.
*/
NT& operator()(int i, int j) ;

/*!
swaps rows \f$ i\f$ and \f$ j\f$.
\pre \f$ 0\le i\le m-1\f$ and \f$ 0\le j\le m-1\f$.

*/
void swap_rows(int i, int j) ;

/*!
swaps columns \f$ i\f$ and
\f$ j\f$.

\pre \f$ 0\le i\le n-1\f$ and \f$ 0\le j\le n-1\f$.
*/
void swap_columns(int i, int j) ;

/*!
an iterator pointing to the
first entry of the \f$ i\f$th row.

\pre \f$ 0\le i\le m-1\f$.
*/
row_iterator row_begin(int i) ;

/*!
an iterator pointing beyond
the last entry of the \f$ i\f$th row.

\pre \f$ 0\le i\le m-1\f$.
*/
row_iterator row_end(int i) ;

/*!
an iterator pointing to the
first entry of the \f$ i\f$th row.

\pre \f$ 0\le i\le m-1\f$.
*/
const_row_iterator row_begin(int i) const;

/*!
an iterator pointing beyond
the last entry of the \f$ i\f$th row.

\pre \f$ 0\le i\le m-1\f$.
*/
const_row_iterator row_end(int i) const;

/*!
an iterator pointing
to the first entry of the \f$ i\f$th column.

\pre \f$ 0\le i\le n-1\f$.

*/
column_iterator column_begin(int i) ;

/*!
an iterator pointing
beyond the last entry of the \f$ i\f$th column.

\pre \f$ 0\le i\le n-1\f$.
*/
column_iterator column_end(int i) ;

/*!
an iterator pointing
to the first entry of the \f$ i\f$th column.

\pre \f$ 0\le i\le n-1\f$.

*/
const_column_iterator column_begin(int i) const;

/*!
an iterator pointing
beyond the last entry of the \f$ i\f$th column.

\pre \f$ 0\le i\le n-1\f$.
*/
const_column_iterator column_end(int i) const;


/*!
an iterator pointing to the first entry
of \f$ M\f$.
*/
iterator begin();

/*!
an iterator pointing beyond the last entry
of \f$ M\f$.
*/
terator end();

/*!
an iterator pointing to the first entry
of \f$ M\f$.
*/
const_iterator begin() const;

/*!
an iterator pointing beyond the last entry
of \f$ M\f$.
*/
const_terator end() const;


/*!
Test for equality.
*/
bool operator==(const Matrix& M1) ;

/*!
Test for inequality.
*/
bool operator!=(const Matrix& M1) ;

/// @}

/// \name Arithmetic Operators
/// @{

/*!
Addition.
\pre `M.row_dimension() == M1.row_dimension()`
\pre `M.column_dimension() == M1.column_dimension()`
*/
Matrix operator+ (const Matrix& M1);

/*!

Subtraction.
\pre `M.row_dimension() == M1.row_dimension()`
\pre `M.column_dimension() == M1.column_dimension()`
*/
Matrix operator- (const Matrix& M1);

/*!
Negation.
*/
Matrix operator-();

/*!
Multiplication.

\pre `M.column_dimension() = M1.row_dimension()`
*/
Matrix operator*(const Matrix& M1)
;

/*!
Multiplication with
vector.

\pre `M.column_dimension() = vec.dimension()`
*/
Vector operator*(const Vector& vec) ;

/*!
Multiplication of every entry with `x`.
*/
Matrix operator*(const NT& x, const Matrix& M);

/*!
Multiplication of every entry with `x`.
*/
Matrix operator*(const Matrix& M, const NT& x) ;

/// @}

}; /* end Matrix */

