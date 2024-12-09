
/*!
\ingroup PkgQPSolverConcepts
\cgalConcept

A model of `NonnegativeQuadraticProgram` describes a convex quadratic program of the form

\f{eqnarray*}{
\mbox{(QP)}& \mbox{minimize}
& \qpx^{T}D\qpx+\qpc^{T}\qpx+c_0 \\
&\mbox{subject to} & A\qpx\qprel \qpb, \\
& & \qpx \geq 0
\f}

in \f$ n\f$ real variables \f$ \qpx=(x_0,\ldots,x_{n-1})\f$.
Here,
<UL>
<LI>\f$ A\f$ is an \f$ m\times n\f$ matrix (the constraint matrix),
<LI>\f$ \qpb\f$ is an \f$ m\f$-dimensional vector (the right-hand side),
<LI>\f$ \qprel\f$ is an \f$ m\f$-dimensional vector of relations
from \f$ \{\leq, =, \geq\}\f$,

<LI>\f$ D\f$ is a symmetric positive-semidefinite \f$ n\times n\f$ matrix (the
quadratic objective function),

<LI>\f$ \qpc\f$ is an \f$ n\f$-dimensional vector (the linear objective
function), and
<LI>\f$ c_0\f$ is a constant.

</UL>

The description is given by appropriate <I>random-access</I>
iterators over the program data, see below. The program therefore
comes in <I>dense</I> representation which includes zero entries.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Quadratic_program<NT>}
\cgalHasModels{CGAL::Quadratic_program_from_mps<NT>}
\cgalHasModels{CGAL::Nonnegative_quadratic_program_from_iterators<A_it, B_it, R_it, FL_it, L_it, FU_it, U_it, D_it, C_it>}
\cgalHasModelsEnd

The value types of all iterator types (nested iterator
types, respectively, for `A_iterator` and `D_iterator`) must be
convertible to some common `IntegralDomain` `ET`.


\sa `QuadraticProgram`
\sa `LinearProgram`
\sa `NonnegativeLinearProgram`

*/

class NonnegativeQuadraticProgram {
public:

/// \name Types
/// @{

/*!
A random access iterator type to go
columnwise over the constraint matrix \f$ A\f$. The value type
is a random access iterator type for an individual column that
goes over the entries in that column.
*/
typedef unspecified_type A_iterator;

/*!
A random access iterator type to go over
the entries of the right-hand side \f$ \qpb\f$.
*/
typedef unspecified_type B_iterator;

/*!
A random access iterator type to go over the
relations \f$ \qprel\f$. The value type of `R_iterator` is
`CGAL::Comparison_result`.
*/
typedef unspecified_type R_iterator;

/*!
A random access iterator type to go rowwise
over the matrix \f$ 2D\f$. The value type
is a random access iterator type for an individual row that
goes over the entries in that row, up to (and including) the
entry on the main diagonal.
*/
typedef unspecified_type D_iterator;

/*!
A random access iterator type to go over the
entries of the linear objective function vector \f$ c\f$.
*/
typedef unspecified_type C_iterator;

/// @}

/// \name Operations
/// @{

/*!
returns the number \f$ n\f$ of variables (number
of columns of \f$ A\f$) in `qp`.
*/
int get_n() const;

/*!
returns the number \f$ m\f$ of constraints
(number of rows of \f$ A\f$) in `qp`.
*/
int get_m() const;

/*!
returns an iterator over the columns
of \f$ A\f$. The corresponding past-the-end iterator is `get_a()+get_n()`.
For \f$ j=0,\ldots,n-1\f$, `*(get_a()+j)` is a random access
iterator for column \f$ j\f$.
*/
A_iterator get_a() const;

/*!
returns an iterator over the entries
of \f$ \qpb\f$. The corresponding past-the-end iterator is
`get_b()+get_m()`.
*/
B_iterator get_b() const;

/*!
returns an iterator over the entries
of \f$ \qprel\f$. The corresponding past-the-end iterator is
`get_r()+get_m()`.
The value `CGAL::SMALLER` stands
for \f$ \leq\f$, `CGAL::EQUAL` stands for \f$ =\f$, and `CGAL::LARGER`
stands for \f$ \geq\f$.
*/
R_iterator get_r() const;

/*!
returns an iterator over the rows of
\f$ 2D\f$. The corresponding past-the-end iterator is `get_d()+get_n()`.
For \f$ i=0,\ldots,n-1\f$, `*(get_d()+i)` is a random access
iterator for the entries in row \f$ i\f$ <I>below or on the diagonal</I>.
The valid range of this iterator is
guaranteed to have length \f$ i+1\f$ but not more. Values to the right
of the diagonal are deduced from the symmetry requirement on \f$ D\f$.
*/
D_iterator get_d() const;

/*!
returns an iterator over the entries
of \f$ \qpc\f$. The corresponding past-the-end iterator is
`get_c()+get_n()`.
*/
C_iterator get_c() const;

/*!
returns the constant term \f$ c_0\f$ of the objective function.
*/
std::iterator_traits<C_iterator>::value_type get_c0() const;

/// @}

}; /* end NonnegativeQuadraticProgram */

