
/*!
\ingroup PkgMatrixSearchConcepts
\cgalConcept


The concept `MonotoneMatrixSearchTraits` is a refinement of
`BasicMatrix` and defines types and operations needed to
compute the maxima for all rows of a totally monotone matrix using
the function `CGAL::monotone_matrix_search`.

\cgalHeading{Notes}

<UL>
<LI>For the sake of efficiency (and in order to achieve the time
bounds claimed for `monotone_matrix_search`), all these
operations have to be realized in constant time - except for
`extract_all_even_rows` which may take linear time.
<LI>There is an adaptor `Dynamic_matrix` that can be used to
add most of the functionality described above to arbitrary
matrix classes.
</UL>

\cgalHasModel `CGAL::Dynamic_matrix<M>`

\sa `CGAL::monotone_matrix_search()`

*/

class MonotoneMatrixSearchTraits {
public:

/// \name Types
/// @{

/*!
The type of a matrix entry.
*/
typedef unspecified_type Value;

/// @}

/// \name Operations
/// @{

/*!
returns the
number of columns.
*/
int number_of_columns() const;

/*!
returns the number
of rows.
*/
int number_of_rows() const;

/*!
returns the entry at position (`row`, `column`).
\pre \f$ 0 \le\f$ `row` \f$ <\f$ `number_of_rows()`, and \f$ 0 \le\f$ `column` \f$ <\f$ `number_of_columns()`.
*/
Entry operator()( int row, int column) const;

/*!
replace
column `old` with column number `new`. \pre \f$ 0 \le\f$ `old`, `new` \f$ <\f$ `number_of_columns()`.
*/
void replace_column( int old, int new);

/*!
returns
a new Matrix consisting of all rows of `m` with even index,
(i.e.\ first row is row \f$ 0\f$ of `m`, second row is row \f$ 2\f$ of
`m` etc.). \pre `number_of_rows()` \f$ > 0\f$.
*/
Matrix* extract_all_even_rows() const;

/*!
deletes the
rightmost columns, such that `m` becomes quadratic.
\pre `number_of_columns()` \f$ \ge\f$ `number_of_rows()`. \post `number_of_rows()` \f$ ==\f$ `number_of_columns()`.
*/
void shrink_to_quadratic_size();

/// @}

}; /* end MonotoneMatrixSearchTraits */

