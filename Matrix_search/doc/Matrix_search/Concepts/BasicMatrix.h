
/*!
\ingroup PkgMatrixSearchConcepts
\cgalConcept


A class `BasicMatrix` has to provide the following
types and operations in order to be a model for
`BasicMatrix`.

\cgalHasModel `CGAL::Dynamic_matrix<M>`

\sa `MonotoneMatrixSearchTraits`
\sa `SortedMatrixSearchTraits`

*/

class BasicMatrix {
public:

/// \name Types
/// @{

/*!
The type of a matrix entry. It has to define
a copy constructor.
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
returns the
number of rows.
*/
int number_of_rows() const;

/*!
returns the entry at position (`row`, `column`).
\pre \f$ 0 \le\f$ `row` \f$ <\f$ `number_of_rows()`, and \f$ 0 \le\f$ `column` \f$ <\f$ `number_of_columns()`.
*/
Entry operator()( int row, int column) const;

/// @}

}; /* end BasicMatrix */

