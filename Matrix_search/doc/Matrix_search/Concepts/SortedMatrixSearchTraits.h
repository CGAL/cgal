
/*!
\ingroup PkgMatrixSearchConcepts
\cgalConcept

The concept `SortedMatrixSearchTraits` defines types and operations
needed to compute the smallest entry in a set of sorted matrices
that fulfills a certain feasibility criterion using the function
`CGAL::sorted_matrix_search`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Sorted_matrix_search_traits_adaptor<F,M>}
\cgalHasModelsEnd

\sa `CGAL::sorted_matrix_search()`
\sa `BasicMatrix`

*/

class SortedMatrixSearchTraits {
public:

/// \name Types
/// @{

/*!
The class used for representing matrices.
It has to be a model for `BasicMatrix`.
*/
typedef unspecified_type Matrix;

/*!
The class used for
representing the matrix elements.
*/
typedef Matrix::Value Value;

/*!
An adaptable binary function
class: `Value` \f$ \times\f$ `Value` \f$ \rightarrow\f$ `bool`
defining a non-reflexive total order on `Value`. This
determines the direction of the search.
*/
typedef unspecified_type Compare_strictly;

/*!
An adaptable binary function
class: `Value` \f$ \times\f$ `Value` \f$ \rightarrow\f$ `bool`
defining the reflexive total order on `Value` corresponding
to `Compare_strictly`.
*/
typedef unspecified_type Compare_non_strictly;

/// @}

/// \name Operations
/// @{

/*!
returns the `Compare_strictly` object to be used for
the search.
*/
Compare_strictly compare_strictly()
const;

/*!
returns the `Compare_non_strictly` object to be used
for the search.
*/
Compare_non_strictly compare_non_strictly()
const;

/*!
The
predicate to determine whether an element `a` is feasible.
It has to be monotone in the sense that `compare(a, b)` and
`is_feasible(a)` imply `is_feasible(b)`.
*/
bool is_feasible( const Value& a);

/// @}

}; /* end SortedMatrixSearchTraits */

