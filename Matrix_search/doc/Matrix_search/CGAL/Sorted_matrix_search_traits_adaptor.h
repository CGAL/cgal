
namespace CGAL {

/*!
\ingroup PkgMatrixSearchRef

The class `Sorted_matrix_search_traits_adaptor` can be used
as an adaptor to create sorted matrix search traits classes for
arbitrary feasibility test and matrix classes `F` resp.\ `M`.

\cgalModels{SortedMatrixSearchTraits}

\tparam M must be a model for `BasicMatrix`
\tparam F must define a copy constructor and a monotone `bool operator()( const Value&)`.

*/
template< typename F, typename M >
class Sorted_matrix_search_traits_adaptor {
public:

/// \name Creation
/// @{

/*!
initializes `t` to use `m` for feasibility
testing.
*/
Sorted_matrix_search_traits_adaptor<F,M>( const F&
m);

/// @}

/// \name Types
/// @{

/*!
typedef to `M`.
*/
typedef M Matrix;

/*!
typedef to `Matrix::Value`.
*/
typedef Matrix::Value Value;

/*!
typedef to
`std::less<Value>`.
*/
typedef std::less<Value> Compare_strictly;

/*!
typedef to
`std::less_equal<Value>`.
*/
typedef std::less_equal<Value> Compare_non_strictly;

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
uses the
feasibility test given during creation.
*/
bool is_feasible(const Value& a);

/// @}

}; /* end Sorted_matrix_search_traits_adaptor */
} /* end namespace CGAL */
