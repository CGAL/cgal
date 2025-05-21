
namespace CGAL {

/*!
\ingroup nt_util

Not all values of a type need to be valid. The function object
class `Is_valid` checks this.

For example, an expression like
`NT(0)/NT(0)` can result in an invalid number.
Routines may have as a precondition that all values are valid.

\cgalModels{AdaptableFunctor}

*/
template< typename T >
class Is_valid {
public:

/// \name Operations
/// @{

/*!
returns if the argument is valid.
*/
bool operator()(const T& x);

/// @}

}; /* end Is_valid */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup nt_util

The function object class `Max` returns the larger of two values.
The order is induced by the second template argument <TT>Less</TT>.
The default value for `Less` is `std::less`.

Note that `T` must be a model of `LessThanComparable`
in case `std::less` is used.

\cgalModels{AdaptableFunctor}

*/
template< typename T, typename Less >
struct Max {

/// \name Creation
/// @{

/*!
default constructor.
*/
Max();

/*!
The constructed object will use `c` to compare the arguments.
*/
Max(Less c);

/// @}

/// \name Operations
/// @{

/*!
returns the larger of `x` and `y`,
with respect to the order induced by `Less`.
*/
T operator()(const T& x, const T& y);

/// @}

}; /* end Max */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup nt_util

The function object class `Min` returns the smaller of two values.
The order is induced by the second template argument <TT>Less</TT>.
The default value for `Less` is `std::less`.

Note that `T` must be a model of `LessThanComparable`
in case `std::less` is used.

\cgalModels{AdaptableFunctor}

*/
template< typename T, typename Less >
struct Min {

/// \name Creation
/// @{

/*!
default constructor.
*/
Min();

/*!
The constructed object will use `c` to compare the arguments.
*/
Min(Less c);

/// @}

/// \name Operations
/// @{

/*!
returns the larger of `x` and `y`,
with respect to the order induced by `Less`.
*/
T operator()(const T& x, const T& y);

/// @}

}; /* end Min */
} /* end namespace CGAL */
