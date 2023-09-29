namespace CGAL {

/*!
\ingroup PkgIntervalSkipListRef

The class `Interval_skip_list_interval` represents intervals with lower and upper
bound of type `Value`. These intervals
can be open or closed at each endpoint.

\cgalHeading{I/O}

The output operator is defined for `std::ostream`.

\cgalModels{Interval}

*/
template< typename Value >
class Interval_skip_list_interval {
public:

/// \name Creation
/// @{

/*!
%Default constructor.
*/
Interval_skip_list_interval();

/*!
Constructs the interval with infimum `i` and supremum `s`.
The arguments `ic` and `uc` have value `true`, iff
the interval is closed at the lower and upper bound, respectively.
*/
Interval_skip_list_interval(const Value& i,
const Value& s,
bool ic = true,
bool uc = true);

/// @}

/// \name Operations
/// @{

/*!
Returns `true`, iff the interval is closed at the lower bound.
*/
bool inf_closed() const;

/*!
Returns `true`, iff the interval is closed at the upper bound.
*/
bool sup_closed() const;

/// @}

}; /* end Interval_skip_list_interval */

/*!
Inserts the interval `i` into the stream `os`.

\pre The output operator for `Value` is defined.
\relates Interval_skip_list_interval
*/
template<typename V>
ostream& operator<<(ostream& os,
const Interval_skip_list_interval<V>& i);

} /* end namespace CGAL */
