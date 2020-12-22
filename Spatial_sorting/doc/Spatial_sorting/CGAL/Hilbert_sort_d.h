namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctionObjects

The function object `Hilbert_sort_d` sorts iterator ranges of
`Traits::Point_d` along a Hilbert curve by recursively subdividing at the median
or the middle depending on the `PolicyTag`.

\tparam Traits must be a model for `SpatialSortingTraits_d`.

\tparam PolicyTag is used to specify the strategy policy.
Possible values are \link CGAL::Hilbert_sort_median_policy `Hilbert_sort_median_policy` \endlink
(the default policy) or \link CGAL::Hilbert_sort_middle_policy `Hilbert_sort_middle_policy` \endlink.

*/
template< typename Traits, typename PolicyTag >
class Hilbert_sort_d {
public:

/// \name Creation
/// @{

/*!
constructs an instance with `traits` as traits class instance.
*/
Hilbert_sort_d(const Traits &traits = Traits());

/// @}

/// \name Operations
/// @{

/*!
It sorts the range `[begin, end)`.
\tparam InputPointIterator must be a model of `RandomAccessIterator` with value type `Traits::Point_d`.
*/
template <class InputPointIterator>
void operator() (InputPointIterator begin, InputPointIterator end) const;

/// @}

}; /* end Hilbert_sort_d */
} /* end namespace CGAL */
