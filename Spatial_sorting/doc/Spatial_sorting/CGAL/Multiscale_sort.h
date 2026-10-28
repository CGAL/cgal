namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctionObjects

The class `Multiscale_sort` represents a sorting algorithm adaptor.
Given a range of `n` points:
<OL>
<LI>it applies `Sort` on the last `(1 - ratio) * n` points,
<LI>it recurses on the first `ratio * n` points,
stopping when there are fewer than `threshold` points.
</OL>

\tparam Sort must be a function object type with a const call operator
that sorts a range of points given by two random access iterators.
The range may be empty. For example, `Hilbert_sort_2`,
`Hilbert_sort_3`, or `Hilbert_sort_d` can be used.

*/
template< typename Sort >
class Multiscale_sort {
public:

/// \name Creation
/// @{

/*!
constructs an instance with `traits` as traits class instance.
*/
Multiscale_sort (const Sort &sort = Sort(), std::ptrdiff_t threshold = 1, double ratio = 0.5);

/// @}

/// \name Operations
/// @{

/*!
sorts the range `[begin, end)`.
`Sort::operator()(InputIterator begin, InputIterator end)` must be defined.
*/
template <class InputIterator>
void operator() (InputIterator begin, InputIterator end) const;

/// @}

}; /* end Multiscale_sort */
} /* end namespace CGAL */
