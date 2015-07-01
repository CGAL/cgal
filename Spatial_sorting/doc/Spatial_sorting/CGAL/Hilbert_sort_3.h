namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctionObjects

The function object `Hilbert_sort_3` sorts iterator ranges of 
`Traits::Point_3` along a Hilbert curve by recursively subdividing at the median 
or the middle depending on the `PolicyTag`. 

\tparam Traits must be a model for `SpatialSortingTraits_3`. 

*/
template< typename Traits, typename PolicyTag >
class Hilbert_sort_3 {
public:

/// \name Creation 
/// @{

/*!
constructs an instance with `traits` as traits class instance. 
*/ 
Hilbert_sort_3(const Traits &traits = Traits()); 

/// @} 

/// \name Operations 
/// @{

/*!
It sorts the range `[begin, end)`. 
\tparam RandomAccessIterator must be an iterator with value type `Traits::Point_3`. 
*/ 
template <class RandomAccessIterator> void operator() (RandomAccessIterator begin, RandomAccessIterator end) const; 

/// @}

}; /* end Hilbert_sort_3 */
} /* end namespace CGAL */
