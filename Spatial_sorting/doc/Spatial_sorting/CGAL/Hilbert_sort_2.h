namespace CGAL {

/*!
\ingroup PkgSpatialSorting

The function object `Hilbert_sort_2` sorts iterator ranges of 
`Traits::Point_2` along a Hilbert curve by recursively subdividing 
at the median 
or the middle depending on the `PolicyTag`. 

### Requirements ###

`Traits` is a model for `SpatialSortingTraits`. 

*/
template< typename Traits, typename PolicyTag >
class Hilbert_sort_2 {
public:

/// \name Creation 
/// @{

/*! 
constructs an instance with `traits` as traits class instance. 
*/ 
Hilbert_sort_2(const Traits &traits = Traits()); 

/// @} 

/// \name Operations 
/// @{

/*! 
sorts the range [`begin`, `end`). 
\requires `RandomAccessIterator::value_type` equals to `Traits::Point_2`. 
*/ 
template <class RandomAccessIterator> void operator() (RandomAccessIterator begin, RandomAccessIterator end) const; 

/// @}

}; /* end Hilbert_sort_2 */
} /* end namespace CGAL */
