namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctionObjects

The function object `Hilbert_sort_2` sorts iterator ranges of 
`Traits::Point_2` along a Hilbert curve by recursively subdividing 
at the median 
or the middle depending on the `PolicyTag`. 


\tparam Traits must be a model for `SpatialSortingTraits`. 

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
It sorts the range `[begin, end)`. 
\tparam RandomAccessIterator must be an iterator with value type `Traits::Point_2`. 
*/ 
template <class RandomAccessIterator> void operator() (RandomAccessIterator begin, RandomAccessIterator end) const; 

/// @}

}; /* end Hilbert_sort_2 */
} /* end namespace CGAL */
