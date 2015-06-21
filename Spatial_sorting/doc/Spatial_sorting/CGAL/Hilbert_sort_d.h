namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctionObjects

The function object `Hilbert_sort_d` sorts iterator ranges of 
`Traits::Point_d` along a Hilbert curve by recursively subdividing at the median 
or the middle depending on the `PolicyTag`. 

\tparam Traits must be a model for `SpatialSortingTraits_d`. 

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
\tparam RandomAccessIterator must be an iteratoe with value type `Traits::Point_d`. 
*/ 
template <class RandomAccessIterator> void operator() (RandomAccessIterator begin, RandomAccessIterator end) const; 

/// @}

}; /* end Hilbert_sort_d */
} /* end namespace CGAL */
