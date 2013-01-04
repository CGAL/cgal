namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctionObjects

The function object `Hilbert_sort_d` sorts iterator ranges of 
`Traits::Point_d` along a Hilbert curve by recursively subdividing at the median 
or the middle depending on the `PolicyTag`. 

### Requirements ###

`Traits` is a model for `SpatialSortingTraits_d`. 

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
sorts the range [`begin`, `end`). 
\cgalRequires `std::iterator_traits<RandomAccessIterator>::value_type` equals to `Traits::Point_d`. 
*/ 
template <class RandomAccessIterator> void operator() (RandomAccessIterator begin, RandomAccessIterator end) const; 

/// @}

}; /* end Hilbert_sort_d */
} /* end namespace CGAL */
