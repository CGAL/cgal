namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctionObjects

The function object `Hilbert_sort_on_sphere_3` sorts iterator ranges of 
`Traits::Point_3` along a Hilbert curve on the unit sphere. Actually, it approximates a Hilbert curve on
the sphere by a Hilbert curve on a cube with faces at \f$x, y, z = \pm 1/\sqrt{3}\f$. For each face of that cube, it calls an appropriate
version of `Hilbert_sort_2` which sorts a subset of the iterator range.
`Hilbert_sort_2` in each face is called with the median or the middle policy depending on the `PolicyTag`.

\tparam Traits must be a model for `SpatialSortingTraits_3`. 

*/
template< typename Traits, typename PolicyTag >
class Hilbert_sort_on_sphere_3 {
public:

/// \name Creation 
/// @{

/*!
constructs an instance with `traits` as traits class instance. 
*/ 
Hilbert_sort_on_sphere_3(const Traits &traits = Traits()); 

/// @} 

/// \name Operations 
/// @{

/*!
sorts the range `[begin, end)`. 
\tparam RandomAccessIterator must be an iterator with value type `Traits::Point_3`. 
*/ 
template <class RandomAccessIterator> void operator() (RandomAccessIterator begin, RandomAccessIterator end) const; 

/// @}

}; /* end Hilbert_sort_on_sphere_3 */
} /* end namespace CGAL */
