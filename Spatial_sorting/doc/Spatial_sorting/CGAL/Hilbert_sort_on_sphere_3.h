namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctionObjects

The function object `Hilbert_sort_on_sphere_3` sorts iterator ranges of 
`Traits::Point_3` along a Hilbert curve on a given sphere. Actually, it approximates a Hilbert curve on
that sphere by a Hilbert curve on a certain cube. For each face of that cube, it calls an appropriate
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
constructs an instance with `traits` as traits class instance, `sq_r` as the squared_radius of the given sphere,
and `p` as the center of the given sphere.
*/ 
Hilbert_sort_on_sphere_3(const Traits &traits = Traits(), 
                         double sq_r = 1.0, 
                         const Traits::Point_3 &p = Traits::Point_3(0,0,0)); 

/// @} 

/// \name Operations 
/// @{

/*!
sorts the range `[begin, end)` along a hilbert curve on the sphere centered at `p` with squared radius `sq_r`; 
these arguments are passed in the construction of the object `Hilbert_sort_on_sphere_3`. 
\tparam RandomAccessIterator must be an iterator with value type `Traits::Point_3`. 
*/ 
template <class RandomAccessIterator> void operator() (RandomAccessIterator begin, RandomAccessIterator end) const; 

/// @}

}; /* end Hilbert_sort_on_sphere_3 */
} /* end namespace CGAL */
