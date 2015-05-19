namespace CGAL {

/*!
\ingroup PkgSpatialSortingFunctions

The function `spatial_sort_on_sphere()` sorts an iterator range of points in a way that 
improves space locality with respect to the intrinsic metric on the sphere given as input. 
Two points close in the order will be close 
on the sphere, and two points close on the sphere will have a high probability 
of being close in the order. The input points are supposed to be close to the input sphere.
If input points are not close to the input sphere, this function still works, but it might not be a good sorting function.

It sorts the range `[begin, end)` in place. 

The default traits class `Default_traits` is the kernel in which the type 
`std::iterator_traits<RandomAccessIterator>::%value_type` is defined. 

The default policy is `Hilbert_sort_median_policy()` and the 
other option is `Hilbert_sort_middle_policy()`. 

The input sphere is given by a squared radius and a center, parameter `sqr_radius` and parameter `center` respectively.
The default squared radius of the sphere is 1.0.
The default center of the sphere is the origin (0,0,0).

\cgalHeading{Requirements}

<OL> 
<LI>`std::iterator_traits<RandomAccessIterator>::%value_type` is convertible to 
`Traits::Point_3`. 
<LI>`Traits` is a model for concept `SpatialSortingTraits_3`. 
</OL> 

\cgalHeading{Precondition}

<OL> 
<LI>`sqr_radius` greater than 0. 
</OL>


\cgalHeading{Implementation}

Creates an instance of `Multiscale_sort<Hilbert_sort_on_sphere_3>` 
where `Hilbert_sort_on_sphere_3` is an Hilbert sorting on the sphere object, 
and calls its `operator()`. 

The `threshold_hilbert` is the minimal size of a point set to be 
subdivided recursively during Hilbert sorting, otherwise random order is used. 
The `threshold_multiscale` value is the minimal size for a sample to 
call the `Hilbert_sort_on_sphere_3` functor, otherwise random order is used. 
The `ratio` value is used to split the original set in two subsets, 
spatial sort is applied on the first subset of size 
`ratio` 
times the original size of the set, `Hilbert_sort_on_sphere_3` functor is applied on the 
second subset. 

*/
template <class RandomAccessIterator, class Traits, class PolicyTag>
void
spatial_sort_on_sphere( RandomAccessIterator begin,
RandomAccessIterator end,
const Traits& traits = Default_traits,
PolicyTag policy = Default_policy,
double sqr_radius = 1.0,
const Traits::Point_3 &center = Default_center,
std::ptrdiff_t threshold_hilbert=default,
std::ptrdiff_t threshold_multiscale=default,
double ratio=default);

} /* namespace CGAL */

