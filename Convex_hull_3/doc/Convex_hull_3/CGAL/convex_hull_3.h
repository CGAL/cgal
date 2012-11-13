namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

\brief computes the convex hull of the set of points in the range
[`first`, `last`). The polyhedron `P` is cleared, then
the convex hull is stored in `P`.
 
\attention This function does not compute the plane equations of the faces of `P`.

\pre There are at least four points in the range 
[`first`, `last`) not all of which are collinear.

The function `::convex_hull_3()` computes the convex hull of a given set of 
three-dimensional points 
Two versions of this function 
are available. The first can be used when it is known that the result 
will be a polyhedron and the second when a degenerate hull 
may also be possible. 

\cgalRequires `InputIterator::value_type` is equivalent to `Traits::Point_3`. 
\cgalRequires `Traits` is a model of the concept `ConvexHullTraits_3`. 
For the purposes of checking the postcondition that the convex hull 
is valid, `Traits` should also be a model of the concept 
`IsStronglyConvexTraits_3`. 
\cgalRequires `Polyhedron_3` is a model of `ConvexHullPolyhedron_3`.


If the kernel `R` of the points determined by `InputIterator::value_type` 
is a kernel with exact predicates but inexact constructions 
(in practice we check `R::Has_filtered_predicates_tag` is `Tag_true` and `R::FT` is a floating point type), 
then the default traits class of `::convex_hull_3()` is `Convex_hull_traits_3<R>`, and `R` otherwise. 

\sa `CGAL::convex_hull_incremental_3()` 

### Implementation ###

The algorithm implemented by these functions is the quickhull algorithm of 
Barnard <I>et al.</I> \cite bdh-qach-96. 

### Example ###

The following program computes the convex hull of a set of 250 random 
points chosen from a sphere of radius 100. It then determines if the resulting 
hull is a segment or a polyhedron. Notice that the traits class is not 
necessary in the call to `::convex_hull_3()` but is used in the definition 
of `Polyhedron_3`. 

\cgalExample{Convex_hull_3/quickhull_3.cpp} 

*/

template <class InputIterator, class Polyhedron_3, class Traits>
void convex_hull_3(InputIterator first, InputIterator last, Polyhedron_3& P, const Traits& ch_traits = Default_traits);

/*!
\ingroup PkgConvexHull3Functions

\brief computes the convex hull of the set of points in the range
[`first`, `last`). The result, which may be a point, a segment,
a triangle, or a polyhedron, is stored in `ch_object`. 

\attention This function does not compute the plane equations of the faces of `P`
in case the result is a polyhedron.

\cgalRequires `InputIterator::value_type` is equivalent to `Traits::Point_3`. 
\cgalRequires `Traits` is a model of the concept `ConvexHullTraits_3`. 
For the purposes of checking the postcondition that the convex hull 
is valid, `Traits` should also be a model of the concept 
`IsStronglyConvexTraits_3`. 
\cgalRequires `Traits` defines a type `Polyhedron_3` that is a model of 
`ConvexHullPolyhedron_3`. 

If the kernel `R` of the points determined by `InputIterator::value_type` 
is a kernel with exact predicates but inexact constructions 
(in practice we check `R::Has_filtered_predicates_tag` is `Tag_true` and `R::FT` is a floating point type), 
then the default traits class of `::convex_hull_3()` is `Convex_hull_traits_3<R>`, and `R` otherwise. 

*/

template <class InputIterator, class Traits>
void convex_hull_3(InputIterator first, InputIterator last,
Object& ch_object,
const Traits& ch_traits = Default_traits);

} /* namespace CGAL */
