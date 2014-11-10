namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

\brief computes the convex hull of the set of points in the range
[`first`, `last`). The polyhedron `P` is cleared, then
the convex hull is stored in `P`. Note that the convex hull will be triangulated,
that is `P` will contain only triangular facets.
 
\attention This function does not compute the plane equations of the faces of `P`.

\pre There are at least four points in the range 
[`first`, `last`) not all of which are collinear.

\tparam InputIterator must be an input iterator with a value type  equivalent to `Traits::Point_3`. 
\tparam Polyhedron_3 must be a model of `ConvexHullPolyhedron_3`.
\tparam Traits must be a model of the concept `ConvexHullTraits_3`. 
For the purposes of checking the postcondition that the convex hull 
is valid, `Traits` should also be a model of the concept 
`IsStronglyConvexTraits_3`. 



If the kernel `R` of the points determined by the value type of `InputIterator` 
is a kernel with exact predicates but inexact constructions 
(in practice we check `R::Has_filtered_predicates_tag` is `Tag_true` and `R::FT` is a floating point type), 
then the default traits class of `::convex_hull_3()` is `Convex_hull_traits_3<R>`, and `R` otherwise. 

\sa `convex_hull_incremental_3()` 

\cgalHeading{Implementation}

The algorithm implemented by these functions is the quickhull algorithm of 
Barnard <I>et al.</I> \cgalCite{bdh-qach-96}. 


*/

template <class InputIterator, class Polyhedron_3, class Traits>
void convex_hull_3(InputIterator first, InputIterator last, Polyhedron_3& P, const Traits& ch_traits = Default_traits);

/*!
\ingroup PkgConvexHull3Functions

\brief computes the convex hull of the set of points in the range
[`first`, `last`). The result, which may be a point, a segment,
a triangle, or a polyhedron, is stored in `ch_object`. 
In the case the result is a polyhedron, the convex hull will be triangulated,
that is the polyhedron will contain only triangular facets.

\attention This function does not compute the plane equations of the faces of `P`
in case the result is a polyhedron.

\tparam InputIterator must be an input iterator with a value type  equivalent to `Traits::Point_3`. 
\tparam Traits must be model of the concept `ConvexHullTraits_3`. 
For the purposes of checking the postcondition that the convex hull 
is valid, `Traits` should also be a model of the concept 
`IsStronglyConvexTraits_3`.   Furthermore, `Traits` must define a type `Polyhedron_3` that is a model of 
`ConvexHullPolyhedron_3`. 

If the kernel `R` of the points determined by the value type  of `InputIterator` 
is a kernel with exact predicates but inexact constructions 
(in practice we check `R::Has_filtered_predicates_tag` is `Tag_true` and `R::FT` is a floating point type), 
then the default traits class of `convex_hull_3()` is `Convex_hull_traits_3<R>`, and `R` otherwise. 

*/

template <class InputIterator, class Traits>
void convex_hull_3(InputIterator first, InputIterator last,
Object& ch_object,
const Traits& ch_traits = Default_traits);

} /* namespace CGAL */
