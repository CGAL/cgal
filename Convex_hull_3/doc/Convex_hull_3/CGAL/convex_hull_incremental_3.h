namespace CGAL {

/*!
\ingroup PkgConvexHull3Functions

The function `convex_hull_incremental_3` computes the convex hull polyhedron from a set 
of given three-dimensional points. 

This function is provided for completeness and educational 
purposes. When an efficient incremental implementation is needed, 
using `CGAL::Delaunay_triangulation_3` together with 
`CGAL::convex_hull_3_to_polyhedron_3` is highly recommended. 

Requirements 
-------------- 

<OL> 
<LI>`Polyhedron` must provide a type `Polyhedron::Traits` 
that defines the following types 
<UL> 
<LI>`Polyhedron::Traits::R`, which is a model of 
the representation class `R` required by 
`CGAL::Convex_hull_d_traits_3<R>` 
<LI>`Polyhedron::Traits::Point` 
</UL> 
<LI>`InputIterator::value_type` must be the same as 
`Polyhedron::Traits::Point` 
</OL> 

\sa `CGAL::convex_hull_3` 
\sa `CGAL::convex_hull_2` 

Implementation 
-------------- 

This function uses the \f$ d\f$-dimensional convex hull incremental construction 
algorithm \cite cms-frric-93 
with \f$ d\f$ fixed to 3. The algorithm requires \f$ O(n^2)\f$ time in the 
worst case and \f$ O(n \log n)\f$ expected time. 

\sa `CGAL::Convex_hull_d<R>` 

Example 
-------------- 

The following example computes the convex hull of a set of 250 random 
points chosen uniformly in a sphere of radius 100. 

\cgalexample{Convex_hull_3/incremental_hull_3.cpp} 

computes the convex hull polyhedron 
of the points in the range [`first`,`beyond`)
and assigns it to `P`. If `test_correctness` is set to
`true`, the tests described in \cite mnssssu-cgpvg-96 are
used to determine the correctness of the resulting polyhedron.

*/
template <class InputIterator, class Polyhedron>
void
convex_hull_incremental_3(InputIterator first, 
InputIterator beyond, 
Polyhedron& P, 
bool test_correctness = false);

} /* namespace CGAL */
