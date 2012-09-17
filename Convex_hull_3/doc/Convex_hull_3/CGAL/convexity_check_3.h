namespace CGAL {

/*!
\ingroup PkgConvexityChecking

The function `is_strongly_convex_3` determines if the vertices of a given polyhedron 
represents a strongly convex set of points or not. A set of points is said 
to be strongly convex if it consists of only extreme points (<I>i.e.</I>, 
vertices of the convex hull). 

The default traits class is the kernel in which the type 
`Polyhedron_3::Point_3` is defined. 

### Requirements ###

<OL> 
<LI>`Polyhedron_3::Point_3` is equivalent to `Traits::Point_3`. 
<LI>`Traits` is a model of the concept `IsStronlyConvexTraits_3` 
<LI>`Polyhedron_3` must define the following types: 
<UL> 
<LI>`Polyhedron_3::Facet_iterator` 
<LI>`Polyhedron_3::Vertex_iterator` 
</UL> 
and the following member functions: 
<UL> 
<LI>`facets_begin()` 
<LI>`facets_end()` 
<LI>`vertices_begin()` 
<LI>`vertices_end()` 
</UL> 
The vertex type of `Polyhedron_3` must be a model of 
`ConvexHullPolyhedronVertex_3`; 
the facet type must be `ConvexHullPolyhedronFacet_3`. 
</OL> 

\sa `CGAL::is_ccw_strongly_convex_2` 
\sa `CGAL::is_cw_strongly_convex_2` 

### Implementation ###

This function implements the tests described in \cite mnssssu-cgpvg-96 to 
determine convexity and requires \f$ O(e + f)\f$ time for a polyhedron with 
\f$ e\f$ edges and \f$ f\f$ faces. 

determines if the set of vertices of the polyhedron `P` represent
a strongly convex set of points or not.
\pre The equations of the facet planes of the polyhedron must have
already been computed.

*/

template<class Polyhedron_3, class Traits>
bool is_strongly_convex_3(Polyhedron_3& P, 
const Traits& traits = Default_traits);

} /* namespace CGAL */
