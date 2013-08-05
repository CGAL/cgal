namespace CGAL {

/*!
\ingroup PkgConvexityChecking

determines if the vertices of a given polyhedron 
represents a strongly convex set of points or not. A set of points is said 
to be strongly convex if it consists of only extreme points (i.e., 
vertices of the convex hull). 

The default traits class is the kernel in which the type 
`Polyhedron::Point_3` is defined. 

\cgalHeading{Requirements}

<OL> 
<LI>`Polyhedron::Point_3` is equivalent to `Traits::Point_3`. 
<LI>`Traits` is a model of the concept `IsStronlyConvexTraits_3` 
<LI>`Polyhedron` must define the following types: 
<UL> 
<LI>`Polyhedron::Facet_iterator` 
<LI>`Polyhedron::Vertex_iterator` 
</UL> 
and the following member functions: 
<UL> 
<LI>`facets_begin()` 
<LI>`facets_end()` 
<LI>`vertices_begin()` 
<LI>`vertices_end()` 
</UL> 
The vertex type of `Polyhedron` must be a model of 
`ConvexHullPolyhedronVertex_3`; 
the facet type must be `ConvexHullPolyhedronFacet_3`. 
</OL> 


\cgalHeading{Implementation}

This function implements the tests described in \cgalCite{mnssssu-cgpvg-96} to 
determine convexity and requires \f$ O(e + f)\f$ time for a polyhedron with 
\f$ e\f$ edges and \f$ f\f$ faces. 

determines if the set of vertices of the polyhedron `P` represent
a strongly convex set of points or not.
\pre The equations of the facet planes of the polyhedron must have
already been computed.

*/

template<class Polyhedron, class Traits>
bool is_strongly_convex_3(Polyhedron& P, 
const Traits& traits = Default_traits);

} /* namespace CGAL */
