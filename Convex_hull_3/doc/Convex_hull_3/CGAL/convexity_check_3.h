namespace CGAL {

/*!
\ingroup PkgConvexityChecking

determines if the vertices of a given polyhedron
represents a strongly convex set of points or not. A set of points is said
to be strongly convex if it consists of only extreme points (i.e.,
vertices of the convex hull).


\tparam PolygonMesh must be a model of the concept`FaceListGraph`.
\tparam Traits must be a model of the concept `IsStronglyConvexTraits_3`.


\cgalHeading{Implementation}

This function implements the tests described in \cgalCite{mnssssu-cgpvg-96} to
determine convexity and requires \f$ O(e + f)\f$ time for a polyhedron with
\f$ e\f$ edges and \f$ f\f$ faces.


*/

template<class PolygonMesh, class Traits>
bool is_strongly_convex_3(PolygonMesh& pm,
const Traits& traits = Default_traits);

} /* namespace CGAL */
