
/*!
\ingroup PkgPolylineSimplification2Concepts
\cgalConcept

The concept gives the required predicates and constructions 
for the polyline simplification.

\cgalRefines `ConstrainedDelaunayTriangulationTraits_2` 

\cgalHasModel All \cgal Kernels

*/

class PolylineSimplificationTraits_2 {
public:

/*!
A function object to compute a squared distance for two points. 
Provides the operator: 

`FT operator()(Point_2 p, Point_2 q)` 
which returns the squared distance between `p` and `q`. 
*/ 
typedef unspecified_type Compute_squared_distance_2;
};



