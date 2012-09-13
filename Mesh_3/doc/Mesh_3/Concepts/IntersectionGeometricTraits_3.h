/*!
\ingroup PkgMesh_3SecondaryConcepts
\cgalconcept

The concept `IntersectionGeometricTraits_3` provides types and functors 
required to implement a model of `MeshDomain_3`, 
when the domain is described by a simplicial surface mesh 
forming its boundary. 
The concept `IntersectionGeometricTraits_3` mainly provides the detection 
and construction of intersections between segments and triangles. 

\hasModel Any `CGAL::Kernel`.

\sa `BisectionGeometricTraits_3` 
\sa `CGAL::Polyhedral_mesh_domain_3<Polyhedron,IGT,TriangleAccessor>` 

*/
class IntersectionGeometricTraits_3 {
public:

/// \name Types 
/// @{

/*! 
Point type. 
*/ 
typedef Hidden_type Point_3; 

/*! 
Segment type. 
*/ 
typedef Hidden_type Segment_3; 

/*! 
Triangle type. 
*/ 
typedef Hidden_type Triangle_3; 

/*! 
Function object that detects an intersection between a 3D segment and a 3D triangle. 
Provides the operators: 

`bool operator()(Segment_3 seg, Triangle_3 tr)` 

`bool operator()(Triangle_3 tr, Segment_3 seg)` 

which return `true`, iff the triangle and the segment 
have a non empty intersection. 
*/ 
typedef Hidden_type Do_intersect_3; 

/*! 
Function object that constructs the intersection 
between a 3D segment and a 3D triangle. 
Provides the operators: 

`CGAL::Object operator()(Segment_3 seg, Triangle_3 tr)` 

`CGAL::Object operator()(Triangle_3 tr, Segment_3 seg)` 

which computes as a `CGAL::Object` 
the intersection between the triangle and the segment. 
`CGAL::Object` is either a point, a segment or 
an empty object. 
*/ 
typedef Hidden_type Intersect_3; 

/// @} 

/// \name Operations 
/// @{

/*! 
Returns the intersection detection functor. 
*/ 
Do_intersect_3 
do_intersect_3_object(); 

/*! 
Returns the intersection constructor. 
*/ 
Intersect_3 
intersect_3_object(); 

/// @}

}; /* end IntersectionGeometricTraits_3 */
