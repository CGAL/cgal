/*!
\ingroup PkgMesh_3SecondaryConcepts
\cgalConcept

The concept `IntersectionGeometricTraits_3` provides types and functors 
required to implement a model of `MeshDomain_3`, 
when the domain is described by a simplicial surface mesh 
forming its boundary. 
The concept `IntersectionGeometricTraits_3` mainly provides the detection 
and construction of intersections between segments and triangles. 

\cgalHasModel Any `CGAL::Kernel`.

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
typedef unspecified_type Point_3; 

/*!
Segment type. 
*/ 
typedef unspecified_type Segment_3; 

/*!
Triangle type. 
*/ 
typedef unspecified_type Triangle_3; 

/*!
Function object that detects an intersection between a 3D segment and a 3D triangle. 
Partial model of `::Kernel::DoIntersect_3`. Provides the operators: 
- `bool operator()(Segment_3 seg, Triangle_3 tr)`

- `bool operator()(Triangle_3 tr, Segment_3 seg)`

which returns `true`, iff the triangle and the segment
have a non empty intersection. 
*/ 
typedef unspecified_type Do_intersect_3; 

/*!
Function object that constructs the intersection 
between a 3D segment and a 3D triangle. 
Partial model of `::Kernel::Intersect_3`. Provides the operators: 

- `boost::optional< boost::variant< Point_3, Segment_3 > > operator()(Segment_3 seg, Triangle_3 tr)`

- `boost::optional< boost::variant< Point_3, Segment_3 > > operator()(Triangle_3 tr, Segment_3 seg)`

which computes the intersection between the triangle and the segment.
*/ 
typedef unspecified_type Intersect_3; 

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
