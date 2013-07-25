namespace CGAL {

/*!
\ingroup PkgMesh_3Domains

The class `Polyhedral_mesh_domain_3` implements 
a domain defined by a simplicial polyhedral surface.
	
The input polyhedral surface must be free of intersection.
It must include (at least) one closed connected component
that defines the boundary of the domain to be meshed.
Inside this bounding component,
the input polyhedral surface may contain
several other components (closed or not)
that will be represented in the final mesh.   
This class is a model of the concept `MeshDomain_3`.

\tparam Polyhedron stands for the type of the input polyhedral surface(s). 
The only requirements for this type is that the triangles of the surfaces 
must be accessible through an object of the class 
`TriangleAccessor`. 

\tparam IGT stands for a geometric traits class 
providing the types and functors required to implement 
the intersection tests and intersection computations 
for polyhedral boundary surfaces. This parameter has to be instantiated 
with a model of the concept `IntersectionGeometricTraits_3`. 

\tparam TriangleAccessor provides access to the triangles 
of the input polyhedral 
surface. It must be a model of the concept 
`TriangleAccessor_3`. It defaults to 
`Triangle_accessor_3<Polyhedron,IGT>`. The type `IGT::Triangle_3` must 
be identical to the type `TriangleAccessor::Triangle_3`. 

\cgalModels `MeshDomain_3` 

\sa `TriangleAccessor_3` 
\sa `IntersectionGeometricTraits_3` 
\sa `CGAL::Triangle_accessor_3<Polyhedron_3<K>,K>` 
\sa `CGAL::make_mesh_3()`. 

*/
template< typename Polyhedron, typename IGT, typename TriangleAccessor >
class Polyhedral_mesh_domain_3 {
public:

/// \name Creation 
/// @{

/*!
Construction from a bouding polyhedral surface which should be closed, and free of intersections.
The inside of `bounding_polyhedron` will be meshed.
*/ 
Polyhedral_mesh_domain_3(Polyhedron bounding_polyhedron); 

/*!
Construction from a polyhedral surface, and a bounding polyhedral surface,.
The first polyhedron should be entirely included inside `bounding_polyhedron`, which has to be closed 
and free of intersections. 
Using this constructor allows to mesh a polyhedral surface which is not closed, or has holes.
The inside of `bounding_polyhedron` will be meshed.
*/
Polyhedral_mesh_domain_3(Polyhedron polyhedron,
				     Polyhedron bounding_polyhedron);
			   

/// @}

}; /* end Polyhedral_mesh_domain_3 */
} /* end namespace CGAL */
