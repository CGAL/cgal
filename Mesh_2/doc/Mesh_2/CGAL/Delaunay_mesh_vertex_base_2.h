
namespace CGAL {

/*!
\ingroup PkgMesh2Ref

The class `Delaunay_mesh_vertex_base_2` is a model for the concept
`DelaunayMeshVertexBase_2`.

This class can be used directly or it can serve as a base to derive other
classes with some additional attributes tuned to a
specific application.

\tparam Traits is the geometric traits class. It must
be the same as the one used for the Delaunay mesh.

\tparam Vb is the base class from which `Delaunay_mesh_vertex_base_2`
derives. It must be a model of the `TriangulationVertexBase_2` concept.

\cgalModels{DelaunayMeshVertexBase_2}

*/
template< typename Traits, typename Vb >
class Delaunay_mesh_vertex_base_2 : Vb {
public:

}; /* end Delaunay_mesh_vertex_base_2 */
} /* end namespace CGAL */
