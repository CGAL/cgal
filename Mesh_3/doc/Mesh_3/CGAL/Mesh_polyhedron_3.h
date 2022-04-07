namespace CGAL {

/*!
\ingroup PkgMesh3Domains

The class `Mesh_polyhedron_3` provides a customized `Polyhedron_3` type. This type uses
as `PolyhedronItems_3` a customized type which adds data to the Vertex, Face and
Halfedge class. Those data are required to use our sharp features
detection algorithm.

\tparam IGT stands for the geometric traits associated
to the meshing process. It should be a model of the two concepts
`PolyhedronTraits_3` and `IntersectionGeometricTraits_3`.

\sa `CGAL::Polyhedron_3<Gt>`
\sa `CGAL::Polyhedral_mesh_domain_with_features_3<IGT>`

*/
template< typename IGT >
struct Mesh_polyhedron_3 {

/// \name Types
/// @{

/*!
`CGAL::Polyhedron_3<IGT>` type with customized `PolyhedronItems_3`
designed to handle sharp feature detection.
*/
typedef unspecified_type type;

/// @}

}; /* end Mesh_polyhedron_3 */
} /* end namespace CGAL */
