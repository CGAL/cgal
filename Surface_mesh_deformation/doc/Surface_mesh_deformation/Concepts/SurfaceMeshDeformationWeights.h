 /// \ingroup PkgSurfaceMeshDeformationConcepts
 /// \cgalConcept
 ///
 /// @brief Concept describing the set of requirements for calculating weights for halfedges.
 ///
 /// \cgalHeading{Example:}
 ///
 /// \code
 /// // a simple model to `SurfaceMeshDeformationWeights` concept, which provides uniform weights
 /// template <class TriangleMesh>
 /// struct Identity_weight
 /// {
 ///   typedef TriangleMesh Triangle_mesh;
 ///   template<class VertexPointMap>
 ///   double operator()(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  /*e*/, const TriangleMesh& /*p*/, VertexPointMap /*v*/)
 ///   { return 1.0; }
 /// };
 /// \endcode
 ///
 ///
class SurfaceMeshDeformationWeights
{
public:
/// \name Types
/// @{
  /// a model of `HalfedgeGraph`
  typedef unspecified_type Triangle_mesh;
/// @}

/// \name Creation
/// @{
  /// Default constructor. Required only if the default parameter is used in the constructor of `CGAL::Surface_mesh_deformation`.
  SurfaceMeshDeformationWeights();
/// @}

/// \name Operations
/// @{
  /// Function computing the halfedge weight of halfedge `he`
  /// \tparam VertexPointMap a model of `ReadWritePropertyMap`</a>  with boost::graph_traits<Triangle_mesh>::vertex_descriptor as key and a 3D point from a \cgal Kernel as value type
  template <class VertexPointMap>
  double operator()(boost::graph_traits<Triangle_mesh>::halfedge_descriptor  he, const Triangle_mesh& triangle_mesh, VertexPointMap vpm);
/// @}
};


