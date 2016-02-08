 /// \ingroup PkgSurfaceMeshDeformationConcepts
 /// \cgalConcept
 ///
 /// @brief Concept describing the set of requirements for calculating weights for halfedges.
 ///
 /// \cgalHeading{Example:}
 ///
 /// \code
 /// // a simple model to SurfaceMeshDeformationWeights concept, which provides uniform weights
 /// template <class HalfedgeGraph>
 /// struct Identity_weight
 /// {
 ///   typedef HalfedgeGraph Halfedge_graph;
 ///   template<class VertexPointMap>
 ///   double operator()(typename boost::graph_traits<HalfedgeGraph>::halfedge_descriptor  /*e*/, const HalfedgeGraph& /*p*/, VertexPointMap /*v*/)
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
  /// a model of HalfedgeGraph
  typedef unspecified_type Halfedge_graph;
/// @}

/// \name Creation
/// @{
  /// Default constructor. Required only if the default parameter is used in the constructor of `CGAL::Surface_mesh_deformation`.
  SurfaceMeshDeformationWeights();
/// @}

/// \name Operations
/// @{
  /// Function computing the halfedge weight of halfedge `he`
  /// \tparam VertexPointMap a model of `ReadWritePropertyMap`</a>  with boost::graph_traits<Halfedge_graph>::vertex_descriptor as key and a 3D point from a \cgal Kernel as value type
  template <class VertexPointMap>
  double operator()(boost::graph_traits<Halfedge_graph>::halfedge_descriptor  he, const Halfedge_graph& halfedge_graph, VertexPointMap vpm);
/// @}
};


