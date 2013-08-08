 /// \ingroup PkgSurfaceModelingConcepts
 /// \cgalConcept
 ///
 /// @brief Concept describing the set of requirements for calculating weights for edges. 
 ///
 /// \code
 /// // a simple model to SurfaceModelingWeightCalculator concept, which provides uniform weights
 /// class Uniform_weight
 /// {
 /// public:
 ///   template<class HalfedgeGraph, class VertexPointMap>
 ///   double operator()(typename boost::graph_traits<HalfedgeGraph>::edge_descriptor  /*e*/, const HalfedgeGraph& /*p*/, VertexPointMap /*v*/)
 ///   { return 1.0; }
 /// };
 /// \endcode
class SurfaceModelingWeightCalculator
{
public:
/// \name Types 
/// @{
  /// a model of HalfedgeGraph
  typedef Hidden_type Halfedge_graph;
/// @} 

/// @{
  /// a model of `ReadWritePropertyMap`</a>  with boost::graph_traits<HalfedgeGraph>::vertex_descriptor as key and `HalfedgeGraph::Point_3` as value type
  typedef Hidden_type VertexPointMap;
/// @} 

/// \name Creation 
/// @{
  /// Default constructor. Required only if the default parameter is used in the constructor of `CGAL::Deform_mesh`.
  SurfaceModelingWeightCalculator();
/// @} 

/// \name Operations 
/// @{
  /// Function computing the edge weight of edge `e`
  double operator()(boost::graph_traits<HalfedgeGraph>::edge_descriptor  e, const HalfedgeGraph& halfedge_graph, VertexPointMap vpm);
/// @}
};


