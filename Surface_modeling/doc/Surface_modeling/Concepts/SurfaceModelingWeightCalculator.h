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
 ///   template<class Polyhedron, class VertexPointMap>
 ///   double operator()(typename boost::graph_traits<Polyhedron>::edge_descriptor  /*e*/, const Polyhedron& /*p*/, VertexPointMap /*v*/)
 ///   { return 1.0; }
 /// };
 /// \endcode
class SurfaceModelingWeightCalculator
{
public:
/// \name Types 
/// @{
  /// a model of HalfedgeGraph
  typedef Hidden_type Polyhedron;
/// @} 

/// @{
  /// a model of `ReadWritePropertyMap`</a>  with boost::graph_traits<Polyhedron>::vertex_descriptor as key and `Polyhedron::Point_3` as value type
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
  double operator()(boost::graph_traits<Polyhedron>::edge_descriptor  e, const Polyhedron& polyhedron, VertexPointMap vpm);
/// @}
};


