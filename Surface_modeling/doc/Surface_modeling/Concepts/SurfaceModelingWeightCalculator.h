 /// \ingroup PkgSurfaceModelingConcepts
 /// \cgalConcept
 ///
 /// @brief Concept describing the set of requirements for calculating weights for edges. 
 ///
 /// \code
 /// // a simple model to SurfaceModelingWeightCalculator concept, which provides uniform weights
 /// template<class Polyhedron>
 /// class Uniform_weight
 /// {
 /// public:
 ///   double operator()(typename boost::graph_traits<Polyhedron>::edge_descriptor  /*e*/, Polyhedron& /*polyhedron*/)
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

/// \name Creation 
/// @{
  /// Default constructor. Required only if the default parameter is used in the constructor of `CGAL::Deform_mesh`.
  SurfaceModelingWeightCalculator();
/// @} 

/// \name Operations 
/// @{
  /// Function computing the edge weight of edge `e`
  double operator()(boost::graph_traits<Polyhedron>::edge_descriptor  e, Polyhedron& polyhedron);
/// @}
};


