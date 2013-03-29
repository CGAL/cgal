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
	typedef Hidden_type Polyhedron;
	typedef Hidden_type edge_descriptor;
/// @} 

/// \name Creation 
/// @{
  /// Default constructor. Required only if default parameter is used in CGAL::Deform_mesh().
  SurfaceModelingWeightCalculator();
/// @} 

/// \name Operations 
/// @{
  /// Function computing edge weight for edge e
  double operator()(edge_descriptor  e, Polyhedron& polyhedron);
/// @}
};


