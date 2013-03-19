 /// \ingroup PkgSurfaceModeling
 /// \cgalConcept
 ///
 /// @brief Concept describing the set of requirements for calculating weights for edges. 
 /// @tparam Polyhedron a model of HalfedgeGraph
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
template<class Polyhedron>
class SurfaceModelingWeightCalculator
{
public:
  /// Function computing edge weight for edge e
  double operator()(edge_descriptor  e, Polyhedron& polyhedron);
};


