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
 ///   typedef typename boost::graph_traits<Polyhedron>::edge_descriptor edge_descriptor;
 ///
 ///   Uniform_weight(Polyhedron& /*polyhedron*/) { } 
 ///
 ///   double operator()(edge_descriptor e)
 ///   { return 1.0; }
 /// };
 /// \endcode
template<class Polyhedron>
class SurfaceModelingWeightCalculator
{
public:
  /// The edge type 
  typedef Hidden_type edge_descriptor;
  /// Constructor accepting polyhedron as parameter
  SurfaceModelingWeightCalculator(Polyhedron& polyhedron);
  /// Function for computing edge weight for edge e
  double operator()(edge_descriptor e);
};


