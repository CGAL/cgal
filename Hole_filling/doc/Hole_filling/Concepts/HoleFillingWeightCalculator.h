 /// \ingroup PkgHoleFillingConcepts
 /// \cgalConcept
 ///
 /// @brief Concept describing the set of requirements for calculating weights for edges and vertices. 

class FairWeightCalculator
{
public:
/// \name Types 
/// @{
  /// a model of Polyhedron
  typedef Hidden_type Polyhedron;
/// @} 


/// \name Operations 
/// @{
  /// Function computing the edge weight of edge `e`
  double w_ij(Polyhedron::Halfedge_handle  e, const Polyhedron& polyhedron);
  
  /// Function computing the vertex weight of vertex `v`
  double w_i(Polyhedron::Vertex_handle  v, const Polyhedron& polyhedron);
/// @}
};


