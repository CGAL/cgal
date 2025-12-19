/*!
\ingroup PkgSurfaceMeshTopologyConcepts
\cgalConcept

The concept `WeightFunctor` defines a functor to calculate the weight of an edge.


  \cgalHasModelsBegin
  \cgalHasModelsBare{\link CGAL::Surface_mesh_topology::Unit_weight_functor `CGAL::Surface_mesh_topology::Unit_weight_functor`\endlink}
  \cgalHasModelsBare{\link CGAL::Surface_mesh_topology::Euclidean_length_weight_functor `CGAL::Surface_mesh_topology::Euclidean_length_weight_functor<Mesh>`\endlink}
  \cgalHasModelsEnd
*/
class WeightFunctor {
public:
/// \name Public types
/// @{

  /*!
    A descriptor to `Dart` for combinatorial/generalized maps, or a halfedge descriptor for models of the `FaceGraph` concept.
  */
  typedef unspecified_type halfedge_descriptor;

  /// Number type of the weights.
  using Weight_t = unspecified_type;
/// @}

/// \name Public member functions
/// @{

  /// Returns the weight of the edge containing `hd`.
  Weight_t operator()(halfedge_descriptor hd) const;
/// @}
};
