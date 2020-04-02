/*!
\ingroup PkgSurfaceMeshTopologyConcepts
\cgalConcept

The concept `WeightFunctor` defines a functor to calculate the weight of an edge
*/

class WeightFunctor {
public:
/// \name Public types
/// @{

  /// Number type of the weights.
  using Weight_t = unspecified_type;

  /// Dart_const_handle of the input mesh.
  using Dart_const_handle = unspecified_type;
/// @}

/// \name Public member functions
/// @{

  /// Returns the weight of the edge containing `dh`.
  Weight_t operator()(Dart_const_handle dh) const;
/// @}
};
