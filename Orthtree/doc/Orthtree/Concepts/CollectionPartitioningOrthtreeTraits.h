/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  Refinement of the OrthtreeTraits concept, adding requirements for the
  traits class of a `CGAL::Orthtree` in order to supports nearest-neighbor searching.

  Nearest neighbor searches expect a tree where `Node_data` is a model of `ForwardRange`.
  The leaf nodes of the tree represent an exclusive partition of the elements contained in the tree.
  This means that no element should be contained by more than one node.

  \cgalRefines{OrthtreeTraits}

  \cgalHasModelsBegin
  \cgalHasModels{CGAL::Orthtree_traits_point<GeomTraits, PointRange, PointMap, DimensionTag>}
  \cgalHasModelsEnd
*/
class CollectionPartitioningOrthtreeTraits {
public:


  /// \name Types
  /// @{

  /*!
   * Sphere type used for the shrinking-sphere approach for neighbor queries
   */
  using Sphere_d = unspecified_type;

  /*!
   * \brief An element of the `Node_data` list-like type.
   *
   * Must be constructible from the type produced by dereferencing a `Node_data` iterator.
   * Typically the same as that type, but can also be an `std::reference_wrapper<>` if the type is not copyable.
   */
  using Node_data_element = unspecified_type;

  /*!
   * \brief Functor with an operator calculate the squared distance of `Node_data_element` from a point.
   *
   * Provides the operator:
   * `FT operator()(const Node_data_element&, const Point_d&)`
   */
  using Squared_distance_of_element = unspecified_type;

  /// @}

  /// \name Operations
  /// @{

  /*!
   * constructs an object of type `Squared_distance_of_element`.
   */
  Squared_distance_of_element get_squared_distance_of_element_object() const;

  /// @}
};
