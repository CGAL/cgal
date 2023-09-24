/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  In addition to the requirements described in the OrthtreeTraits concept,
  the concept `CollectionPartitioningOrthtreeTraits` defines the requirements for the
  traits class of a `CGAL::Orthtree` which supports nearest-neighbor searching.

  Nearest neighbor searches expect a tree with nodes which contain list types.
  The leaf nodes of the tree represent an exclusive partition of the elements contained in the tree.
  This means that no element should be contained by more than one node.

  \cgalRefines{OrthtreeTraits}

  \cgalHasModelsBegin
  \cgalHasModels{CGAL::Orthtree_traits_point<GeomTraits, PointSet, PointMap, DimensionTag>}
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
   * \brief Functor with an operator to produce a geometric object from a `Node_data_element`.
   *
   * The return type of the functor must be a valid argument to `CGAL::squared_distance`.
   */
  using Get_geometric_object_for_element = unspecified_type;

  /// @}

  /// \name Operations
  /// @{

  /*!
   * Function used to construct an object of type `Get_geometric_object_for_element`.
   */
  Get_geometric_object_for_element get_geometric_object_for_element_object() const;

  /// @}
};
