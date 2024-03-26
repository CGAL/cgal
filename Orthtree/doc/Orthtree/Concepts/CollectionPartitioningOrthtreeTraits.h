/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  Refinement of the `OrthtreeTraitsWithData` concept, adding requirements for the
  traits class of a `CGAL::Orthtree` in order to supports nearest-neighbor searching.

  Nearest neighbor searches expect a tree where `Node_data` is a model of `ForwardRange`.
  The leaf nodes of the tree represent an exclusive partition of the elements contained in the tree.
  This means that no element should be contained by more than one node.

  \cgalRefines{OrthtreeTraitsWithData}

  \cgalHasModelsBegin
  \cgalHasModels{CGAL::Orthtree_traits_point<GeomTraits, PointRange, PointMap, DimensionTag>}
  \cgalHasModelsEnd
*/
class CollectionPartitioningOrthtreeTraits {
public:


  /// \name Types
  /// @{

  /*!
   * Sphere Type used for the shrinking-sphere approach for neighbor queries; needs to be copy assignable.
   */
  using Sphere_d = unspecified_type;

  /*!
   * \brief The data type contained by each node; must be a model of `ForwardRange` in addition to default constructible, copy constructible and copy assignable.
   */
  using Node_data = unspecified_type;

  /*!
   * \brief An element of the `Node_data` list-like type.
   *
   * Must be constructible from the value type of a `Node_data::iterator`.
   * Typically the same as that type, but can also be an `std::reference_wrapper<>` if the type is not copyable.
   */
  using Node_data_element = unspecified_type;

  /*!
   * \brief Functor with an operator that calculates the squared distance of a `Node_data_element` from a point.
   *
   * Provides the operator:
   * `FT operator()(const Node_data_element&, const Point_d&)`
   */
  using Squared_distance_of_element = unspecified_type;

  /*!
   * \brief Functor with an operator that constructs a `Sphere_d` from a provided center and squared radius.
   *
   * Provides the operator:
   * `Sphere_d operator()(const Point_d&, const FT&)`
  */
  using Construct_sphere_d = unspecified_type;

  /*!
   * \brief Functor with an operator that provides the center of a `Sphere_d`.
   *
   * Provides the operator:
   * `Point_d operator()(const Sphere_d&)`
  */
  using Construct_center_d = unspecified_type;

  /*!
   * \brief Functor with an operator that provides the squared radius of a `Sphere_d`.
   *
   * Provides the operator:
   * `FT operator()(const Sphere_d&)`
  */
  using Compute_squared_radius_d = unspecified_type;

  /// @}

  /// \name Operations
  /// @{

  /*!
   * constructs an object of type `ConstructSphere_d`.
   */
  Construct_sphere_d construct_sphere_d_object() const;

  /*!
   * constructs an object of type `ConstructCenter_d`.
   */
  Construct_center_d construct_center_d_object() const;

  /*!
   * constructs an object of type `ComputeSquaredRadius_d`.
   */
  Compute_squared_radius_d compute_squared_radius_d_object() const;

  /*!
   * constructs an object of type `Squared_distance_of_element`.
   */
  Squared_distance_of_element squared_distance_of_element_object() const;

  /// @}
};
