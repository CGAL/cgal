/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  The concept `OrthtreeTraits` defines the requirements for the
  template parameter of the `CGAL::Orthtree` class.

  \cgalHasModelsBegin
  \cgalHasModels{CGAL::Orthtree_traits_point<GeomTraits, PointRange, PointMap, dimension >}
  \cgalHasModels{CGAL::Orthtree_traits_face_graph<PolygonMesh, VPM>}
  \cgalHasModels{CGAL::Orthtree_traits_base< K, dimension >}
  \cgalHasModelsEnd
*/
class OrthtreeTraits
{
public:

  /// \name Types
  /// @{
  using Node_index = unspecified_type; ///< An integer type for nodes
  constexpr int dimension; ///< Dimension.
  using FT = unspecified_type; ///< The number type of the %Cartesian coordinates of types `Point_d`
  using Point_d = unspecified_type; ///< Point type.
  using Bbox_d = unspecified_type; ///< Bounding box type. Must be constructible from a pair of `Point_d` types.

  /*!
    A random access iterator type to enumerate the
    %Cartesian coordinates of a point of type `Point_d`.
  */
  using Cartesian_const_iterator_d = unspecified_type;

  /*!
   * \brief The data type contained by each node.
   */
  using Node_data = unspecified_type;

  /*!
   * \brief Integral number type which can take on values indicating adjacency directions.
   *
   * Must be able to take on values ranging from 0 to the number of faces of the (hyper)rectangle type, equivalent to 2 * D.
   */
  using Adjacency = unspecified_type; ///< Specify the adjacency directions

  /*!
   * \brief Functor with an operator to create the bounding box of the root node.
   *
   * Provides the operator:
   * `Bbox_d operator()()`
   *
   * The bounding box must enclose all elements contained by the tree.
   * It may be tight-fitting. The orthtree constructor produces a bounding box surrounding this region.
   * For traits which assign no data to each node, this can be defined to return a fixed region.
   */
  using Construct_root_node_bbox = unspecified_type;

  /*!
   * \brief Functor which initializes the contained elements for the root node.
   *
   * Each node of a tree has an associated `Node_data` value.
   * For most nodes, this is set by `Distribute_node_contents`, but that is not possible for the root node.
   * Instead, this functor initializes the `Node_data` of the root node.
   * It takes no arguments, and returns an instance of `Node_data`.
   *
   * Provides the operator:
   * `Node_data operator()()`
   *
   * Typically, the `Node_data` of the root node contains all the elements in the tree.
   * For a tree in which each node contains a span (such as `std::span()`) this function would return the span containing all items.
   *
   */
  using Construct_root_node_contents = unspecified_type;

  /*!
   * \brief functor which fills the contents of the nodes children.
   *
   * Provides the operator:
   * `void operator()(typename Tree::Node_index, Tree&, const Point_d&)`
   *
   * It can use `tree.children(node_index)` to access the children of the node, and `tree.data(node_index)`
   * to access its children and the contents of the node.
   * It must distribute the contents of the node to each of its children.
   * For a tree in which each node contains a span, this may mean rearranging the contents of the original node
   * and producing spans containing a subset of its contents for each of its children.
   * For compatibility with locate, the center of the node is considered to be part of the upper half.
   */
  using Distribute_node_contents = unspecified_type;

  /*!
   * \brief Functor with an operator to construct a `Point_d` from an initializer list.
   *
   * For trees which use a different kernel for the bounding box type,
   * the return type of this functor must match the kernel used by the bounding box type and not that of the contents.
   */
  using Construct_point_d = unspecified_type;

  /// @}

  /// \name Operations
  /// @{

  /*!
   * constructs an object of type `Construct_root_node_bbox`.
   */
  Construct_root_node_bbox construct_root_node_bbox_object() const;

  /*!
   * constructs an object of type `Construct_root_node_contents`.
   */
  Construct_root_node_contents construct_root_node_contents_object() const;

  /*!
   * constructs an object of type `Distribute_node_contents`.
   */
  Distribute_node_contents distribute_node_contents_object() const;

  /*!
   * constructs an object of type `Construct_point_d`.
   */
  Construct_point_d construct_point_d_object() const;

  /// @}
};
