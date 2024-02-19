/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  The concept `OrthtreeTraitsWithData` defines the requirements for the
  template parameter of the `CGAL::Orthtree` class for a node type that stores data.

  \cgalRefines{OrthtreeTraits}

  \cgalHasModelsBegin
  \cgalHasModels{CGAL::Orthtree_traits_point<GeomTraits, PointRange, PointMap, dimension>}
  \cgalHasModels{CGAL::Orthtree_traits_face_graph<PolygonMesh, VPM>}
  \cgalHasModels{CGAL::Orthtree_traits_base<K, dimension>}
  \cgalHasModelsEnd
*/
class OrthtreeTraitsWithData
{
public:

  /// \name Types
  /// @{
  /*!
   * \brief The data type contained by each node.
   */
  using Node_data = unspecified_type;

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
   * \brief Functor which fills the contents of the nodes children.
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

  /// @}

  /// \name Operations
  /// @{

  /*!
   * constructs an object of type `Construct_root_node_contents`.
   */
  Construct_root_node_contents construct_root_node_contents_object() const;

  /*!
   * constructs an object of type `Distribute_node_contents`.
   */
  Distribute_node_contents distribute_node_contents_object() const;

  /// @}
};
