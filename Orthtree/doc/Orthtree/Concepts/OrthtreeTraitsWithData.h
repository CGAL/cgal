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
   * \brief The data type contained by each node. Must be default constructible, copy constructible and copy assignable.
   */
  using Node_data = unspecified_type;

  /*!
   * \brief Functor which initializes elements contained by the root node.
   *
   * Each node of a tree has an associated `Node_data` value.
   * This functor initializes the `Node_data` of the root node.
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
   * `void operator()(Node_index, Orthtree<Traits>&, const Point_d&)`
   *
   * The functor is called during refinement of the `Orthtree` on a node after it has been split. The purpose of the functor is to
   * distribute the `Node_data`, accessible via `tree.data()`, to the data of the nodes children, accessible via `tree.children()`.
   * The first parameter is the `Node_index` of the node. The second parameter provides the instance of the `Orthtree`
   * and the last parameter is the barycenter of the node which will be used as shared corner amongst the children of the node.
   *
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
