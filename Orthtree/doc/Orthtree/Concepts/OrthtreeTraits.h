/*!
  \ingroup PkgOrthtreeConcepts
  \cgalConcept

  The concept `OrthtreeTraits` defines the requirements for the
  template parameter of the `CGAL::Orthtree` class.

  todo: update list of models
  \cgalHasModel `CGAL::Orthtree_traits_2<GeomTraits>`
  \cgalHasModel `CGAL::Orthtree_traits_3<GeomTraits>`
  \cgalHasModel `CGAL::Orthtree_traits_d<GeomTraits,Dimension>`
*/
class OrthtreeTraits
{
public:

  /// \name Types
  /// @{

  typedef unspecified_type Dimension; ///< Dimension type (see `CGAL::Dimension_tag`).
  typedef unspecified_type FT; ///< The number type of the %Cartesian coordinates of types `Point_d`
  typedef unspecified_type Point_d; ///< Point type.
  typedef unspecified_type Bbox_d; ///< Bounding box type. Must be constructible from a pair of Point_d types.

  /*!
    A random access iterator type to enumerate the
    %Cartesian coordinates of a point.
  */
  typedef unspecified_type Cartesian_const_iterator_d;
  typedef std::array<FT, Dimension::value> Array; ///< Array used for easy point constructions.


  /*!
   * \brief The data type contained by each node.
   */
  typedef unspecified_type Node_data;

  typedef unspecified_type Adjacency; ///< Specify the adjacency directions

  /*!
    Functor with an operator to construct a `Point_d` from an `Array` object.
  */
  typedef unspecified_type Construct_point_d_from_array;

  /// @}

  /// \name Operations
  /// @{

  /*!
    Function used to construct an object of type `Construct_point_d_from_array`.
  */
  Construct_point_d_from_array construct_point_d_from_array_object() const;

  /*!
   * \brief Produces a bounding box which encloses the contents of the tree
   *
   * The bounding box must enclose all elements contained by the tree.
   * It may be tight-fitting, the orthtree constructor produces a bounding cube surrounding this region.
   * For traits which assign no data to each node, this can be defined to return a fixed region.
   *
   * @return std::pair<min, max>, where min and max represent cartesian corners which define a bounding box
   */
  std::pair<Array, Array> root_node_bbox() const;

  /*!
   * \brief Initializes the contained elements for the root node.
   *
   * Typically produces a `Node_data` which contains all the elements in the tree.
   * e.g. For a tree where each node contains a set of points,
   * root_node_contents() will produce the list of all points.
   *
   * @return The `Node_data` instance to be contained by the root node
   */
  Node_data root_node_contents() const;

  /*!
   * \brief Distributes the `Node_data` contents of a node to its immediate children.
   *
   * Invoked after a node is split.
   * Adds the contents of the node n to each of its children.
   * May rearrange or modify n's `Node_data`, but generally expected not to reset n.
   * After distributing n's contents, n should still have an list of elements it encloses.
   * Each of n's children should have an accurate list of the subset of elements within n they enclose.
   *
   * For an empty tree, this can be a null-op.
   *
   * @tparam Node_index The index type used by an orthtree implementation
   * @tparam Tree An Orthree implementation
   *
   * @param n The index of the node who's contents must be distributed.
   * @param tree The Orthtree which n belongs to
   * @param center The coordinate center of the node n, which its contents should be split around.
   */
  template <typename Node_index, typename Tree>
  void distribute_node_contents(Node_index n, Tree& tree, const Point_d& center);

  /// @}
};
