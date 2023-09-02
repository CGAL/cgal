// Copyright (c) 2007-2020  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro, Simon Giraudot, Cédric Portaneri, Tong Zhao

#ifndef CGAL_ORTHTREE_H
#define CGAL_ORTHTREE_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Orthtree/Cartesian_ranges.h>
#include <CGAL/Orthtree/Split_predicates.h>
#include <CGAL/Orthtree/Traversals.h>
#include <CGAL/Orthtree/Traversal_iterator.h>
#include <CGAL/Orthtree/IO.h>

#include <CGAL/Property_container.h>
#include <CGAL/property_map.h>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Dimension.h>

#include <boost/function.hpp>
#include <boost/core/span.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>

#include <iostream>
#include <fstream>
#include <ostream>
#include <functional>

#include <bitset>
#include <stack>
#include <queue>
#include <vector>
#include <math.h>

namespace CGAL {

/*!
  \ingroup PkgOrthtreeClasses

  \brief A data structure using an axis-aligned hybercubic
  decomposition of dD space for efficient point access and
  computations.

  \details It builds a hierarchy of nodes which subdivide the space
  based on a collection of points.  Each node represents an
  axis-aligned hypercubic region of space.  A node contains the range
  of points that are present in the region it defines, and it may
  contain \f$2^{dim}\f$ other nodes which further subdivide the
  region.

  \sa `CGAL::Quadtree`
  \sa `CGAL::Octree`

  \tparam Traits_ must be a model of `OrthtreeTraits`
 */
template <typename Traits_>
class Orthtree {

public:

  /// \name Template Types
  /// @{
  typedef Traits_ Traits; ///< Geometry traits
  /// @}

  /// \name Traits Types
  /// @{
  typedef typename Traits::Dimension Dimension; ///< Dimension of the tree
  typedef typename Traits::FT FT; ///< Number type.
  typedef typename Traits::Point_d Point; ///< Point type.
  typedef typename Traits::Bbox_d Bbox; ///< Bounding box type.
  typedef typename Traits::Sphere_d Sphere; ///< Sphere type.
  typedef typename Traits::Adjacency Adjacency; ///< Adjacency type.

  typedef typename Traits::Node_data Node_data;
  // todo: Node_data_element will only exist for certain Traits types, so I don't know if it can be re-exported

  /// \cond SKIP_IN_MANUAL
  typedef typename Traits::Array Array; ///< Array type.
  /// \endcond
  /// @}

  /// \name Public Types
  /// @{

  /*!
   * \brief Self typedef for convenience.
   */
  typedef Orthtree<Traits> Self;

  /*!
   * \brief Degree of the tree (number of children of non-leaf nodes).
   */
  typedef Dimension_tag<(2 << (Dimension::value - 1))> Degree;

  /*!
   * \brief Index of a given node in the tree; the root always has index 0.
   */
  typedef std::size_t Node_index;

  /*!
   * \brief Optional index of a node in the tree.
   */
  typedef std::optional<Node_index> Maybe_node_index;

  // todo: maybe this could be private?
  typedef Properties::Property_container<Node_index> Node_property_container;

  /*!
    \brief Set of bits representing this node's relationship to its parent.

    Equivalent to an array of Booleans, where index[0] is whether `x`
    is greater, index[1] is whether `y` is greater, index[2] is whether
    `z` is greater, and so on for higher dimensions if needed.
    Used to represent a node's relationship to the center of its parent.
   */
  typedef std::bitset<Dimension::value> Local_coordinates;

  /*!
    \brief Coordinates representing this node's relationship
    with the rest of the tree.

    Each value `(x, y, z, ...)` of global coordinates is calculated by doubling
    the parent's global coordinates and adding the local coordinates.
   */
  typedef std::array<std::uint32_t, Dimension::value> Global_coordinates;

  /*!
   * \brief A predicate that determines whether a node must be split when refining a tree.
   */
  typedef std::function<bool(Node_index, const Self&)> Split_predicate;

  /*!
   * \brief A model of `ConstRange` whose value type is `Node_index`.
   */
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Node_range;
  typedef unspecified_type Node_index_range;
#else
  typedef boost::iterator_range<Index_traversal_iterator<Self>> Node_index_range;
#endif

  /// \cond SKIP_IN_MANUAL
  typedef Orthtrees::internal::Cartesian_ranges<Traits> Cartesian_ranges;
  /// \endcond

  /// @}

private: // data members :

  Traits m_traits; /* the tree traits */

  Node_property_container m_node_properties;
  Node_property_container::Array <Node_data>& m_node_points;
  Node_property_container::Array <std::uint8_t>& m_node_depths;
  Node_property_container::Array <Global_coordinates>& m_node_coordinates;
  Node_property_container::Array <Maybe_node_index>& m_node_parents;
  Node_property_container::Array <Maybe_node_index>& m_node_children;

  Point m_bbox_min;                  /* input bounding box min value */

  std::vector<FT> m_side_per_depth;      /* side length per node's depth */

  Cartesian_ranges cartesian_range; /* a helper to easily iterator on coordinates of points */

public:

  /// \name Constructor
  /// @{

  /*!
    \brief creates an orthtree from a collection of points.

    The constructed orthtree has a root node, with no children, that
    contains the points passed. That root node has a bounding box that
    encloses all of the points passed, with padding according to the
    `enlarge_ratio`.

    That root node is built as a `n`-cube (square in 2D, cube in 3D,
    etc.) whose edge size is the longest bounding box edge multiplied
    by `enlarge_ratio`. Using an `enlarge_ratio>1.0` prevents some
    points from being exactly on the border of some cells, which can
    lead to over-refinement.

    This single-node orthtree is valid and compatible
    with all Orthtree functionality, but any performance benefits are
    unlikely to be realized until `refine()` is called.

    \warning The input point range is not copied. It is used directly
    and is rearranged by the `Orthtree`. Altering the point range
    after creating the orthtree might leave it in an invalid state.

    \param point_range input point range.
    \param point_map property map to access the input points.
    \param enlarge_ratio ratio to which the bounding box should be enlarged.
    \param traits the traits object.
  */
  explicit Orthtree(Traits traits, const FT enlarge_ratio = 1.2) :
    m_traits(traits),
    m_node_points(m_node_properties.add_property<Node_data>("points")),
    m_node_depths(m_node_properties.add_property<std::uint8_t>("depths", 0)),
    m_node_coordinates(m_node_properties.add_property<Global_coordinates>("coordinates")),
    m_node_parents(m_node_properties.add_property<Maybe_node_index>("parents")),
    m_node_children(m_node_properties.add_property<Maybe_node_index>("children")) {

    m_node_properties.emplace();


    // init bbox with first values found
    auto [bbox_min, bbox_max] = m_traits.root_node_bbox_object()();

    // Dilate the bounding box
    Array bbox_centroid;
    FT max_length = FT(0);
    for (std::size_t i = 0; i < Dimension::value; ++i) {
      bbox_centroid[i] = (bbox_min[i] + bbox_max[i]) / FT(2);
      max_length = (std::max)(max_length, bbox_max[i] - bbox_min[i]);
    }
    max_length *= enlarge_ratio / FT(2);
    for (std::size_t i = 0; i < Dimension::value; ++i) {
      bbox_min[i] = bbox_centroid[i] - max_length;
      bbox_max[i] = bbox_centroid[i] + max_length;
    }

    auto construct_point_d_from_array
      = m_traits.construct_point_d_from_array_object();

    // save orthtree attributes
    m_bbox_min = construct_point_d_from_array(bbox_min);
    m_side_per_depth.push_back(bbox_max[0] - bbox_min[0]);
    data(root()) = m_traits.root_node_contents_object()();
  }

  /// @}

  /// \cond SKIP_IN_MANUAL

  // copy constructor
  Orthtree(const Orthtree& other) :
    m_traits(other.m_traits),
    m_bbox_min(other.m_bbox_min), m_side_per_depth(other.m_side_per_depth),
    m_node_properties(other.m_node_properties),
    m_node_points(m_node_properties.get_property<Node_data>("points")),
    m_node_depths(m_node_properties.get_property<std::uint8_t>("depths")),
    m_node_coordinates(m_node_properties.get_property<Global_coordinates>("coordinates")),
    m_node_parents(m_node_properties.get_property<Maybe_node_index>("parents")),
    m_node_children(m_node_properties.get_property<Maybe_node_index>("children")) {}

  // move constructor
  Orthtree(Orthtree&& other) :
    m_traits(other.m_traits),
    m_bbox_min(other.m_bbox_min), m_side_per_depth(other.m_side_per_depth),
    m_node_properties(std::move(other.m_node_properties)),
    m_node_points(m_node_properties.get_property<Node_data>("points")),
    m_node_depths(m_node_properties.get_property<std::uint8_t>("depths")),
    m_node_coordinates(m_node_properties.get_property<Global_coordinates>("coordinates")),
    m_node_parents(m_node_properties.get_property<Maybe_node_index>("parents")),
    m_node_children(m_node_properties.get_property<Maybe_node_index>("children")) {

    // todo: makes sure moved-from is still valid. Maybe this shouldn't be necessary.
    other.m_node_properties.emplace();
  }

  // Non-necessary but just to be clear on the rule of 5:

  // assignment operators deleted (PointRange is a ref)
  Orthtree& operator=(const Orthtree& other) = delete;

  Orthtree& operator=(Orthtree&& other) = delete;

  // move constructor
  /// \endcond

  /// \name Tree Building
  /// @{

  /*!
    \brief recursively subdivides the orthtree until it meets the given criteria.

    The split predicate is an `std::function` that takes a `Node_index` and an Orthtree reference, and
    returns a Boolean value (where `true` implies that the corresponding node needs to
    be split, `false` that the node should be a leaf).

    This function may be called several times with different
    predicates: in that case, nodes already split are left unaltered,
    while nodes that were not split and for which `split_predicate`
    returns `true` are split.

    \param split_predicate determines whether or not a node needs to
    be subdivided.
   */
  void refine(const Split_predicate& split_predicate) {

    // Initialize a queue of nodes that need to be refined
    std::queue<Node_index> todo;
    todo.push(0);

    // Process items in the queue until it's consumed fully
    while (!todo.empty()) {

      // Get the next element
      auto current = todo.front();
      todo.pop();

      // Check if this node needs to be processed
      if (split_predicate(current, *this)) {

        // Split the node, redistributing its points to its children
        split(current);

      }

      // Check if the node has children which need to be processed
      if (!is_leaf(current)) {

        // Process each of its children
        for (int i = 0; i < Degree::value; ++i)
          todo.push(child(current, i));
      }
    }
  }

  /*!

    \brief Convenience overload that refines an orthtree using a
    maximum depth and maximum number of points in a node as split
    predicate.

    This is equivalent to calling
    `refine(Orthtrees::Maximum_depth_and_maximum_number_of_inliers(max_depth,
    bucket_size))`.

    The refinement is stopped as soon as one of the conditions is
    violated: if a node has more inliers than `bucket_size` but is
    already at `max_depth`, it is not split. Similarly, a node that is
    at a depth smaller than `max_depth` but already has fewer inliers
    than `bucket_size`, it is not split.

    \param max_depth deepest a tree is allowed to be (nodes at this depth will not be split).
    \param bucket_size maximum points a node is allowed to contain.
   */
  void refine(size_t max_depth = 10, size_t bucket_size = 20) {
    refine(Orthtrees::Maximum_depth_and_maximum_number_of_inliers(max_depth, bucket_size));
  }

  /*!
    \brief refines the orthtree such that the difference of depth
    between two immediate neighbor leaves is never more than 1.
   */
  void grade() {

    // Collect all the leaf nodes
    std::queue<Node_index> leaf_nodes;
    for (Node_index leaf: traverse(Orthtrees::Leaves_traversal<Self>(*this))) {
      leaf_nodes.push(leaf);
    }

    // Iterate over the nodes
    while (!leaf_nodes.empty()) {

      // Get the next node
      Node_index node = leaf_nodes.front();
      leaf_nodes.pop();

      // Skip this node if it isn't a leaf anymore
      if (!is_leaf(node))
        continue;

      // Iterate over each of the neighbors
      for (int direction = 0; direction < 6; ++direction) {

        // Get the neighbor
        auto neighbor = adjacent_node(node, direction);

        // If it doesn't exist, skip it
        if (!neighbor)
          continue;

        // Skip if this neighbor is a direct sibling (it's guaranteed to be the same depth)
        // TODO: This check might be redundant, if it doesn't affect performance maybe I could remove it
        if (parent(*neighbor) == parent(node))
          continue;

        // If it's already been split, skip it
        if (!is_leaf(*neighbor))
          continue;

        // Check if the neighbor breaks our grading rule
        // TODO: could the rule be parametrized?
        if ((depth(node) - depth(*neighbor)) > 1) {

          // Split the neighbor
          split(*neighbor);

          // Add newly created children to the queue
          for (int i = 0; i < Degree::value; ++i) {
            leaf_nodes.push(child(*neighbor, i));
          }
        }
      }
    }
  }

  /// @}

  /// \name Accessors
  /// @{

  /*!
   * \brief Provides direct read-only access to the tree Traits.
   *
   * @return a const reference to the Traits instantiation.
   */
  const Traits& traits() const { return m_traits; }

  /*!
    \brief provides read-only access to the root node, and by
    extension the rest of the tree.

    \return a const reference to the root node of the tree.
   */
  Node_index root() const { return 0; }

  /*!
    \brief returns the deepest level reached by a leaf node in this tree (root being level 0).
   */
  std::size_t depth() const { return m_side_per_depth.size() - 1; }

  /*!
    \brief constructs a node index range using a tree-traversal function.

    This method allows iteration over the nodes of the tree with a
    user-selected order (preorder, postorder, leaves-only, etc.).

    \tparam Traversal model of `OrthtreeTraversal` that provides functions
    compatible with the type of the orthtree

    \param traversal the instance of `Traversal` used

    \return a forward input iterator over the node indices of the tree
   */
  template <typename Traversal>
  Node_index_range traverse(Traversal traversal) const {

    Node_index first = traversal.first_index();

    auto next = [=](const Self& tree, Node_index index) -> Maybe_node_index {
      return traversal.next_index(index);
    };

    return boost::make_iterator_range(Index_traversal_iterator<Self>(*this, first, next),
                                      Index_traversal_iterator<Self>());
  }


  /*!
    \brief Convenience method for using a traversal without constructing it yourself

    \tparam Traversal model of `OrthtreeTraversal` that provides functions
    compatible with the type of the orthtree

    \param args Arguments to to pass to the traversal's constructor, excluding the first (always an orthtree reference)

    \return a forward input iterator over the node indices of the tree
   */
  template <typename Traversal, typename ...Args>
  Node_index_range traverse(Args&& ...args) const {
    return traverse(Traversal{*this, std::forward<Args>(args)...});
  }

  /*!
    \brief constructs the bounding box of a node.

    \note The object constructed is not the bounding box of the point
    subset inside the node, but the bounding box of the node itself
    (cubic).
   */
  Bbox bbox(Node_index n) const {

    // Determine the side length of this node
    FT size = m_side_per_depth[depth(n)];

    // Determine the location this node should be split
    Array min_corner;
    Array max_corner;
    for (int i = 0; i < Dimension::value; i++) {

      min_corner[i] = m_bbox_min[i] + (global_coordinates(n)[i] * size);
      max_corner[i] = min_corner[i] + size;
    }

    // Create the bbox
    return {m_traits.construct_point_d_from_array_object()(min_corner),
            m_traits.construct_point_d_from_array_object()(max_corner)};
  }

  /// @}

  /// \name Custom Properties
  /// @{

  template <typename T>
  std::pair<std::reference_wrapper<Node_property_container::Array < T>>, bool>
  get_or_add_node_property(const std::string& name, const T default_value = T()) {
    return m_node_properties.get_or_add_property(name, default_value);
  }

  template <typename T>
  Node_property_container::Array <T>& add_node_property(const std::string& name, const T default_value = T()) {
    return m_node_properties.add_property(name, default_value);
  }

  template <typename T>
  Node_property_container::Array <T>& get_node_property(const std::string& name) {
    return m_node_properties.get_property<T>(name);
  }

  template <typename T>
  std::optional<std::reference_wrapper<Node_property_container::Array<T>>>
  get_node_property_if_exists(const std::string& name) {
    return m_node_properties.get_property_if_exists<T>(name);
  }

  // todo: is it ever useful to be able to delete/reset properties?

  /// @}

  /// \name Queries
  /// @{

  /*!
    \brief finds the leaf node which contains the point `p`.

    Traverses the orthtree and finds the deepest cell that has a
    domain enclosing the point passed. The point passed must be within
    the region enclosed by the orthtree (bbox of the root node).

    \param point query point.
    \return the node which contains the point.
   */
  const Node_index locate(const Point& point) const {

    // Make sure the point is enclosed by the orthtree
    CGAL_precondition (CGAL::do_intersect(point, bbox(root())));

    // Start at the root node
    Node_index node_for_point = root();

    // Descend the tree until reaching a leaf node
    while (!is_leaf(node_for_point)) {

      // Find the point to split around
      Point center = barycenter(node_for_point);

      // Find the index of the correct sub-node
      Local_coordinates local_coords;
      std::size_t dimension = 0;
      for (const auto& r: cartesian_range(center, point))
        local_coords[dimension++] = (get < 0 > (r) < get < 1 > (r));

      // Find the correct sub-node of the current node
      node_for_point = child(node_for_point, local_coords.to_ulong());
    }

    // Return the result
    return node_for_point;
  }

  /*!
    \brief finds the leaf nodes that intersect with any primitive.

    \note this function requires the function
    `bool CGAL::do_intersect(QueryType, Traits::Bbox_d)` to be defined.

    This function finds all the intersecting nodes and returns them as const pointers.

    \tparam Query the primitive class (e.g. sphere, ray)
    \tparam OutputIterator a model of `OutputIterator` that accepts `Node_index` types
    \param query the intersecting primitive.
    \param output output iterator.
   */
  template <typename Query, typename OutputIterator>
  OutputIterator intersected_nodes(const Query& query, OutputIterator output) const {
    return intersected_nodes_recursive(query, root(), output);
  }

  /// @}

  /// \name Operators
  /// @{

  /*!
    \brief compares the topology of the orthtree with that of `rhs`.

    Trees may be considered equivalent even if they contain different points.
    Equivalent trees must have the same bounding box and the same node structure.
   */
  bool operator==(const Self& rhs) const {

    // Identical trees should have the same bounding box
    if (rhs.m_bbox_min != m_bbox_min || rhs.m_side_per_depth[0] != m_side_per_depth[0])
      return false;

    // Identical trees should have the same depth
    if (rhs.depth() != depth())
      return false;

    // If all else is equal, recursively compare the trees themselves
    return is_topology_equal(*this, rhs);
  }

  /*!
    \brief compares the topology of the orthtree with that of `rhs`.
   */
  bool operator!=(const Self& rhs) const {
    return !operator==(rhs);
  }

  /// @}

  /// \name Node Access
  /// @{

  bool is_leaf(Node_index n) const {
    return !m_node_children[n].has_value();
  }

  bool is_root(Node_index n) const {
    return n == 0;
  }

  std::size_t depth(Node_index n) const {
//    std::cerr << n
//              << " " << m_node_depths.size()
//              << std::endl;
    return m_node_depths[n];
  }

  Node_data& data(Node_index n) {
    return m_node_points[n];
  }

  const Node_data& data(Node_index n) const {
    return m_node_points[n];
  }

  Global_coordinates global_coordinates(Node_index n) const {
    return m_node_coordinates[n];
  }

  Local_coordinates local_coordinates(Node_index n) const {
    Local_coordinates result;
    for (std::size_t i = 0; i < Dimension::value; ++i)
      result[i] = global_coordinates(n)[i] & 1;
    return result;
  }

  /*!
    \brief returns this node's parent.
    \pre `!is_root()`
   */
  Node_index parent(Node_index node) const {
    CGAL_precondition (!is_root(node));
    return *m_node_parents[node];
  }

  Node_index child(Node_index node, std::size_t i) const {
    CGAL_precondition (!is_leaf(node));
    return *m_node_children[node] + i;
  }

  Node_index descendant(Node_index node, std::size_t i) { return child(node, i); }

  template <typename... Indices>
  Node_index descendant(Node_index node, std::size_t i, Indices... remaining_indices) {
    return descendant(child(node, i), remaining_indices...);
  }

  template <typename... Indices>
  Node_index node(Indices... indices) {
    return descendant(root(), indices...);
  }

  const Maybe_node_index next_sibling(Node_index n) const {

    // Root node has no siblings
    if (is_root(n)) return {};

    // Find out which child this is
    std::size_t local_coords = local_coordinates(n).to_ulong(); // todo: add local_coordinates(n) helper

    // The last child has no more siblings
    if (int(local_coords) == Degree::value - 1)
      return {};

    // The next sibling is the child of the parent with the following local coordinates
    return child(parent(n), local_coords + 1);
  }

  const Maybe_node_index next_sibling_up(Node_index n) const {

    // the root node has no next sibling up
    if (n == 0) return {};

    auto up = Maybe_node_index{parent(n)};
    while (up) {

      if (next_sibling(*up)) return {next_sibling(*up)};

      up = is_root(*up) ? Maybe_node_index{} : Maybe_node_index{parent(*up)};
    }

    return {};
  }

  Node_index deepest_first_child(Node_index n) const {

    auto first = n;
    while (!is_leaf(first))
      first = child(first, 0);

    return first;
  }

  Maybe_node_index first_child_at_depth(Node_index n, std::size_t d) const {

    std::queue<Node_index> todo;
    todo.push(n);

    while (!todo.empty()) {
      Node_index node = todo.front();
      todo.pop();

      if (depth(node) == d)
        return node;

      if (!is_leaf(node))
        for (int i = 0; i < Degree::value; ++i)
          todo.push(child(node, i));
    }

    return {};
  }

  /*!
  \brief splits the node into subnodes.

  Only leaf nodes should be split.
  When a node is split it is no longer a leaf node.
  A number of `Degree::value` children are constructed automatically, and their values are set.
  Contents of this node are _not_ propagated automatically.
  It is the responsibility of the caller to redistribute the points contained by a node after splitting
 */
  void split(Node_index n) {

    // Make sure the node hasn't already been split
    CGAL_precondition (is_leaf(n));

    // Split the node to create children
    using Local_coordinates = Local_coordinates;
    m_node_children[n] = m_node_properties.emplace_group(Degree::value);
    for (std::size_t i = 0; i < Degree::value; i++) {

      Node_index c = *m_node_children[n] + i;

      // Make sure the node isn't one of its own children
      CGAL_assertion(n != *m_node_children[n] + i);

      Local_coordinates local_coordinates{i};
      for (int i = 0; i < Dimension::value; i++)
        m_node_coordinates[c][i] = (2 * m_node_coordinates[n][i]) + local_coordinates[i];
      m_node_depths[c] = m_node_depths[n] + 1;
      m_node_parents[c] = n;
    }

    // Check if we've reached a new max depth
    if (depth(n) + 1 == m_side_per_depth.size()) {
      // Update the side length map
      m_side_per_depth.push_back(*(m_side_per_depth.end() - 1) / 2);
    }

    // Find the point around which the node is split
    Point center = barycenter(n);

    // Add the node's points to its children
    m_traits.distribute_node_contents_object()(n, *this, center);
  }

  /*!
   * \brief eliminates this node's children, making it a leaf node.
   *
   * When a node is un-split, its children are automatically deleted.
   * After un-splitting a node it will be considered a leaf node.
   * Idempotent, un-splitting a leaf node has no effect.
   */
  void unsplit(Node_index n) {
    // todo: the child nodes should be de-allocated!
  }

  Point barycenter(Node_index n) const {

    // Determine the side length of this node
    FT size = m_side_per_depth[depth(n)];

    // Determine the location this node should be split
    Array bary;
    std::size_t i = 0;
    for (const FT& f: cartesian_range(m_bbox_min)) {
      bary[i] = FT(global_coordinates(n)[i]) * size + size / FT(2) + f;
      ++i;
    }

    // Convert that location into a point
    auto construct_point_d_from_array
      = m_traits.construct_point_d_from_array_object();
    return construct_point_d_from_array(bary);
  }

  static bool is_topology_equal(Node_index lhsNode, const Self& lhsTree, Node_index rhsNode, const Self& rhsTree) {

    // If one node is a leaf, and the other isn't, they're not the same
    if (lhsTree.is_leaf(lhsNode) != rhsTree.is_leaf(rhsNode))
      return false;

    // If both nodes are non-leaf
    if (!lhsTree.is_leaf(lhsNode)) {

      // Check all the children
      for (int i = 0; i < Degree::value; ++i) {
        // If any child cell is different, they're not the same
        if (!is_topology_equal(lhsTree.child(lhsNode, i), lhsTree,
                               rhsTree.child(rhsNode, i), rhsTree))
          return false;
      }
    }

    return (lhsTree.global_coordinates(lhsNode) == rhsTree.global_coordinates(rhsNode));
  }

  static bool is_topology_equal(const Self& lhs, const Self& rhs) {
    return is_topology_equal(lhs.root(), lhs, rhs.root(), rhs);
  }

  /*!
    \brief finds the directly adjacent node in a specific direction

    \pre `!is_null()`
    \pre `direction.to_ulong < 2 * Dimension::value`

    Adjacent nodes are found according to several properties:
    - adjacent nodes may be larger than the seek node, but never smaller
    - a node has at most `2 * Dimension::value` different adjacent nodes (in 3D: left, right, up, down, front, back)
    - adjacent nodes are not required to be leaf nodes

    Here's a diagram demonstrating the concept for a Quadtree:

    ```
    +---------------+---------------+
    |               |               |
    |               |               |
    |               |               |
    |       A       |               |
    |               |               |
    |               |               |
    |               |               |
    +-------+-------+---+---+-------+
    |       |       |   |   |       |
    |   A   |  (S)  +---A---+       |
    |       |       |   |   |       |
    +---+---+-------+---+---+-------+
    |   |   |       |       |       |
    +---+---+   A   |       |       |
    |   |   |       |       |       |
    +---+---+-------+-------+-------+
    ```

    - (S) : Seek node
    - A  : Adjacent node

    Note how the top adjacent node is larger than the seek node.  The
    right adjacent node is the same size, even though it contains
    further subdivisions.

    This implementation returns the adjacent node if it's found.  If
    there is no adjacent node in that direction, it returns a null
    node.

    \param direction which way to find the adjacent node relative to
    this one. Each successive bit selects the direction for the
    corresponding dimension: for an Octree in 3D, 010 means: negative
    direction in X, position direction in Y, negative direction in Z.

    \return the index of the adjacent node if it exists, nothing otherwise.
  */
  Maybe_node_index adjacent_node(Node_index n, Local_coordinates direction) const {

    // Direction:   LEFT  RIGHT  DOWN    UP  BACK FRONT
    // direction:    000    001   010   011   100   101

    // Nodes only have up to 2*dim different adjacent nodes (since cubes have 6 sides)
    CGAL_precondition(direction.to_ulong() < Dimension::value * 2);

    // The root node has no adjacent nodes!
    if (is_root(n)) return {};

    // The least significant bit indicates the sign (which side of the node)
    bool sign = direction[0];

    // The first two bits indicate the dimension/axis (x, y, z)
    uint8_t dimension = uint8_t((direction >> 1).to_ulong());

    // Create an offset so that the bit-significance lines up with the dimension (e.g. 1, 2, 4 --> 001, 010, 100)
    int8_t offset = (uint8_t) 1 << dimension;

    // Finally, apply the sign to the offset
    offset = (sign ? offset : -offset);

    // Check if this child has the opposite sign along the direction's axis
    if (local_coordinates(n)[dimension] != sign) {
      // This means the adjacent node is a direct sibling, the offset can be applied easily!
      return {child(parent(n), local_coordinates(n).to_ulong() + offset)};
    }

    // Find the parent's neighbor in that direction, if it exists
    auto adjacent_node_of_parent = adjacent_node(parent(n), direction);

    // If the parent has no neighbor, then this node doesn't have one
    if (!adjacent_node_of_parent) return {};

    // If the parent's adjacent node has no children, then it's this node's adjacent node
    if (is_leaf(*adjacent_node_of_parent))
      return adjacent_node_of_parent;

    // Return the nearest node of the parent by subtracting the offset instead of adding
    return {child(*adjacent_node_of_parent, local_coordinates(n).to_ulong() - offset)};
  }

  /*!
    \brief equivalent to `adjacent_node()`, with an adjacency direction rather than a bitset.
   */
  Maybe_node_index adjacent_node(Node_index n, Adjacency adjacency) const {
    return adjacent_node(n, std::bitset<Dimension::value>(static_cast<int>(adjacency)));
  }

  /// @}

private: // functions :

  bool do_intersect(Node_index n, const Sphere& sphere) const {

    // Create a cubic bounding box from the node
    Bbox node_cube = bbox(n);

    // Check for intersection between the node and the sphere
    return CGAL::do_intersect(node_cube, sphere);
  }

  template <typename Query, typename Node_output_iterator>
  Node_output_iterator intersected_nodes_recursive(const Query& query, Node_index node,
                                                   Node_output_iterator output) const {

    // Check if the current node intersects with the query
    if (CGAL::do_intersect(query, bbox(node))) {

      // if this node is a leaf, then it's considered an intersecting node
      if (is_leaf(node)) {
        // todo: output iterator should hold indices instead of pointers
        *output++ = node;
        return output;
      }

      // Otherwise, each of the children need to be checked
      for (int i = 0; i < Degree::value; ++i) {
        intersected_nodes_recursive(query, child(node, i), output);
      }
    }
    return output;
  }

public:

  /// \cond SKIP_IN_MANUAL
  void dump_to_polylines(std::ostream& os) const {
    for (const auto n: traverse<Orthtrees::Preorder_traversal>())
      if (is_leaf(n)) {
        Bbox box = bbox(n);
        dump_box_to_polylines(box, os);
      }
  }

  void dump_box_to_polylines(const Bbox_2& box, std::ostream& os) const {
    // dump in 3D for visualisation
    os << "5 "
       << box.xmin() << " " << box.ymin() << " 0 "
       << box.xmin() << " " << box.ymax() << " 0 "
       << box.xmax() << " " << box.ymax() << " 0 "
       << box.xmax() << " " << box.ymin() << " 0 "
       << box.xmin() << " " << box.ymin() << " 0" << std::endl;
  }

  void dump_box_to_polylines(const Bbox_3& box, std::ostream& os) const {
    // Back face
    os << "5 "
       << box.xmin() << " " << box.ymin() << " " << box.zmin() << " "
       << box.xmin() << " " << box.ymax() << " " << box.zmin() << " "
       << box.xmax() << " " << box.ymax() << " " << box.zmin() << " "
       << box.xmax() << " " << box.ymin() << " " << box.zmin() << " "
       << box.xmin() << " " << box.ymin() << " " << box.zmin() << std::endl;

    // Front face
    os << "5 "
       << box.xmin() << " " << box.ymin() << " " << box.zmax() << " "
       << box.xmin() << " " << box.ymax() << " " << box.zmax() << " "
       << box.xmax() << " " << box.ymax() << " " << box.zmax() << " "
       << box.xmax() << " " << box.ymin() << " " << box.zmax() << " "
       << box.xmin() << " " << box.ymin() << " " << box.zmax() << std::endl;

    // Traversal edges
    os << "2 "
       << box.xmin() << " " << box.ymin() << " " << box.zmin() << " "
       << box.xmin() << " " << box.ymin() << " " << box.zmax() << std::endl;
    os << "2 "
       << box.xmin() << " " << box.ymax() << " " << box.zmin() << " "
       << box.xmin() << " " << box.ymax() << " " << box.zmax() << std::endl;
    os << "2 "
       << box.xmax() << " " << box.ymin() << " " << box.zmin() << " "
       << box.xmax() << " " << box.ymin() << " " << box.zmax() << std::endl;
    os << "2 "
       << box.xmax() << " " << box.ymax() << " " << box.zmin() << " "
       << box.xmax() << " " << box.ymax() << " " << box.zmax() << std::endl;
  }

  std::string to_string(Node_index node) {
    std::stringstream stream;
    internal::print_orthtree_node(stream, node, *this);
    return stream.str();
  }

  friend std::ostream& operator<<(std::ostream& os, const Self& orthtree) {
    // Iterate over all nodes
    for (auto n: orthtree.traverse(Orthtrees::Preorder_traversal<Self>(orthtree))) {
      // Show the depth
      for (int i = 0; i < orthtree.depth(n); ++i)
        os << ". ";
      // Print the node
      internal::print_orthtree_node(os, n, orthtree);
      os << std::endl;
    }
    return os;
  }

  /// \endcond

}; // end class Orthtree

} // namespace CGAL

#endif // CGAL_ORTHTREE_H
