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
#include <CGAL/span.h>

#include <boost/function.hpp>
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
  decomposition of dD space for efficient access and
  computation.

  \details It builds a hierarchy of nodes which subdivices the space.
  Each node represents an axis-aligned hypercubic region of space.
  The contents of nodes depend on the Traits class, non-leaf nodes also
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
  using Traits = Traits_; ///< Geometry traits
  /// @}

  /// \name Traits Types
  /// @{
  using Dimension = typename Traits::Dimension; ///< Dimension of the tree
  using FT = typename Traits::FT; ///< Number type.
  using Point = typename Traits::Point_d; ///< Point type.
  using Bbox = typename Traits::Bbox_d; ///< Bounding box type.
  using Sphere = typename Traits::Sphere_d; ///< Sphere type.
  using Adjacency = typename Traits::Adjacency; ///< Adjacency type.

  using Node_data = typename Traits::Node_data;
  // todo: Node_data_element will only exist for certain Traits types, so I don't know if it can be re-exported

  /// @}

  /// \name Public Types
  /// @{

  /*!
   * \brief Self typedef for convenience.
   */
  using Self = Orthtree<Traits>;

  /*!
   * \brief Degree of the tree (number of children of non-leaf nodes).
   */
  using Degree = Dimension_tag<(2 << (Dimension::value - 1))>;

  /*!
   * \brief Index of a given node in the tree; the root always has index 0.
   */
  using Node_index = std::size_t;

  /*!
   * \brief Optional index of a node in the tree.
   */
  using Maybe_node_index = std::optional<Node_index>;

  // todo: maybe this could be private?
  using Node_property_container = Properties::Experimental::Property_container<Node_index>;

  /*!
    \brief Set of bits representing this node's relationship to its parent.

    Equivalent to an array of Booleans, where index[0] is whether `x`
    is greater, index[1] is whether `y` is greater, index[2] is whether
    `z` is greater, and so on for higher dimensions if needed.
    Used to represent a node's relationship to the center of its parent.
   */
  using Local_coordinates = std::bitset<Dimension::value>;

  /*!
    \brief Coordinates representing this node's relationship
    with the rest of the tree.

    Each value `(x, y, z, ...)` of global coordinates is calculated by doubling
    the parent's global coordinates and adding the local coordinates.
   */
  using Global_coordinates = std::array<std::uint32_t, Dimension::value>;

  /*!
   * \brief A predicate that determines whether a node must be split when refining a tree.
   */
  using Split_predicate = std::function<bool(Node_index, const Self&)>;

  /*!
   * \brief A model of `ConstRange` whose value type is `Node_index`.
   */
#ifdef DOXYGEN_RUNNING
  using Node_range = unspecified_type;
  using Node_index_range = unspecified_type;
#else
  using Node_index_range = boost::iterator_range<Index_traversal_iterator<Self>>;
#endif

  /// \cond SKIP_IN_MANUAL
  using Cartesian_ranges = Orthtrees::internal::Cartesian_ranges<Traits>;
  /// \endcond

  /// @}

private: // data members :

  Traits m_traits; /* the tree traits */

  Node_property_container m_node_properties;
  Node_property_container::Array <Node_data>& m_node_contents;
  Node_property_container::Array <std::uint8_t>& m_node_depths;
  Node_property_container::Array <Global_coordinates>& m_node_coordinates;
  Node_property_container::Array <Maybe_node_index>& m_node_parents;
  Node_property_container::Array <Maybe_node_index>& m_node_children;

  Point m_bbox_min;                  /* input bounding box min value */

  using Bbox_dimensions = std::array<FT, Dimension::value>;
  std::vector<Bbox_dimensions> m_side_per_depth;      /* side length per node's depth */

  Cartesian_ranges cartesian_range; /* a helper to easily iterate over coordinates of points */

public:

  /// \name Constructor
  /// @{

  /*!
    \brief creates an orthtree for a traits instance.

    The constructed orthtree has a root node with no children,
    containing the contents determined by `construct_root_node_contents_object` from the Traits class.
    That root node has a bounding box determined by `construct_root_node_bbox_object` from the Traits class,
    which typically encloses its contents.

    This single-node orthtree is valid and compatible
    with all Orthtree functionality, but any performance benefits are
    unlikely to be realized until `refine()` is called.

    \param traits the traits object.
  */
  explicit Orthtree(Traits traits) :
    m_traits(traits),
    m_node_contents(m_node_properties.add_property<Node_data>("contents")),
    m_node_depths(m_node_properties.add_property<std::uint8_t>("depths", 0)),
    m_node_coordinates(m_node_properties.add_property<Global_coordinates>("coordinates")),
    m_node_parents(m_node_properties.add_property<Maybe_node_index>("parents")),
    m_node_children(m_node_properties.add_property<Maybe_node_index>("children")) {

    m_node_properties.emplace();

    // init bbox with first values found
    auto bbox = m_traits.construct_root_node_bbox_object()();

    // Determine dimensions of the root bbox
    Bbox_dimensions size;
    for (int i = 0; i < Dimension::value; ++i)
      size[i] = (bbox.max)()[i] - (bbox.min)()[i];

    // save orthtree attributes
    m_bbox_min = (bbox.min)();
    m_side_per_depth.push_back(size);
    data(root()) = m_traits.construct_root_node_contents_object()();
  }

  /// @}

  /// \cond SKIP_IN_MANUAL

  // copy constructor
  Orthtree(const Orthtree& other) :
    m_traits(other.m_traits),
    m_bbox_min(other.m_bbox_min), m_side_per_depth(other.m_side_per_depth),
    m_node_properties(other.m_node_properties),
    m_node_contents(m_node_properties.get_property<Node_data>("contents")),
    m_node_depths(m_node_properties.get_property<std::uint8_t>("depths")),
    m_node_coordinates(m_node_properties.get_property<Global_coordinates>("coordinates")),
    m_node_parents(m_node_properties.get_property<Maybe_node_index>("parents")),
    m_node_children(m_node_properties.get_property<Maybe_node_index>("children")) {}

  // move constructor
  Orthtree(Orthtree&& other) :
    m_traits(other.m_traits),
    m_bbox_min(other.m_bbox_min), m_side_per_depth(other.m_side_per_depth),
    m_node_properties(std::move(other.m_node_properties)),
    m_node_contents(m_node_properties.get_property<Node_data>("contents")),
    m_node_depths(m_node_properties.get_property<std::uint8_t>("depths")),
    m_node_coordinates(m_node_properties.get_property<Global_coordinates>("coordinates")),
    m_node_parents(m_node_properties.get_property<Maybe_node_index>("parents")),
    m_node_children(m_node_properties.get_property<Maybe_node_index>("children")) {

    // todo: makes sure moved-from is still valid. Maybe this shouldn't be necessary.
    other.m_node_properties.emplace();
  }

  // Non-necessary but just to be clear on the rule of 5:

  // assignment operators deleted
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

        // Split the node, redistributing its contents to its children
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
    maximum depth and maximum number of inliers in a node as split
    predicate.

    This is equivalent to calling
    `refine(Orthtrees::Maximum_depth_and_maximum_number_of_inliers(max_depth,
    bucket_size))`.

    The refinement is stopped as soon as one of the conditions is
    violated: if a node has more inliers than `bucket_size` but is
    already at `max_depth`, it is not split. Similarly, a node that is
    at a depth smaller than `max_depth` but already has fewer inliers
    than `bucket_size`, it is not split.

    \warning This convenience method is only appropriate for trees with traits classes where
    `Node_data` is a list-like type with a `size()` method.

    \param max_depth deepest a tree is allowed to be (nodes at this depth will not be split).
    \param bucket_size maximum number of items a node is allowed to contain.
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

    auto next = [=](const Self&, Node_index index) -> Maybe_node_index {
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

    \note The object constructed is not the bounding box of the node's contents,
    but the bounding box of the node itself.
    For a cubic orthtree, this will always be cubic.
   */
  Bbox bbox(Node_index n) const {

    using Cartesian_coordinate = std::array<FT, Dimension::value>;
    Cartesian_coordinate min_corner, max_corner;
    Bbox_dimensions size = m_side_per_depth[depth(n)];
    for (int i = 0; i < Dimension::value; i++) {
      min_corner[i] = m_bbox_min[i] + (global_coordinates(n)[i] * size[i]);
      max_corner[i] = m_bbox_min[i] + ((global_coordinates(n)[i] + 1) * size[i]);
    }
    return {std::apply(m_traits.construct_point_d_object(), min_corner),
            std::apply(m_traits.construct_point_d_object(), max_corner)};
  }

  /// @}

  /// \name Custom Properties
  /// @{

  /*!
   * \brief pass-through to `get_or_add_property` for node properties.
   */
  template <typename T>
  std::pair<std::reference_wrapper<Node_property_container::Array < T>>, bool>
  get_or_add_node_property(const std::string& name, const T default_value = T()) {
    return m_node_properties.get_or_add_property(name, default_value);
  }

  /*!
   * \brief pass-through to `add_property` for node properties.
   */
  template <typename T>
  Node_property_container::Array <T>& add_node_property(const std::string& name, const T default_value = T()) {
    return m_node_properties.add_property(name, default_value);
  }

  /*!
   * \brief pass-through to `get_property` for node properties.
   */
  template <typename T>
  Node_property_container::Array <T>& get_node_property(const std::string& name) {
    return m_node_properties.get_property<T>(name);
  }

  /*!
   * \brief pass-through to `get_property_if_exists` for node properties.
   */
  template <typename T>
  std::optional<std::reference_wrapper<Node_property_container::Array < T>>>
  get_node_property_if_exists(const std::string& name) {
    return m_node_properties.get_property_if_exists<T>(name);
  }

  // todo: is it ever useful to be able to delete/reset properties?

  /// @}

  /// \name Queries
  /// @{

  /*!
    \brief finds the leaf node which contains a particular point in space.

    Traverses the orthtree and finds the deepest cell that has a
    domain enclosing the point passed. The point passed must be within
    the region enclosed by the orthtree (bbox of the root node).

    \param point query point.
    \return the index of the node which contains the point.
   */
    Node_index locate(const Point& point) const {

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

    This function finds all the intersecting nodes and writes their indices to the ouput iterator.

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

    Trees may be considered equivalent even if they have different contents.
    Equivalent trees must have the same root bounding box and the same node structure.
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

  /*!
   * \brief Determines whether the node specified by index `n` is a leaf node.
   *
   * @param n index of the node to check.
   * @return true of the node is a leaf, false otherwise.
   */
  bool is_leaf(Node_index n) const {
    return !m_node_children[n].has_value();
  }

  /*!
   * \brief Determines whether the node specified by index `n` is a root node.
   *
   * @param n index of the node to check.
   * @return true of the node is a root, false otherwise.
   */
  bool is_root(Node_index n) const {
    return n == 0;
  }

  /*!
   * \brief Determines the depth of the node specified.
   *
   * The root node has depth 0, its children have depth 1, and so on.
   *
   * @param n index of the node to check.
   * @return the depth of the node within its tree.
   */
  std::size_t depth(Node_index n) const {
    return m_node_depths[n];
  }

  /*!
   * \brief Retrieves a reference to the Node_data associated with the node specified by `n`.
   *
   * @param n index of the node to retrieve data for.
   * @return the data associated with the node.
   */
  Node_data& data(Node_index n) {
    return m_node_contents[n];
  }

  /*!
   * \brief const version of `data()`
   */
  const Node_data& data(Node_index n) const {
    return m_node_contents[n];
  }

  /*!
   * \brief Retrieves the global coordinates of the node.
   */
  Global_coordinates global_coordinates(Node_index n) const {
    return m_node_coordinates[n];
  }

  /*!
   * \brief Retrieves the local coordinates of the node.
   */
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

  /*!
    \brief returns this node's `i`th child.
    \pre `!is_leaf()`
   */
  Node_index child(Node_index node, std::size_t i) const {
    CGAL_precondition (!is_leaf(node));
    return *m_node_children[node] + i;
  }

  /*!
   * \brief Retrieves an arbitrary descendant of the node specified by `node`.
   *
   * Convenience function to avoid the need to call `orthtree.child(orthtree.child(node, 0), 1)`.
   *
   * Each index in `indices` specifies which child to enter as descending the tree from `node` down.
   * Indices are evaluated in the order they appear as parameters, so
   * `descendant(root, 0, 1)` returns the second child of the first child of the root.
   *
   * @param node the node to descend
   * @param indices the descent to perform
   * @return the index of the specified descendant node.
   */
  template <typename... Indices>
  Node_index descendant(Node_index node, Indices... indices) {
    return recursive_descendant(node, indices...);
  }

  /*!
   * \brief Convenience function for retrieving arbitrary nodes, equivalent to `tree.descendant(tree.root(), indices...)`.
   */
  template <typename... Indices>
  Node_index node(Indices... indices) {
    return descendant(root(), indices...);
  }

  /*!
   * \brief Finds the next sibling in the parent of the node specified by the index `n`.
   *
   * Traverses the tree in increasing order of local index (e.g. 000, 001, 010, etc.)
   *
   * @param n the node to find the sibling of.
   * @return the next sibling of `n` if `n` is not the last node in its parent, otherwise nothing.
   */
  const Maybe_node_index next_sibling(Node_index n) const {

    // Root node has no siblings
    if (is_root(n)) return {};

    // Find out which child this is
    std::size_t local_coords = local_coordinates(n).to_ulong();

    // The last child has no more siblings
    if (int(local_coords) == Degree::value - 1)
      return {};

    // The next sibling is the child of the parent with the following local coordinates
    return child(parent(n), local_coords + 1);
  }

  /*!
   * \brief Finds the next sibling of the parent of the node specified by `n` if it exists.
   *
   * @param n the node to find the sibling up of.
   * @return The next sibling of the parent of `n` if `n` is not the root and its parent has a sibling, otherwise nothing.
   */
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

  /*!
   * \brief Finds the leaf node reached when descending the tree and always choosing child 0.
   *
   * This is the starting point of a depth-first traversal.
   *
   * @param n the node to find the deepest first child of.
   * @return the index of the deepest first child.
   */
  Node_index deepest_first_child(Node_index n) const {

    auto first = n;
    while (!is_leaf(first))
      first = child(first, 0);

    return first;
  }

  /*!
   * \brief Finds node reached when descending the tree to a depth `d` and always choosing child 0.
   *
   * Similar to `deepest_first_child`, but does not go to a fixed depth.
   *
   * @param n the node to find the `d`th first child of.
   * @param d the depth to descend to.
   * @return the index of the `d`th first child, nothing if the tree is not deep enough.
   */
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
  Contents of this node are _not_ propagated automatically, this is responsibility of the
  `distribute_node_contents_object` in the Traits class.
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
      // Update the side length map with the dimensions of the children
      Bbox_dimensions size = m_side_per_depth.back();
      Bbox_dimensions child_size;
      for (int i = 0; i < Dimension::value; ++i)
        child_size[i] = size[i] / FT(2);
      m_side_per_depth.push_back(child_size);
    }

    // Find the point around which the node is split
    Point center = barycenter(n);

    // Add the node's contents to its children
    m_traits.distribute_node_contents_object()(n, *this, center);
  }


  /*!
   * \brief Finds the center point of the node specified by `n`.
   *
   * @param n The node to find the center point for.
   * @return the center point of node `n`.
   */
  Point barycenter(Node_index n) const {

    // Determine the side length of this node
    Bbox_dimensions size = m_side_per_depth[depth(n)];

    // Determine the location this node should be split
    Bbox_dimensions bary;
    std::size_t i = 0;
    for (const FT& f: cartesian_range(m_bbox_min)) {
      bary[i] = FT(global_coordinates(n)[i]) * size[i] + size[i] / FT(2) + f;
      ++i;
    }

    // Convert that location into a point
    return std::apply(m_traits.construct_point_d_object(), bary);
  }

  /*!
   * \brief Determines whether a pair of subtrees have the same topology.
   *
   * @param lhsNode a node in lhsTree
   * @param lhsTree an Orthtree
   * @param rhsNode a node in rhsTree
   * @param rhsTree another Orthtree
   * @return true if lhsNode and rhsNode have the same topology, false otherwise
   */
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

  /*!
   * \brief Helper function for calling `is_topology_equal` on the root nodes of two trees.
   */
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

    \param n
    query node.

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

  Node_index recursive_descendant(Node_index node, std::size_t i) { return child(node, i); }

  template <typename... Indices>
  Node_index recursive_descendant(Node_index node, std::size_t i, Indices... remaining_indices) {
    return recursive_descendant(child(node, i), remaining_indices...);
  }

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
