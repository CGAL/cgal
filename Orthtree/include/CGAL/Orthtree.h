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

#include <CGAL/property_map.h>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Dimension.h>

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
  \tparam PointRange_ must be a model of range whose value type is the key type of `PointMap_`
  \tparam PointMap_ must be a model of `ReadablePropertyMap` whose value type is `Traits_::Point_d`
 */
template<typename Traits_, typename PointRange_,
         typename PointMap_ = Identity_property_map<typename Traits_::Point_d> >
class Orthtree
{

public:

  /// \name Template Types
  /// @{
  typedef Traits_ Traits; ///< Geometry traits
  typedef PointRange_ PointRange; ///< Point range
  typedef PointMap_ PointMap; ///< Point map
  /// @}

  /// \name Traits Types
  /// @{
  typedef typename Traits::Dimension Dimension; ///< Dimension of the tree
  typedef typename Traits::FT FT; ///< Number type.
  typedef typename Traits::Point_d Point; ///< Point type.
  typedef typename Traits::Bbox_d Bbox; ///< Bounding box type.
  typedef typename Traits::Sphere_d Sphere; ///< Sphere type.

  /// \cond SKIP_IN_MANUAL
  typedef typename Traits::Array Array; ///< Array type.
  typedef typename Traits::Construct_point_d_from_array
  Construct_point_d_from_array;
  typedef typename Traits::Construct_bbox_d
  Construct_bbox_d;
  /// \endcond
  /// @}

  /// \name Public Types
  /// @{

  /*!
   * \brief Self typedef for convenience.
   */
  typedef Orthtree<Traits, PointRange, PointMap> Self;

  /*!
   * \brief Degree of the tree (number of children of non-leaf nodes).
   */
  typedef Dimension_tag<(2 << (Dimension::value-1))> Degree;

  /*!
   * \brief The Sub-tree / Orthant type.
   */
  class Node;

  /*!
   * \brief A predicate that determines whether a node must be split when refining a tree.
   */
  typedef std::function<bool(Node)> Split_predicate;

  /*!
   * \brief A model of `ConstRange` whose value type is `Node`.
   */
#ifdef DOXYGEN_RUNNING
  typedef unspecified_type Node_range;
#else
  typedef boost::iterator_range<Traversal_iterator<Node> > Node_range;
#endif

  /// \cond SKIP_IN_MANUAL

  /*!
   * \brief A function that determines the next node in a traversal given the current one.
   */
  typedef std::function<Node(Node)> Node_traversal_method_const;

  /// \endcond

  /// \cond SKIP_IN_MANUAL
  typedef typename PointRange::iterator Range_iterator;
  typedef typename std::iterator_traits<Range_iterator>::value_type Range_type;
  typedef Orthtrees::internal::Cartesian_ranges<Traits> Cartesian_ranges;
  /// \endcond

  /// @}

private: // data members :

  Traits m_traits; /* the tree traits */
  PointRange& m_range;              /* input point range */
  PointMap m_point_map;          /* property map: `value_type of InputIterator` -> `Point` (Position) */

  Node m_root;                      /* root node of the orthtree */

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
  Orthtree(PointRange& point_range,
           PointMap point_map = PointMap(),
           const FT enlarge_ratio = 1.2,
           Traits traits = Traits())
    : m_traits (traits)
    , m_range (point_range)
    , m_point_map (point_map)
    , m_root(Node(), 0)
  {
    Array bbox_min;
    Array bbox_max;

    // init bbox with first values found
    {
      const Point& p = get (m_point_map, *(point_range.begin()));
      std::size_t i = 0;
      for (const FT& x : cartesian_range(p))
      {
        bbox_min[i] = x;
        bbox_max[i] = x;
        ++ i;
      }
    }

    for (const Range_type& r : point_range)
    {
      const Point& p = get (m_point_map, r);
      std::size_t i = 0;
      for (const FT& x : cartesian_range(p))
      {
        bbox_min[i] = (std::min)(x, bbox_min[i]);
        bbox_max[i] = (std::max)(x, bbox_max[i]);
        ++ i;
      }
    }

    Array bbox_centroid;
    FT max_length = FT(0);
    for (std::size_t i = 0 ; i < Dimension::value; ++ i)
    {
      bbox_centroid[i] = (bbox_min[i] + bbox_max[i]) / FT(2);
      max_length = (std::max)(max_length, bbox_max[i] - bbox_min[i]);
    }
    max_length *= enlarge_ratio / FT(2);

    for (std::size_t i = 0 ; i < Dimension::value; ++ i)
    {
      bbox_min[i] = bbox_centroid[i] - max_length;
      bbox_max[i] = bbox_centroid[i] + max_length;
    }

    Construct_point_d_from_array construct_point_d_from_array
      = m_traits.construct_point_d_from_array_object();

    // save orthtree attributes
    m_bbox_min = construct_point_d_from_array(bbox_min);
    m_side_per_depth.push_back(bbox_max[0] - bbox_min[0]);
    m_root.points() = {point_range.begin(), point_range.end()};
  }

  /// @}

  /// \cond SKIP_IN_MANUAL

  // copy constructor
  Orthtree (const Orthtree& other)
    : m_traits (other.m_traits)
    , m_range (other.m_range)
    , m_point_map (other.m_point_map)
    , m_root (other.m_root.deep_copy())
    , m_bbox_min (other.m_bbox_min)
    , m_side_per_depth(other.m_side_per_depth)
  { }

  // move constructor
  Orthtree (Orthtree&& other)
    : m_traits (other.m_traits)
    , m_range (other.m_range)
    , m_point_map (other.m_point_map)
    , m_root (other.m_root)
    , m_bbox_min (other.m_bbox_min)
    , m_side_per_depth(other.m_side_per_depth)
  {
    other.m_root = Node(Node(), 0);
  }

  // Non-necessary but just to be clear on the rule of 5:

  // assignement operators deleted (PointRange is a ref)
  Orthtree& operator= (const Orthtree& other) = delete;
  Orthtree& operator= (Orthtree&& other) = delete;
  // Destructor
  ~Orthtree()
  {
    std::queue<Node> nodes;
    nodes.push(m_root);
    while (!nodes.empty())
    {
      Node node = nodes.front();
      nodes.pop();
      if (!node.is_leaf())
        for (std::size_t i = 0; i < Degree::value; ++ i)
          nodes.push (node[i]);
      node.free();
    }
  }

  // move constructor
  /// \endcond

  /// \name Tree Building
  /// @{

  /*!
    \brief recursively subdivides the orthtree until it meets the given criteria.

    The split predicate is a `std::function` that takes a `Node` and
    returns a Boolean value (where `true` implies that a `Node` needs to
    be split, `false` that the `Node` should be a leaf).

    This function may be called several times with different
    predicates: in that case, nodes already split are left unaltered,
    while nodes that were not split and for which `split_predicate`
    returns `true` are split.

    \param split_predicate determines whether or not a node needs to
    be subdivided.
   */
  void refine(const Split_predicate& split_predicate) {

    // If the tree has already been refined, reset it
    if (!m_root.is_leaf())
      m_root.unsplit();

    // Reset the side length map, too
    m_side_per_depth.resize(1);

    // Initialize a queue of nodes that need to be refined
    std::queue<Node> todo;
    todo.push(m_root);

    // Process items in the queue until it's consumed fully
    while (!todo.empty()) {

      // Get the next element
      Node current = todo.front();
      todo.pop();

      // Check if this node needs to be processed
      if (split_predicate(current)) {

        // Check if we've reached a new max depth
        if (current.depth() == depth()) {

          // Update the side length map
          m_side_per_depth.push_back(*(m_side_per_depth.end() - 1) / 2);
        }

        // Split the node, redistributing its points to its children
        split(current);

        // Process each of its children
        for (int i = 0; i < Degree::value; ++i)
          todo.push(current[i]);

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
    std::queue<Node> leaf_nodes;
    for (Node leaf : traverse(Orthtrees::Leaves_traversal())) {
      // TODO: I'd like to find a better (safer) way of doing this
      leaf_nodes.push(leaf);
    }

    // Iterate over the nodes
    while (!leaf_nodes.empty()) {

      // Get the next node
      Node node = leaf_nodes.front();
      leaf_nodes.pop();

      // Skip this node if it isn't a leaf anymore
      if (!node.is_leaf())
        continue;

      // Iterate over each of the neighbors
      for (int direction = 0; direction < 6; ++direction) {

        // Get the neighbor
        Node neighbor = node.adjacent_node(direction);

        // If it doesn't exist, skip it
        if (neighbor.is_null())
          continue;

        // Skip if this neighbor is a direct sibling (it's guaranteed to be the same depth)
        // TODO: This check might be redundant, if it doesn't affect performance maybe I could remove it
        if (neighbor.parent() == node.parent())
          continue;

        // If it's already been split, skip it
        if (!neighbor.is_leaf())
          continue;

        // Check if the neighbor breaks our grading rule
        // TODO: could the rule be parametrized?
        if ((node.depth() - neighbor.depth()) > 1) {

          // Split the neighbor
          split(neighbor);

          // Add newly created children to the queue
          for (int i = 0; i < Degree::value; ++i) {
            leaf_nodes.push(neighbor[i]);
          }
        }
      }
    }
  }

  /// @}

  /// \name Accessors
  /// @{

  /*!
    \brief returns the root node.
   */
  Node root() const { return m_root; }

  /*!
    \brief Convenience function to access the child nodes of the root
    node by their indices.

    `my_tree[5]` is equivalent to `my_tree.root()[5]`.

    \sa `Node::operator[]()`

    \param index the index of the child node.
    \return the accessed node.
   */
  Node operator[](std::size_t index) const { return m_root[index]; }

  /*!
    \brief returns the deepest level reached by a leaf node in this tree (root being level 0).
   */
  std::size_t depth() const { return m_side_per_depth.size() - 1; }

  /*!
    \brief constructs a node range using a tree-traversal function.

    This method allows to iterate on the nodes of the tree with a
    user-selected order (preorder, postorder, leaves-only, etc.).

    \tparam Traversal model of `OrthtreeTraversal` that provides functions
    compatible with the type of the orthree

    \param traversal the instance of `Traversal` used

    \return a forward input iterator over the nodes of the tree
   */
  template<typename Traversal>
  Node_range traverse(const Traversal &traversal = Traversal()) const {

    Node first = traversal.first(m_root);

    Node_traversal_method_const next
      = [&](const Node& n) -> Node { return traversal.next(n); };

    return boost::make_iterator_range(Traversal_iterator<Node>(first, next),
                                      Traversal_iterator<Node>());
  }

  /*!
    \brief constructs the bounding box of a node.

    \note The object constructed is not the bounding box of the point
    subset inside the node, but the bounding box of the node itself
    (cubic).
   */
  Bbox bbox(const Node &node) const {
    // Determine the side length of this node
    FT size = m_side_per_depth[node.depth()];

    // Determine the location this node should be split
    Array min_corner;
    Array max_corner;
    for (int i = 0; i < Dimension::value; i++) {

      min_corner[i] = m_bbox_min[i] + (node.global_coordinates()[i] * size);
      max_corner[i] = min_corner[i] + size;
    }

    // Create the bbox
    Construct_bbox_d construct_bbox
      = m_traits.construct_bbox_d_object();
    return construct_bbox(min_corner, max_corner);
  }

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
  Node locate(const Point &point) const {

    // Make sure the point is enclosed by the orthtree
    CGAL_precondition (CGAL::do_intersect(point, bbox(m_root)));

    // Start at the root node
    auto node_for_point = m_root;

    // Descend the tree until reaching a leaf node
    while (!node_for_point.is_leaf()) {

      // Find the point to split around
      Point center = barycenter(node_for_point);

      // Find the index of the correct sub-node
      typename Node::Local_coordinates index;
      std::size_t dimension = 0;
      for (const auto& r : cartesian_range(center, point))
        index[dimension ++] = (get<0>(r) < get<1>(r));

      // Find the correct sub-node of the current node
      node_for_point = node_for_point[index.to_ulong()];
    }

    // Return the result
    return node_for_point;
  }

  /*!
    \brief finds the `k` nearest neighbors of `query`.

    Nearest neighbors are outputted in order of increasing distance to
    `query`.

    \tparam OutputIterator a model of `OutputIterator` that accept `Point_d` objects.
    \param query a query point.
    \param k the number of neighbors.
    \param output the output iterator.
   */
  template<typename OutputIterator>
  OutputIterator nearest_neighbors (const Point& query,
                                    std::size_t k,
                                    OutputIterator output) const {
    Sphere query_sphere (query, (std::numeric_limits<FT>::max)());
    return nearest_k_neighbors_in_radius(query_sphere, k, output);
  }

  /*!
    \brief finds the points in `sphere`.

    Nearest neighbors are outputted in order of increasing distance to
    the center of `sphere`.

    \tparam OutputIterator a model of `OutputIterator` that accept `Point_d` objects.
    \param query query sphere.
    \param output output iterator.
   */
  template<typename OutputIterator>
  OutputIterator nearest_neighbors (const Sphere& query, OutputIterator output) const {
    Sphere query_sphere = query;
    return nearest_k_neighbors_in_radius(query_sphere,
                                         (std::numeric_limits<std::size_t>::max)(), output);
  }

  /*!
    \brief finds the leaf nodes that intersect with any primitive.

    \note this function requires the function
    `bool CGAL::do_intersect(QueryType, Traits::Bbox_d)` to be defined.

    This function finds all the intersecting nodes and returns them as const pointers.

    \tparam Query the primitive class (e.g. sphere, ray)
    \tparam OutputIterator a model of `OutputIterator` that accepts `Node` objects
    \param query the intersecting primitive.
    \param output output iterator.
   */
  template<typename Query, typename OutputIterator>
  OutputIterator intersected_nodes (const Query& query, OutputIterator output) const {
    return intersected_nodes_recursive(query, root(), output);
  }

  /// @}

  /// \name Operators
  /// @{

  /*!
    \brief compares the topology of the orthtree with that of `rhs`.

    Trees may be considered equivalent even if they contain different points.
    Equivalent trees must have the same bounding box and the same node structure.
    Node structure is evaluated by comparing the root nodes using the node equality operator.
   */
  bool operator==(const Self &rhs) const {

    // Identical trees should have the same bounding box
    if (rhs.m_bbox_min != m_bbox_min || rhs.m_side_per_depth[0] != m_side_per_depth[0])
      return false;

    // Identical trees should have the same depth
    if (rhs.depth() != depth())
      return false;

    // If all else is equal, recursively compare the trees themselves
    return Node::is_topology_equal(rhs.m_root, m_root);
  }

  /*!
    \brief compares the topology of the orthtree with that of `rhs`.
   */
  bool operator!=(const Self &rhs) const {
    return !operator==(rhs);
  }

  /// @}

  // TODO: Document this
  // TODO: Could this method name be reduced to just "center" ?
  Point barycenter(const Node& node) const {

    // Determine the side length of this node
    FT size = m_side_per_depth[node.depth()];

    // Determine the location this node should be split
    Array bary;
    std::size_t i = 0;
    for (const FT& f : cartesian_range(m_bbox_min))
    {
      bary[i] = FT(node.global_coordinates()[i]) * size + size / FT(2) + f;
      ++ i;
    }

    // Convert that location into a point
    Construct_point_d_from_array construct_point_d_from_array
      = m_traits.construct_point_d_from_array_object();
    return construct_point_d_from_array(bary);
  }

private: // functions :

  void reassign_points(Node &node, Range_iterator begin, Range_iterator end, const Point &center,
                       std::bitset<Dimension::value> coord = {},
                       std::size_t dimension = 0) {

    // Root case: reached the last dimension
    if (dimension == Dimension::value) {

      node[coord.to_ulong()].points() = {begin, end};

      return;
    }

    // Split the point collection around the center point on this dimension
    Range_iterator split_point = std::partition
      (begin, end,
       [&](const Range_type &a) -> bool {
        // This should be done with cartesian iterator but it seems
        // complicated to do efficiently
         return (get(m_point_map, a)[int(dimension)] < center[int(dimension)]);
      });

    // Further subdivide the first side of the split
    std::bitset<Dimension::value> coord_left = coord;
    coord_left[dimension] = false;
    reassign_points(node, begin, split_point, center, coord_left, dimension + 1);

    // Further subdivide the second side of the split
    std::bitset<Dimension::value> coord_right = coord;
    coord_right[dimension] = true;
    reassign_points(node, split_point, end, center, coord_right, dimension + 1);

  }

  void split(Node& node) {

    // Make sure the node hasn't already been split
    CGAL_precondition (node.is_leaf());

    // Split the node to create children
    node.split();

    // Find the point to around which the node is split
    Point center = barycenter(node);

    // Add the node's points to its children
    reassign_points(node, node.points().begin(), node.points().end(), center);
  }

  bool do_intersect(const Node &node, const Sphere &sphere) const {

    // Create a cubic bounding box from the node
    Bbox node_cube = bbox(node);

    // Check for overlap between the node's box and the sphere as a box, to quickly catch some cases
    // FIXME: Activating this causes slower times
//    if (!do_overlap(node_cube, sphere.bbox()))
//      return false;

    // Check for intersection between the node and the sphere
    return CGAL::do_intersect(node_cube, sphere);
  }

  // TODO: There has to be a better way than using structs like these!
  struct Point_with_distance {
    Point point;
    FT distance;
  };
  struct Node_index_with_distance {
    typename Node::Local_coordinates index;
    FT distance;

    Node_index_with_distance (const typename Node::Local_coordinates& index,
                              const FT& distance)
      : index(index), distance(distance)
    { }
  };

  void nearest_k_neighbors_recursive(Sphere& search_bounds, const Node &node,
                                     std::vector<Point_with_distance> &results, FT epsilon = 0) const {

    // Check whether the node has children
    if (node.is_leaf()) {

      // Base case: the node has no children

      // Loop through each of the points contained by the node
      // Note: there might be none, and that should be fine!
      for (auto point_index : node.points()) {

        // Retrieve each point from the orthtree's point map
        auto point = get(m_point_map, point_index);

        // Pair that point with its distance from the search point
        Point_with_distance current_point_with_distance =
                {point, squared_distance(point, search_bounds.center())};

        // Check if the new point is within the bounds
        if (current_point_with_distance.distance < search_bounds.squared_radius()) {

          // Check if the results list is full
          if (results.size() == results.capacity()) {

            // Delete a point if we need to make room
            results.pop_back();
          }

          // Add the new point
          results.push_back(current_point_with_distance);

          // Sort the list
          std::sort(results.begin(), results.end(), [=](auto &left, auto &right) {
            return left.distance < right.distance;
          });

          // Check if the results list is full
          if (results.size() == results.capacity()) {

            // Set the search radius
            search_bounds = Sphere(search_bounds.center(), results.back().distance + epsilon);
          }
        }
      }
    } else {

      // Recursive case: the node has children

      // Create a list to map children to their distances
      std::vector<Node_index_with_distance> children_with_distances;
      children_with_distances.reserve(Degree::value);

      // Fill the list with child nodes
      for (int index = 0; index < Degree::value; ++index) {
        Node child_node = node[index];

        // Add a child to the list, with its distance
        children_with_distances.emplace_back(typename Node::Local_coordinates(index),
                                             CGAL::squared_distance(search_bounds.center(), barycenter(child_node)));
      }

      // Sort the children by their distance from the search point
      std::sort(children_with_distances.begin(), children_with_distances.end(), [=](auto &left, auto &right) {
        return left.distance < right.distance;
      });

      // Loop over the children
      for (auto child_with_distance : children_with_distances) {
        Node child_node = node[child_with_distance.index.to_ulong()];

        // Check whether the bounding box of the child intersects with the search bounds
        if (do_intersect(child_node, search_bounds)) {

          // Recursively invoke this function
          nearest_k_neighbors_recursive(search_bounds, child_node, results);
        }
      }
    }
  }

  template<typename Query, typename Node_output_iterator>
  Node_output_iterator intersected_nodes_recursive(const Query &query, const Node &node,
                                                   Node_output_iterator output) const {

    // Check if the current node intersects with the query
    if (CGAL::do_intersect(query, bbox(node))) {

      // if this node is a leaf, than it's considered an intersecting node
      if (node.is_leaf()) {
        *output++ = node;
        return output;
      }

      // Otherwise, each of the children need to be checked
      for (int i = 0; i < Degree::value; ++i) {
        intersected_nodes_recursive(query, node[i], output);
      }
    }
    return output;
  }

  /*!
    \brief finds the `k` points within a specific radius that are nearest to `query`.

    This function guarantees that there are no closer points than the ones returned,
    but it does not guarantee that it will return at least `k` points.
    For a query where the search radius encloses `k` or fewer points, all enclosed points will be returned.
    If the search radius passed is too small, no points may be returned.
    This function is useful when the user already knows how sparse the points are,
    or if they do not care about points that are too far away.
    Setting a small radius may have performance benefits.

    \tparam OutputIterator must be a model of `OutputIterator` that accepts points
    \param search_point the location to find points near
    \param search_radius_squared the size of the region to search within
    \param k the number of points to find
    \param output the output iterator to add the found points to (in order of increasing distance)
   */
  template<typename OutputIterator>
  OutputIterator nearest_k_neighbors_in_radius
  (Sphere& query_sphere,
   std::size_t k, OutputIterator output) const {

    // Create an empty list of points
    std::vector<Point_with_distance> points_list;
    if (k != (std::numeric_limits<std::size_t>::max)())
      points_list.reserve(k);

    // Invoking the recursive function adds those points to the vector (passed by reference)
    nearest_k_neighbors_recursive(query_sphere, m_root, points_list);

    // Add all the points found to the output
    for (auto &item : points_list)
      *output++ = item.point;

    return output;
  }

public:

  /// \cond SKIP_IN_MANUAL
  void dump_to_polylines (std::ostream& os) const
  {
    for (const Node& n : traverse<Orthtrees::Preorder_traversal>())
      if (n.is_leaf())
      {
        Bbox box = bbox(n);
        dump_box_to_polylines (box, os);
      }
  }

  void dump_box_to_polylines (const Bbox_2& box, std::ostream& os) const
  {
    // dump in 3D for visualisation
    os << "5 "
       << box.xmin() << " " << box.ymin() << " 0 "
       << box.xmin() << " " << box.ymax() << " 0 "
       << box.xmax() << " " << box.ymax() << " 0 "
       << box.xmax() << " " << box.ymin() << " 0 "
       << box.xmin() << " " << box.ymin() << " 0" << std::endl;
  }
  void dump_box_to_polylines (const Bbox_3& box, std::ostream& os) const
  {
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

  friend std::ostream& operator<< (std::ostream& os, const Self& orthtree)
  {
    // Create a range of nodes
    auto nodes = orthtree.traverse(Orthtrees::Preorder_traversal());
    // Iterate over the range
    for (auto &n : nodes) {
      // Show the depth
      for (int i = 0; i < n.depth(); ++i)
        os << ". ";
      // Print the node
      os << n << std::endl;
    }
    return os;
  }

  /// \endcond

}; // end class Orthtree

} // namespace CGAL

#include <CGAL/Orthtree/Node.h>

#endif // CGAL_ORTHTREE_H
