// Copyright (c) 2007-2020  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro, Cédric Portaneri, Tong Zhao

#ifndef CGAL_ORTHTREE_NODE_H
#define CGAL_ORTHTREE_NODE_H

#include <CGAL/license/Orthtree.h>


#include <boost/range/iterator_range.hpp>

#include <array>
#include <memory>
#include <bitset>
#include <cassert>
#include <iostream>

namespace CGAL {

/*!

  \brief represents a single node of the tree. Alternatively referred
  to as a cell, orthant, or subtree.

  A `Node` is a lightweight object and thus generally passed by
  copy. It is also a `Range` whose value type is `Traits::Point_d`.

  \cgalModels Range
 */
template<typename Traits, typename PointRange, typename PointMap>
class Orthtree<Traits, PointRange, PointMap>::Node
{

public:

  /// \name Types
  /// @{

  typedef Orthtree<Traits, PointRange, PointMap> Enclosing; ///< Orthtree type (enclosing class).
  typedef typename Enclosing::Dimension Dimension; ///< Dimension type.
  typedef typename Enclosing::Degree Degree; ///< Degree type.

  /*!
    \brief self typedef for convenience
   */
  typedef Orthtree<Traits, PointRange, PointMap>::Node Self;


  /// \cond SKIP_IN_MANUAL
  /*!
   * \brief array for containing the child nodes of this node
   */
  typedef std::array<Self, Degree::value> Children;
  /// \endcond

  /*!
    \brief set of bits representing this node's relationship to its parent.

    Equivalent to an array of Booleans, where index[0] is whether x
    is greater, index[1] is whether y is greater, index[2] is whether
    z is greater, and so on for higher dimensions if needed.
    Used to represent a node's relationship to the center of its parent.
   */
  typedef std::bitset<Dimension::value> Local_coordinates;

  /*!
    \brief coordinates representing this node's relationship
    with the rest of the tree.

    Each value (x, y, z, ...) of global coordinates is calculated by doubling
    the parent's global coordinates and adding the local coordinates.
   */
  typedef std::array<std::uint32_t, Dimension::value> Global_coordinates;


  typedef typename PointRange::iterator iterator; ///< iterator type.
  typedef typename PointRange::const_iterator const_iterator; ///< constant iterator type.

  /// \cond SKIP_IN_MANUAL
  /*!
    \brief point range type.
   */
  typedef boost::iterator_range<iterator> Point_range;
  /// \endcond

  /*!
    \brief easy access to adjacency directions.
   */
  typedef typename Traits::Adjacency Adjacency;

  /// @}

private:

  // make Node trivially copiabled
  struct Data
  {
    Point_range points;
    Self parent;
    std::uint8_t depth;
    Global_coordinates global_coordinates;
    std::unique_ptr<Children> children;

    Data (Self parent)
      : parent (parent), depth (0) { }
  };

  std::shared_ptr<Data> m_data;

public:

  /// \cond SKIP_IN_MANUAL

  /// \name Construction
  /// @{

  // Default creates null node
  Node() { }

  /*!
    \brief creates a new node, optionally as the child of a parent

    If no parent is provided, the node created is assumed to be the root of a tree.
    This means that the parent reference is a nullptr, and the depth is zero.
    If a parent is provided, the node becomes the child of that parent.
    In that case, an index should be passed, telling this node its relationship to its parent.
    Depth and global coordinates are automatically determined in the constructor,
    and should generally be considered immutable after construction.

    \param parent A reference to the node containing this one
    \param index This node's relationship to its parent
  */
  explicit Node(Self parent, Local_coordinates local_coordinates)
    : m_data (new Data(parent)) {

    if (!parent.is_null()) {

      m_data->depth = parent.m_data->depth + 1;

      for (int i = 0; i < Dimension::value; i++)
        m_data->global_coordinates[i] = (2 * parent.m_data->global_coordinates[i]) + local_coordinates[i];

    }
    else
      for (int i = 0; i < Dimension::value; i++)
        m_data->global_coordinates[i] = 0;
  }

  /// @}

  /// \name Mutators
  /// @{

  /*!
    \brief split a node into subnodes

    Only leaf nodes should be split.
    When a node is split it is no longer a leaf node.
    8 Children are constructed automatically, and their values are set.
    Contents of this node are _not_ propagated automatically.
    It's the responsibility of the caller to redistribute the points contained by a node after splitting
   */
  void split() {

    assert(is_leaf());

    m_data->children = std::make_unique<Children>();
    for (int index = 0; index < Degree::value; index++) {

      (*m_data->children)[index] = std::move(Self(*this, {Local_coordinates(index)}));

    }
  }

  /*!
   * \brief eliminate this node's children, making it a leaf node
   *
   * When a node is un-split, its children are automatically deleted.
   * After un-splitting a node it will be considered a leaf node
   */
  void unsplit() {

    m_data->children.reset();
  }

  /// @}

  /// \endcond

  /// \name Type & Location
  /// @{

  /*!
    \brief returns `true` if the node is null, `false` otherwise.
  */
  bool is_null() const { return (m_data == nullptr); }

  /*!
    \brief returns `true` if the node has no parent, `false` otherwise.
    \pre `!is_null()`
  */
  bool is_root() const
  {
    CGAL_precondition(!is_null());
    return m_data->parent.is_null();
  }

  /*!
    \brief returns `true` if the node has no children, `false` otherwise.
    \pre `!is_null()`
  */
  bool is_leaf() const
  {
    CGAL_precondition(!is_null());
    return (!m_data->children);
  }

  /*!
    \brief returns this node's depth.
    \pre `!is_null()`
  */
  std::uint8_t depth() const
  {
    CGAL_precondition (!is_null());
    return m_data->depth;
  }

  /*!
    \brief returns this node's local coordinates (in relation to its parent).
    \pre `!is_null()`
  */
  Local_coordinates local_coordinates() const {

    CGAL_precondition (!is_null());
    // TODO: There must be a better way of doing this!

    Local_coordinates result;

    for (std::size_t i = 0; i < Dimension::value; ++ i)
      result[i] = global_coordinates()[i] & 1;

    return result;
  }

  /*!
    \brief returns this node's global coordinates.
    \pre `!is_null()`
  */
  Global_coordinates global_coordinates() const
  {
    CGAL_precondition (!is_null());
    return m_data->global_coordinates;
  }


  /// \name Adjacency
  /// @{

  /*!
    \brief returns this node's parent.
    \pre `!is_null()`
  */
  Self parent() const
  {
    CGAL_precondition (!is_null());
    return m_data->parent;
  }

  /*!
    \brief returns the nth child fo this node.

    \pre `!is_null()`
    \pre `!is_leaf()`
    \pre `0 <= index && index < Degree::value`

    The operator can be chained. For example, `n[5][2][3]` returns the
    third child of the second child of the fifth child of a node `n`.
  */
  Self operator[](std::size_t index) const {

    CGAL_precondition (!is_null());
    CGAL_precondition (!is_leaf());
    CGAL_precondition (0 <= index && index < Degree::value);

    return (*m_data->children)[index];
  }

  /*!
    \brief find the directly adjacent node in a specific direction

    \pre `!is_null()`
    \pre `direction.to_ulong < 2 * Dimension::value`

    Adjacent nodes are found according to several properties:
    - adjacent nodes may be larger than the seek node, but never smaller
    - a node has at most `2*Dimension::value` different adjacent nodes (in 3D: left, right, up, down, front, back)
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

    \return the adjacent node if it exists, a null node otherwise.
  */
  Self adjacent_node (Local_coordinates direction) const
  {
    CGAL_precondition(!is_null());

    // Direction:   LEFT  RIGHT  DOWN    UP  BACK FRONT
    // direction:    000    001   010   011   100   101

    // Nodes only have up to 2*dim different adjacent nodes (since cubes have 6 sides)
    CGAL_precondition(direction.to_ulong() < Dimension::value * 2);

    // The root node has no adjacent nodes!
    if (is_root())
      return Self();

    // The least significant bit indicates the sign (which side of the node)
    bool sign = direction[0];

    // The first two bits indicate the dimension/axis (x, y, z)
    uint8_t dimension = (direction >> 1).to_ulong();

    // Create an offset so that the bit-significance lines up with the dimension (e.g. 1, 2, 4 --> 001, 010, 100)
    int8_t offset = (uint8_t) 1 << dimension;

    // Finally, apply the sign to the offset
    offset = (sign ? offset : -offset);

    // Check if this child has the opposite sign along the direction's axis
    if (local_coordinates()[dimension] != sign) {

      // This means the adjacent node is a direct sibling, the offset can be applied easily!
      return parent()[local_coordinates().to_ulong() + offset];
    }

    // Find the parent's neighbor in that direction if it exists
    Self adjacent_node_of_parent = parent().adjacent_node(direction);

    // If the parent has no neighbor, then this node doesn't have one
    if (adjacent_node_of_parent.is_null())
      return Node();

    // If the parent's adjacent node has no children, then it's this node's adjacent node
    if (adjacent_node_of_parent.is_leaf())
      return adjacent_node_of_parent;

    // Return the nearest node of the parent by subtracting the offset instead of adding
    return adjacent_node_of_parent[local_coordinates().to_ulong() - offset];

  }

  /*!
    \brief equivalent to adjacent_node, with an adjacency direction
    rather than a bitset.
   */
  Self adjacent_node(Adjacency adjacency) const {
    return adjacent_node(std::bitset<Dimension::value>(static_cast<int>(adjacency)));
  }

  /// @}

  /// \name Point Range
  /// @{

  /// \cond SKIP_IN_MANUAL
  /*!
   * \brief access to the content held by this node
   * \return a reference to the collection of point indices
   */
  Point_range &points() { return m_data->points; }
  /*!
   * \brief read-only access to the content held by this node
   * \return a read-only reference to the collection of point indices
   */
  const Point_range &points() const { return m_data->points; }

  /// \endcond

  /*!
    \brief returns the number of points of this node.
   */
  bool empty() const {
    return m_data->points.empty();
  }

  /*!
    \brief returns the number of points of this node.
   */
  std::size_t size() const {
    return std::distance(m_data->points.begin(), m_data->points.end());
  }

  /*!
    \brief returns the iterator at the start of the collection of
    points held by this node.
   */
  iterator begin() { return m_data->points.begin(); }

  /*!
    \brief returns the iterator at the end of the collection of
    points held by this node.
   */
  iterator end() { return m_data->points.end(); }

  /*!
    \brief returns the iterator at the start of the collection of
    points held by this node.
   */
  const_iterator begin() const { return m_data->points.begin(); }

  /*!
    \brief returns the iterator at the end of the collection of
    points held by this node.
   */
  const_iterator end() const { return m_data->points.end(); }

  /// @}

  /// \cond SKIP_IN_MANUAL
  /// \name Operators
  /// @{

  /*!
   * \brief compare the topology of this node to another node
   *
   * \todo
   *
   * \param rhs node to compare with
   * \return whether the nodes have different topology
   */
  bool operator==(const Self &rhs) const {
    return m_data == rhs.m_data;
  }

  static bool is_topology_equal (const Self& a, const Self& b)
  {
    CGAL_assertion (!a.is_null() && !b.is_null());

    // If one node is a leaf, and the other isn't, they're not the same
    if (a.is_leaf() != b.is_leaf())
      return false;

    // If both nodes are non-leaf nodes
    if (!a.is_leaf()) {

      // Check all the children
      for (int i = 0; i < Degree::value; ++i) {
        // If any child cell is different, they're not the same
        if (!is_topology_equal(a[i], b[i]))
          return false;
      }
    }

    // If both nodes are leaf nodes, they must be in the same location
    return (a.global_coordinates() == b.global_coordinates());
  }

  friend std::ostream& operator<< (std::ostream& os, const Self& node)
  {
    return internal::print_orthtree_node(os, node);
  }
  /// \endcond
};

}

#endif //CGAL_ORTHTREE_NODE_H
