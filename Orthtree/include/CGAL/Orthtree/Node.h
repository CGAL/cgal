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
  \brief represents a single node of the tree. Alternatively referred to as a cell, orthant, or subtree

  \details The role of the node isn't fully stable yet

  \tparam Point_index is the datatype the node will contain
 */
template<typename Traits, typename PointRange, typename PointMap>
class Orthtree<Traits, PointRange, PointMap>::Node
{

public:

  /// \name Types
  /// @{

  typedef Orthtree<Traits, PointRange, PointMap> Parent; ///< Orthtree type.
  typedef typename Parent::Dimension Dimension; ///< Dimension type.
  typedef typename Parent::Degree Degree; ///< Degree type.

  /*!
   * \brief self typedef for convenience
   */
  typedef Orthtree<Traits, PointRange, PointMap>::Node Self;

  /*!
   * \brief array for containing the child nodes of this node
   */
  typedef std::array<Self, Degree::value> Children;

  /*!
   * \brief set of bits representing this node's relationship to its parent
   *
   * Equivalent to an array of booleans, where index[0] is whether x
   * is greater, index[1] is whether y is greater, index[2] is whether
   * z is greater, and so on for higher dimensions if neede.
   * Used to represent a node's relationship to the center of its parent.
   */
  typedef std::bitset<Dimension::value> Index;

  /*!
   * \brief coordinate location representing this node's relationship with the rest of the tree
   *
   * Each value (x, y, z, ...) of a location is calculated by doubling
   * the parent's location and adding the Index.
   * \todo Maybe I should add an example?
   */
  typedef std::array<std::uint32_t, Dimension::value> Int_location;

  /*!
   * \brief a collection of point indices represented by begin and end iterators
   */
  typedef boost::iterator_range<typename PointRange::iterator> Point_range;

  typedef typename Traits::Adjacency Adjacency;

  /// @}

private:

  // make Node trivially copiabled
  struct Data
  {
    Point_range points;
    Self parent;
    std::uint8_t depth;
    Int_location location;
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
    Depth and location are automatically determined in the constructor,
    and should generally be considered immutable after construction.

    \param parent A reference to the node containing this one
    \param index This node's relationship to its parent
  */
  explicit Node(Self parent, Index index)
    : m_data (new Data(parent)) {

    if (!parent.is_null()) {

      m_data->depth = parent.m_data->depth + 1;

      for (int i = 0; i < Dimension::value; i++)
        m_data->location[i] = (2 * parent.m_data->location[i]) + index[i];

    }
    else
      for (int i = 0; i < Dimension::value; i++)
        m_data->location[i] = 0;
  }

  /// @}

  /// \name Mutators
  /// @{

  /*!
   * \brief split a node into subnodes
   *
   * Only leaf nodes should be split.
   * When a node is split it is no longer a leaf node.
   * 8 Children are constructed automatically, and their values are set.
   * Contents of this node are _not_ propagated automatically.
   * It's the responsibility of the caller to redistribute the points contained by a node after splitting
   */
  void split() {

    assert(is_leaf());

    m_data->children = std::make_unique<Children>();
    for (int index = 0; index < Degree::value; index++) {

      (*m_data->children)[index] = std::move(Self(*this, {Index(index)}));

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

  /// \name Accessors
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

  /// @}

  /// \name Property Accessors
  /// @{

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
    \brief returns this node's location.
    \pre `!is_null()`
  */
  Int_location location() const
  {
    CGAL_precondition (!is_null());
    return m_data->location;
  }

  /*!
    \brief returns this node's index in relation to its parent.
    \pre `!is_null()`
  */
  Index index() const {

    CGAL_precondition (!is_null());
    // TODO: There must be a better way of doing this!

    Index result;

    for (std::size_t i = 0; i < Dimension::value; ++ i)
      result[i] = location()[i] & 1;

    return result;
  }

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
    \brief find the directly adjacent node in a specific direction

    \pre `!is_null()`
    \pre `direction.to_ulong < 2 * Dimension::value`

    Adjacent nodes are found according to several properties:
    - Adjacent nodes may be larger than the seek node, but never smaller
    - A node can have no more than 6 different adjacent nodes (left, right, up, down, front, back)
    - A node is free to have fewer than 6 adjacent nodes
    (e.g. edge nodes have no neighbors in some directions, the root node has none at all).
    - Adjacent nodes are not required to be leaf nodes

    Here's a diagram demonstrating the concept for a quadtree.
    Because it's in 2d space, the seek node has only four neighbors (up, down, left, right)

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

    (S) : Seek node
    A  : Adjacent node

    Note how the top adjacent node is larger than the seek node.
    The right adjacent node is the same size, even though it contains further subdivisions.

    This implementation returns a pointer to the adjacent node if it's found.
    If there is no adjacent node in that direction, it returns nullptr.

    \todo explain how direction is encoded

    \param direction which way to find the adjacent node relative to this one
    \return a pointer to the adjacent node if it exists
  */
  Self adjacent_node (std::bitset<Dimension::value> direction) const
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
    if (index()[dimension] != sign) {

      // This means the adjacent node is a direct sibling, the offset can be applied easily!
      return parent()[index().to_ulong() + offset];
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
    return adjacent_node_of_parent[index().to_ulong() - offset];

  }

  /*!
   * \brief equivalent to adjacent_node, with a Direction rather than a bitset
   */
  Self adjacent_node(Adjacency adjacency) const {
    return adjacent_node(std::bitset<Dimension::value>(static_cast<int>(adjacency)));
  }

  /// @}

  /// \name Value Accessors
  /// @{

  /*!
   * \brief access to the content held by this node
   * \return a reference to the collection of point indices
   */
  Point_range &points() { return m_data->points; }

  /*!
   * \brief access to the points via iterator
   * \return the iterator at the start of the collection of point indices held by this node
   */
  typename Point_range::iterator begin() { return m_data->points.begin(); }

  /*!
   * \brief access to the points via iterator
   * \return the iterator at the end of the collection of point indices held by this node
   */
  typename Point_range::iterator end() { return m_data->points.end(); }

  /*!
   * \brief read-only access to the content held by this node
   * \return a read-only reference to the collection of point indices
   */
  const Point_range &points() const { return m_data->points; }

  /*!
   * \brief read-only access to the points via iterator
   * \return the iterator at the start of the collection of point indices held by this node
   */
  typename Point_range::iterator begin() const { return m_data->points.begin(); }

  /*!
   * \brief read-only access to the points via iterator
   * \return the iterator at the end of the collection of point indices held by this node
   */
  typename Point_range::iterator end() const { return m_data->points.end(); }

  /*!
   * \brief check whether this node contains any points
   * \return if this node contains no points
   */
  bool empty() const {
    return m_data->points.empty();
  }

  /*!
   * \brief count the points contained by this node
   * \return the number of points this node owns
   */
  std::size_t size() const {
    return std::distance(m_data->points.begin(), m_data->points.end());
  }

  /// @}

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

    // TODO: Should I compare the values they contain
//          if (m_data->points != rhs.m_data->points)
//            return false;

    if (is_null() || rhs.is_null())
      return (is_null() == rhs.is_null());

    // If one node is a leaf, and the other isn't, they're not the same
    if (is_leaf() != rhs.is_leaf())
      return false;

    // If both nodes are non-leaf nodes
    if (!is_leaf()) {

      // Check all the children
      for (int i = 0; i < Degree::value; ++i) {

        // If any child cell is different, they're not the same
        if ((*m_data->children)[i] != rhs[i])
          return false;
      }
    }

    // If both nodes are leaf nodes, they must be in the same location
    return (location() == rhs.location());
  }

  /*!
   * \brief compares the topology of this node to another node
   *
   * \todo
   *
   * \param rhs node to compare with
   * \return whether the trees have different topology
   */
  bool operator!=(const Self &rhs) const {
    return !operator==(rhs);
  }

  /// @}

  /// \cond SKIP_IN_MANUAL
  friend std::ostream& operator<< (std::ostream& os, const Self& node)
  {
    return internal::print_orthtree_node(os, node);
  }
  /// \endcond
};

}

#endif //CGAL_ORTHTREE_NODE_H
