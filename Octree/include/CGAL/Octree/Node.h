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

#ifndef CGAL_OCTREE_NODE_H
#define CGAL_OCTREE_NODE_H

#include <CGAL/license/Octree.h>

#include <boost/range/iterator_range.hpp>

#include <array>
#include <memory>
#include <bitset>
#include <cassert>
#include <iostream>

namespace CGAL {
namespace Octree {

/*!
 * \ingroup PkgOctreeClasses
 *
 * \brief represents a single node of the tree. Alternatively referred to as a cell, octant, or subtree
 *
 * \details The role of the node isn't fully stable yet
 *
 * \tparam Point_index is the datatype the node will contain
 */
template<typename Point_index>
class Node {

public:

  /// \name Types
  /// @{

  typedef Node<Point_index> Self;

  /*!
   * \brief array for containing the child nodes of this node
   */
  typedef std::array<Self, 8> Children;

  /*!
   * \brief set of bits representing this node's relationship to its parent
   *
   * Equivalent to an array of three booleans,
   * where index[0] is whether x is greater,
   * index[1] is whether y is greater,
   * and index[2] is whether z is greater.
   * Used to represent a node's relationship to the center of its parent.
   */
  typedef std::bitset<3> Index;

  /*!
   * \brief coordinate location representing this node's relationship with the rest of the tree
   *
   * Each value (x, y, z) of a location is calculated by doubling the parent's location
   * and adding the Index.
   * \todo Maybe I should add an example?
   */
  typedef std::array<uint32_t, 3> Int_location;

  /*!
   * \brief a collection of point indices represented by begin and end iterators
   */
  typedef boost::iterator_range<Point_index> Point_range;

  /// @}

  // TODO: There's probably a better name for this
  // TODO: Should I use an enum class?
  enum Child {
    LEFT_BOTTOM_BACK,
    RIGHT_BOTTOM_BACK,
    LEFT_TOP_BACK,
    RIGHT_TOP_BACK,
    LEFT_BOTTOM_FRONT,
    RIGHT_BOTTOM_FRONT,
    LEFT_TOP_FRONT,
    RIGHT_TOP_FRONT
  };

  enum Direction {
    LEFT,
    RIGHT,
    DOWN,
    UP,
    BACK,
    FRONT
  };


private:

  Point_range m_points;

  const Self *m_parent;

  uint8_t m_depth;

  Int_location m_location;

  std::unique_ptr<Children> m_children;

public:

  /// \name Construction
  /// @{

  /*!
   * \brief creates a new node, optionally as the child of a parent
   *
   * If no parent is provided, the node created is assumed to be the root of a tree.
   * This means that the parent reference is a nullptr, and the depth is zero.
   * If a parent is provided, the node becomes the child of that parent.
   * In that case, an index should be passed, telling this node its relationship to its parent.
   * Depth and location are automatically determined in the constructor,
   * and should generally be considered immutable after construction.
   *
   * \param parent A reference to the node containing this one
   * \param index This node's relationship to its parent
   */
  explicit Node(Self *parent = nullptr, Index index = 0) : m_parent(parent), m_depth(0),
                                                           m_location({0, 0, 0}) {

    if (parent) {

      m_depth = parent->m_depth + 1;

      for (int i = 0; i < 3; i++)
        m_location[i] = (2 * parent->m_location[i]) + index[i];

    }
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

    m_children = std::make_unique<Children>();
    for (int index = 0; index < 8; index++) {

      (*m_children)[index] = std::move(Self(this, {Index(index)}));
    }
  }

  /*!
   * \brief eliminate this node's children, making it a leaf node
   *
   * When a node is un-split, its children are automatically deleted.
   * After un-splitting a node it will be considered a leaf node
   */
  void unsplit() {

    m_children.reset();
  }

  /// @}

  /// \name Child Accessors
  /// @{

  /*!
   * \brief access the child nodes of this node by their indices
   *
   * \todo Explain how index values map to the Index type
   *
   * Retrieves a reference to the child node described by the index.
   * The operator can be chained.
   * for example, to access the third child of the second child of the fifth child of a node `n`
   *
   *     n[5][2][3];
   *
   * \param index The index of the child node, as an int
   * \return A reference to the node
   */
  Self &operator[](int index) {

    assert(!is_leaf());
    assert(0 <= index && index < 8);

    return (*m_children)[index];
  }

  /*!
   * \brief read-only access the child nodes of this node by their indices
   *
   * \param index The index of the child node, as an int
   * \return A const reference to the node
   */
  const Self &operator[](int index) const {

    assert(!is_leaf());
    assert(0 <= index && index < 8);

    return (*m_children)[index];
  }

  /// @}

  /// \name Property Accessors
  /// @{

  /*!
   * \brief read-only access to this node's parent
   *
   * Ownership of a node is not equivalent to ownership of the entire tree,
   * so it's not possible to obtain write access to a node's parent,
   * only its children.
   * Note that the return type is nullable,
   * attempting to find the parent of a root node will return null.
   * \todo Should I instead assert the node isn't root? (that would make this undefined behavior)
   *
   * \return A const pointer to the parent of this node (possibly nullptr)
   */
  const Self *parent() const { return m_parent; }

  /*!
   * \brief retrieve this node's depth in the tree
   * \return the depth of this node, where root has a depth of 0
   */
  const uint8_t &depth() const { return m_depth; }

  /*!
   * \brief retrieve this node's location in the tree
   *
   * \todo Should I link to an explanation of the location type?
   *
   * \return this node's location
   */
  Int_location location() const {

    return m_location;
  }

  /*!
   * \brief retrieve this node's index in relation to its parent
   * \return the index of this nod3
   */
  Index index() const {

    // TODO: There must be a better way of doing this!

    Index result;

    result[0] = location()[0] & 1;
    result[1] = location()[1] & 1;
    result[2] = location()[2] & 1;

    return result;
  }

  /*!
   * \brief determine whether this node is a leaf node
   * \return whether this node has no children
   */
  bool is_leaf() const { return (!m_children); }

  /*!
   * \brief determine whether this node is the root node
   * \return whether this node has no parent
   */
  bool is_root() const { return (!m_parent); }

  /// @}

  /// \name Value Accessors
  /// @{

  /*!
   * \brief access to the content held by this node
   * \return a reference to the collection of point indices
   */
  Point_range &points() { return m_points; }

  /*!
   * \brief read-only access to the content held by this node
   * \return a read-only reference to the collection of point indices
   */
  const Point_range &points() const { return m_points; }

  /*!
   * \brief check whether this node contains any points
   * \return if this node contains no points
   */
  bool empty() const {
    return m_points.empty();
  }

  /*!
   * \brief count the points contained by this node
   * \return the number of points this node owns
   */
  std::size_t size() const {
    return std::distance(m_points.begin(), m_points.end());
  }

  /// @}

  /// \name Operators
  /// @{

  /*!
   * \brief compare the topology of this node to another node
   *
   * \todo
   *
   * \param rhs
   * \return
   */
  bool operator==(const Node &rhs) const {

    // TODO: Should I compare the values they contain
//          if (m_points != rhs.m_points)
//            return false;

    // If one node is a leaf, and the other isn't, they're not the same
    if (is_leaf() != rhs.is_leaf())
      return false;

    // If both nodes are non-leaf nodes
    if (!is_leaf()) {

      // Check all the children
      for (int i = 0; i < 8; ++i) {

        // If any child cell is different, they're not the same
        if ((*m_children)[i] != rhs[i])
          return false;
      }
    }

    // If both nodes are leaf nodes, they must be in the same location
    return (location() == rhs.location());
  }

  bool operator!=(const Node &rhs) const {
    return !operator==(rhs);
  }

  const Self *adjacent(Direction direction) const {
    return adjacent(std::bitset<3>(static_cast<int>(direction)));
  }

  const Self *adjacent(std::bitset<3> direction) const {

    // Nodes only have up to 6 different adjacent nodes (since cubes have 6 sides)
    assert(direction.to_ulong() < 6);

    // The root node has no adjacent nodes!
    if (is_root())
      return nullptr;

    // Direction:   LEFT  RIGHT  DOWN    UP  BACK FRONT
    // direction:    000    001   010   011   100   101

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
      return &(*parent())[index().to_ulong() + offset];
    }

    // Find the parent's neighbor in that direction if it exists
    auto *adjacent_node_of_parent = parent()->adjacent(direction);

    // If the parent has no neighbor, then this node doesn't have one
    if (!adjacent_node_of_parent)
      return nullptr;

    // If the parent's adjacent node has no children, then it's this node's adjacent node
    if (adjacent_node_of_parent->is_leaf())
      return adjacent_node_of_parent;

    // Return the nearest node of the parent by subtracting the offset instead of adding
    return &(*adjacent_node_of_parent)[index().to_ulong() - offset];

  }

  /// @}
};

}
}

#endif //CGAL_OCTREE_NODE_H
