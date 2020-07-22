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

namespace CGAL {
namespace Octree {
namespace Node {

/*!
 * \ingroup PkgOctreeClasses
 *
 * \brief A class representing a single node of the tree. Alternatively referred to as a cell, octant, or subtree
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

  /*!
   * \brief Array for containing the children of this node
   */
  typedef std::array<Node<Point_index>, 8> Child_list;

  /*!
   * \brief Set of bits representing this node's relationship to its parent
   *
   * Equivalent to an array of three booleans,
   * where index[0] is whether x is greater,
   * index[1] is whether y is greater,
   * and index[2] is whether z is greater.
   * Used to represent a node's relationship to the center of its parent.
   */
  typedef std::bitset<3> Index;

  /*!
   * \brief Coordinate location representing this node's relationship with the rest of the tree
   *
   * Each value (x, y, z) of a location is calculated by doubling the parent's location
   * and adding the Index.
   * \todo Maybe I should add an example?
   */
  typedef std::array<uint32_t, 3> Int_location;

  /*!
   * \brief A collection of point indices represented by begin and end iterators
   */
  typedef boost::iterator_range<Point_index> Point_range;

  /// @}

private:

  Point_range m_points;

  const Node<Point_index> *m_parent;

  uint8_t m_depth;

  Int_location m_location;

  std::unique_ptr<Child_list> m_children;

public:

  /// \name Construction
  /// @{

  /*!
   * \brief Creates a new node, optionally as the child of a parent
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
  explicit Node(Node<Point_index> *parent = nullptr, Index index = 0) : m_parent(parent), m_depth(0), m_location({0, 0, 0}) {

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
   * \brief Split a node into subnodes
   *
   * Only leaf nodes should be split.
   * When a node is split it is no longer a leaf node.
   * 8 Children are constructed automatically, and their values are set.
   * Contents of this node are _not_ propagated automatically.
   * It's the responsibility of the caller to redistribute the points contained by a node after splitting
   */
  void split() {

    assert(is_leaf());

    m_children = std::make_unique<Child_list>();
    for (int index = 0; index < 8; index++) {

      (*m_children)[index] = std::move(Node<Point_index>(this, {Index(index)}));
    }
  }

  /*!
   * \brief Eliminate this node's children, making it a leaf node
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
   * \brief Access the child nodes of this node by their indices
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
  Node<Point_index> &operator[](int index) {

    assert(!is_leaf());
    assert(0 <= index && index < 8);

    return (*m_children)[index];
  }

  /*!
   * \brief Read-only access the child nodes of this node by their indices
   *
   * \param index The index of the child node, as an int
   * \return A const reference to the node
   */
  const Node<Point_index> &operator[](int index) const {

    assert(!is_leaf());
    assert(0 <= index && index < 8);

    return (*m_children)[index];
  }

  /// @}

  /// \name Property Accessors
  /// @{

  /*!
   * \brief Read-only access to this node's parent
   *
   * \return
   */
  const Node<Point_index> *parent() const { return m_parent; }

  /*!
   * \brief Retrieve this node's depth in the tree
   * \return
   */
  const uint8_t &depth() const { return m_depth; }

  /*!
   * \brief Retrieve this node's location in the tree
   * \return
   */
  Int_location location() const {

    return m_location;
  }

  /*!
   * \brief Retrieve this node's index in relation to its parent
   * \return
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
   * \brief Determine whether this node is a leaf node
   * \return
   */
  bool is_leaf() const { return (!m_children); }

  /*!
   * \brief Determine whether this node is the root node
   * \return
   */
  bool is_root() const { return (!m_parent); }

  /// @}

  /// \name Value Accessors
  /// @{

  /*!
   * \brief Access to the content held by this node
   * \return
   */
  Point_range &points() { return m_points; }

  /*!
   * \brief Read-only access to the content held by this node
   * \return
   */
  const Point_range &points() const { return m_points; }

  /*!
   * \brief Check whether this node contains any points
   * \return
   */
  bool is_empty() const {
    return m_points.empty();
  }

  /*!
   * \brief Count the points contained by this node
   * \return
   */
  std::size_t number_of_points() const {
    return std::distance(m_points.begin(), m_points.end());
  }

  /// @}

  /// \name Operators
  /// @{

  /*!
   * \brief Compare the topology of this node to another node
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
        if (!((*m_children)[i] == rhs[i]))
          return false;
      }
    }

    // If both nodes are leaf nodes, they must be in the same location
    return (location() == rhs.location());
  }

  bool operator!=(const Node &rhs) const {
    return !operator==(rhs);
  }

  /// @}
};
}

}
}

#endif //CGAL_OCTREE_NODE_H
