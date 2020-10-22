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

#include <CGAL/Orthtree/IO.h>

#include <boost/range/iterator_range.hpp>

#include <array>
#include <memory>
#include <bitset>
#include <cassert>
#include <iostream>

namespace CGAL {

/*!
 * \ingroup PkgOrthtreeClasses
 *
 * \brief represents a single node of the tree. Alternatively referred to as a cell, orthant, or subtree
 *
 * \details The role of the node isn't fully stable yet
 *
 * \tparam Point_index is the datatype the node will contain
 */
template<class Traits, class PointRange, class PointMap>
class Orthtree<Traits, PointRange, PointMap>::Node {

public:

  /// \cond SKIP_IN_MANUAL
  typedef Orthtree<Traits, PointRange, PointMap> Parent;
  typedef typename Parent::Dimension Dimension;
  typedef typename Parent::Degree Degree;
  /// \endcond

  /// \name Types
  /// @{

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
  typedef std::array<uint32_t, Dimension::value> Int_location;

  /*!
   * \brief a collection of point indices represented by begin and end iterators
   */
  typedef boost::iterator_range<typename PointRange::iterator> Point_range;

  // TODO: Should I use enum classes?

  /*!
   * \brief the index of a node relative to its parent (a position defined by the corners of a cube)
   *
   * Corners are mapped to numbers as 3-bit integers, in "zyx" order.
   *
   * For example:
   * > right-top-back --> x=1, y=1, z=0 --> zyx = 011 --> 3
   *
   * The following diagram may be a useful reference:
   *
   *           6          7
   *            +--------+
   *           /|       /|             y+
   *          / |      / |             *  z+
   *       2 +--------+ 3|             | *
   *         |  |     |  |             |/
   *         |  +-----|--+             +-----* x+
   *         | / 4    | / 5
   *         |/       |/
   *         +--------+
   *       0           1
   *
   * This lookup table may also be helpful:
   *
   * | Child                 | bitset | number | Enum                  |
   * | --------------------- | ------ | ------ | --------------------- |
   * | left, bottom, back    | 000    | 0      | LEFT_BOTTOM_BACK      |
   * | right, bottom, back   | 001    | 1      | RIGHT_BOTTOM_BACK     |
   * | left, top, back       | 010    | 2      | LEFT_TOP_BACK         |
   * | right, top, back      | 011    | 3      | RIGHT_TOP_BACK        |
   * | left, bottom, front   | 100    | 4      | LEFT_BOTTOM_FRONT     |
   * | right, bottom, front  | 101    | 5      | RIGHT_BOTTOM_FRONT    |
   * | left, top, front      | 110    | 6      | LEFT_TOP_FRONT        |
   * | right, top, front     | 111    | 7      | RIGHT_TOP_FRONT       |
   */
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

  /*!
   * \brief two directions along each axis in cartesian space, relative to a node
   *
   * Directions are mapped to numbers as 3-bit integers,
   * though the numbers 6 and 7 are not used because there are only 6 different directions.
   *
   * The first two bits indicate the axis (00 = x, 01 = y, 10 = z),
   * the third bit indicates the direction along that axis (0 = -, 1 = +).
   *
   * The following diagram may be a useful reference:
   *
   *            3 *
   *              |  * 5
   *              | /                  y+
   *              |/                   *  z+
   *     0 *------+------* 1           | *
   *             /|                    |/
   *            / |                    +-----* x+
   *         4 *  |
   *              * 2
   *
   * This lookup table may also be helpful:
   *
   * | Direction | bitset | number | Enum  |
   * | --------- | ------ | ------ | ----- |
   * | `-x`      | 000    | 0      | LEFT  |
   * | `+x`      | 001    | 1      | RIGHT |
   * | `-y`      | 010    | 2      | DOWN  |
   * | `+y`      | 011    | 3      | UP    |
   * | `-z`      | 100    | 4      | BACK  |
   * | `+z`      | 101    | 5      | FRONT |
   */
  enum Direction {
    LEFT,
    RIGHT,
    DOWN,
    UP,
    BACK,
    FRONT
  };

  /// @}

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
  explicit Node(Self *parent = nullptr, Index index = 0) : m_parent(parent), m_depth(0) {

    if (parent) {

      m_depth = parent->m_depth + 1;

      for (int i = 0; i < Dimension::value; i++)
        m_location[i] = (2 * parent->m_location[i]) + index[i];

    }
    else
      for (int i = 0; i < Dimension::value; i++)
        m_location[i] = 0;
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
    for (int index = 0; index < Degree::value; index++) {

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
    assert(0 <= index && index < Degree::value);

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
    assert(0 <= index && index < Degree::value);

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
   * \return the index of this node
   */
  Index index() const {

    // TODO: There must be a better way of doing this!

    Index result;

    for (std::size_t i = 0; i < Dimension::value; ++ i)
      result[i] = location()[i] & 1;

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

  /*!
   * \brief find the directly adjacent node in a specific direction
   *
   * Adjacent nodes are found according to several properties:
   * - Adjacent nodes may be larger than the seek node, but never smaller
   * - A node can have no more than 6 different adjacent nodes (left, right, up, down, front, back)
   * - A node is free to have fewer than 6 adjacent nodes
   *   (e.g. edge nodes have no neighbors in some directions, the root node has none at all).
   * - Adjacent nodes are not required to be leaf nodes
   *
   *
   * Here's a diagram demonstrating the concept for a quadtree.
   * Because it's in 2d space, the seek node has only four neighbors (up, down, left, right)
   *
   *     +---------------+---------------+
   *     |               |               |
   *     |               |               |
   *     |               |               |
   *     |       A       |               |
   *     |               |               |
   *     |               |               |
   *     |               |               |
   *     +-------+-------+---+---+-------+
   *     |       |       |   |   |       |
   *     |   A   |  (S)  +---A---+       |
   *     |       |       |   |   |       |
   *     +---+---+-------+---+---+-------+
   *     |   |   |       |       |       |
   *     +---+---+   A   |       |       |
   *     |   |   |       |       |       |
   *     +---+---+-------+-------+-------+
   *
   *         (S) : Seek node
   *          A  : Adjacent node
   *
   * Note how the top adjacent node is larger than the seek node.
   * The right adjacent node is the same size, even though it contains further subdivisions.
   *
   * This implementation returns a pointer to the adjacent node if it's found.
   * If there is no adjacent node in that direction, it returns nullptr.
   *
   * \todo explain how direction is encoded
   *
   * \param direction which way to find the adjacent node relative to this one
   * \return a pointer to the adjacent node if it exists
   */
  const Self *adjacent_node(std::bitset<Dimension::value> direction) const {

    // Direction:   LEFT  RIGHT  DOWN    UP  BACK FRONT
    // direction:    000    001   010   011   100   101

    // Nodes only have up to 6 different adjacent nodes (since cubes have 6 sides)
    assert(direction.to_ulong() < 6);

    // The root node has no adjacent nodes!
    if (is_root())
      return nullptr;

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
    auto *adjacent_node_of_parent = parent()->adjacent_node(direction);

    // If the parent has no neighbor, then this node doesn't have one
    if (!adjacent_node_of_parent)
      return nullptr;

    // If the parent's adjacent node has no children, then it's this node's adjacent node
    if (adjacent_node_of_parent->is_leaf())
      return adjacent_node_of_parent;

    // Return the nearest node of the parent by subtracting the offset instead of adding
    return &(*adjacent_node_of_parent)[index().to_ulong() - offset];

  }

  /*!
   * \brief equivalent to adjacent_node, with a Direction rather than a bitset
   */
  const Self *adjacent_node(Direction direction) const {
    return adjacent_node(std::bitset<Dimension::value>(static_cast<int>(direction)));
  }

  /*!
   * \brief equivalent to adjacent_node, except non-const
   */
  Self *adjacent_node(std::bitset<Dimension::value> direction) {
    return const_cast<Self *>(const_cast<const Self *>(this)->adjacent_node(direction));
  }

  /*!
   * \brief equivalent to adjacent_node, with a Direction rather than a bitset and non-const
   */
  Self *adjacent_node(Direction direction) {
    return adjacent_node(std::bitset<Dimension::value>(static_cast<int>(direction)));
  }

  /// @}

  /// \name Value Accessors
  /// @{

  /*!
   * \brief access to the content held by this node
   * \return a reference to the collection of point indices
   */
  Point_range &points() { return m_points; }

  /*!
   * \brief access to the points via iterator
   * \return the iterator at the start of the collection of point indices held by this node
   */
  typename Point_range::iterator begin() { return m_points.begin(); }

  /*!
    * \brief access to the points via iterator
    * \return the iterator at the end of the collection of point indices held by this node
    */
  typename Point_range::iterator end() { return m_points.end(); }

  /*!
   * \brief read-only access to the content held by this node
   * \return a read-only reference to the collection of point indices
   */
  const Point_range &points() const { return m_points; }

  /*!
    * \brief read-only access to the points via iterator
    * \return the iterator at the start of the collection of point indices held by this node
    */
  typename Point_range::iterator begin() const { return m_points.begin(); }

  /*!
    * \brief read-only access to the points via iterator
    * \return the iterator at the end of the collection of point indices held by this node
    */
  typename Point_range::iterator end() const { return m_points.end(); }

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
   * \param rhs node to compare with
   * \return whether the nodes have different topology
   */
  bool operator==(const Self &rhs) const {

    // TODO: Should I compare the values they contain
//          if (m_points != rhs.m_points)
//            return false;

    // If one node is a leaf, and the other isn't, they're not the same
    if (is_leaf() != rhs.is_leaf())
      return false;

    // If both nodes are non-leaf nodes
    if (!is_leaf()) {

      // Check all the children
      for (int i = 0; i < Degree::value; ++i) {

        // If any child cell is different, they're not the same
        if ((*m_children)[i] != rhs[i])
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
