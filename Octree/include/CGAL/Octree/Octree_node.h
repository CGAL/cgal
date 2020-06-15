// Copyright (c) 2007-2008  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Tong Zhao, Cédric Portaneri

#ifndef OCTREE_OCTREE_NODE_H
#define OCTREE_OCTREE_NODE_H

#include <CGAL/bounding_box.h>
#include <boost/iterator/transform_iterator.hpp>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

/*
 * These headers were not included here originally
 * Adding them was necessary to make this header self sufficient
 */
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <iostream>
#include <fstream>

#include <stack>
#include <queue>
#include <vector>
#include <math.h>

namespace CGAL {

  struct IntPoint_3 {
    uint32_t m_coords[3] = {0, 0, 0};

    IntPoint_3() {}

    IntPoint_3(int x, int y, int z) {
      m_coords[0] = x;
      m_coords[1] = y;
      m_coords[2] = z;
    }

    uint32_t &x() { return m_coords[0]; }

    const uint32_t &x() const { return m_coords[0]; }

    uint32_t &y() { return m_coords[1]; }

    const uint32_t &y() const { return m_coords[1]; }

    uint32_t &z() { return m_coords[2]; }

    const uint32_t &z() const { return m_coords[2]; }

    uint32_t &operator[](int i) {
      assert(!(i < 0 || i > 2));
      return m_coords[i];
    }

    const uint32_t &operator[](int i) const {
      assert(i < 0 || i > 2);
      return m_coords[i];
    }

    bool operator==(const IntPoint_3 &other) const {
      return m_coords[0] == other.x() &&
             m_coords[1] == other.y() &&
             m_coords[2] == other.z();
    }

    bool operator<(const IntPoint_3 &other) const {
      if (x() != other.x()) {
        return (x() < other.x());
      }
      if (y() != other.y()) {
        return (y() < other.y());
      }
      return (z() < other.z());
    }
  };

  template<class Kernel,
          class PointRange>
  class Octree_node {
  public: // types :
    typedef Octree_node<Kernel, PointRange> Node;
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;
    typedef IntPoint_3 IntPoint;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename PointRange::const_iterator InputIterator;
    typedef typename std::list<InputIterator> IterList;

  private: // members :

    // TODO: Replace with pointer to fixed size array, or std::optional<std::array<Node, 8>> ?
    Node *m_children; /* pointer the the 8 possible child nodes. Leaf if NULL */

    Node *m_parent;    /* pointer the the single parent node. Root if NULL */

    IntPoint m_location;    /* integer location of current node (x,y,z) on the current depth grid */

    uint8_t m_depth;    /* current depth inside the octree */
    IterList m_points;   /* list of iterators of the input pwn contained in the node */

  public: // functions :
    Octree_node() :
            m_children(NULL),
            m_parent(NULL),
            m_location(IntPoint(0, 0, 0)),
            m_depth(0) {}

    ~Octree_node() {
      unsplit();
    }

    void unsplit() {
      if (m_children != NULL) {
        for (int i = 0; i < 8; i++)
          m_children[i].unsplit();

        delete[] m_children;
        m_children = NULL;
      }
    }

    void split() {
      m_children = new Node[8];
      for (int child_id = 0; child_id < 8; child_id++) {
        m_children[child_id].set_parent(this);
        m_children[child_id].depth() = m_depth + 1;

        for (int j = 0; j < 3; j++) {
          m_children[child_id].location()[j] = 2 * m_location[j] + ((child_id >> j) & 1);
        }
      }
    }

    bool is_leaf() const { return (m_children == NULL); }

    Node *children() { return m_children; }

    const Node *children() const { return m_children; }

    Node *child(const unsigned int index) const {
      if (m_children == NULL || index > 7)
        return NULL;
      else
        return &(m_children[index]);
    }

    Node *parent() { return m_parent; }

    const Node *parent() const { return m_parent; }

    void set_parent(Node *parent) { m_parent = parent; }

    IterList &points() { return m_points; }

    const IterList &points() const { return m_points; }

    void add_point(InputIterator point) { m_points.push_back(point); }

    size_t num_points() const { return m_points.size(); }

    bool is_empty() const { return (m_points.size() == 0); }

    IntPoint &location() { return m_location; }

    const IntPoint &location() const { return m_location; }

    uint8_t &depth() { return m_depth; }

    const uint8_t &depth() const { return m_depth; }

    bool is_sibling(Node *neighbor) const {
      return (m_parent == neighbor->parent());
    }

    // dir: LEFT = 000, RIGHT = 001, DOWN = 010, UP = 011, BACK = 100, FRONT= 101
    Node *find_greater_or_equal_neighbor(int dir) const {
      if (m_parent == NULL) return NULL;

      unsigned int dir_axis = dir & 1;  // 0, 1, 0, 1, 0, 1
      unsigned int bit_axis = dir >> 1; // 0, 0, 1, 1, 2, 2
      unsigned int dir_idx_offset = 1;  // -1, 1, -2, 2, -4, 4
      dir_idx_offset <<= bit_axis;
      if (!dir_axis) dir_idx_offset = -dir_idx_offset;

      for (int child_id = 0; child_id < 8; child_id++) {
        // is 'this' an opposite 'dir' child?
        if (((child_id >> bit_axis) & 1) != dir_axis && m_parent->child(child_id) == this) {
          // return 'dir' sibling child
          return m_parent->child(child_id + dir_idx_offset);
        }
      }

      Node *parent_neighbor = m_parent->find_greater_or_equal_neighbor(dir);
      if (parent_neighbor == NULL || parent_neighbor->is_leaf()) return parent_neighbor;
      for (int child_id = 0; child_id < 8; child_id++) {
        // 'this' is guaranted to be a 'dir' child
        if (((child_id >> bit_axis) & 1) == dir_axis && m_parent->child(child_id) == this) {
          // return opposite 'dir' neighbor child
          return parent_neighbor->child(child_id - dir_idx_offset);
        }
      }

      return NULL;
    }

    // dir: LEFT = 000, RIGHT = 001, DOWN = 010, UP = 011, BACK = 100, FRONT= 101
    std::list<Node *> find_smaller_neighbors(Node *ge_neighbor, int dir) const {
      std::list<Node *> le_neighbors;
      unsigned int dir_axis = dir & 1;  // 0, 1, 0, 1, 0, 1
      unsigned int bit_axis = dir >> 1; // 0, 0, 1, 1, 2, 2

      std::queue<Node *> possible_neighbors;
      if (ge_neighbor != NULL) possible_neighbors.push(ge_neighbor);
      while (!possible_neighbors.empty()) {
        Node *node = possible_neighbors.front();

        if (node->is_leaf()) {
          le_neighbors.push_back(node);
        } else {
          for (int child_id = 0; child_id < 8; child_id++) {
            if (((child_id >> bit_axis) & 1) != dir_axis) {
              // add to queue the opposite 'dir' neighbor child
              possible_neighbors.push(node->child(child_id));
            }
          }
        }
        possible_neighbors.pop();
      }
      return le_neighbors;
    }

    bool is_balanced() const {
      for (int dir = 0; dir < 6; dir++) {
        Node *ge_neighbor = find_greater_or_equal_neighbor(dir);
        std::list<Node *> neighbors = find_smaller_neighbors(ge_neighbor, dir);
        for (Node *neighbor : neighbors) {
          if (neighbor != NULL && !is_sibling(neighbor)
              && (std::abs(this->depth() - neighbor->depth()) > 1)) {
            return false;
          }
        }
      }
      return true;
    }

    std::list<Node *> find_unbalanced_neighbors_to_split() const {
      std::list<Node *> neighbors_to_split;
      for (int dir = 0; dir < 6; dir++) {
        Node *neighbor = find_greater_or_equal_neighbor(dir);
        if (neighbor != NULL && !is_sibling(neighbor) && neighbor->is_leaf()
            && ((this->depth() - neighbor->depth()) > 1)) {
          neighbors_to_split.push_back(neighbor);
        }
      }
      return neighbors_to_split;
    }

  }; // end class Octree_node
}

#endif //OCTREE_OCTREE_NODE_H
