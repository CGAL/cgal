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
#include <bitset>

namespace CGAL {

  template<class Kernel,
          class PointRange>
  class Octree_node {
  public: // types :
    typedef Octree_node<Kernel, PointRange> Node;
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename PointRange::const_iterator InputIterator;
    typedef typename std::list<InputIterator> IterList;

    typedef std::array<uint32_t, 3> IntPoint;
    typedef std::array<Node, 8> ChildList;

  private: // members :

    //Node *m_children; /* pointer the the 8 possible child nodes. Leaf if NULL */
    std::unique_ptr<ChildList> m_children;

    Node *m_parent;    /* pointer the the single parent node. Root if NULL */
    IntPoint m_location;    /* integer location of current node (x,y,z) on the current depth grid */
    uint8_t m_depth;    /* current depth inside the octree */
    IterList m_points;   /* list of iterators of the input pwn contained in the node */

  public: // functions :
    Octree_node() :
            //m_children(NULL),
            m_parent(NULL),
            m_location(IntPoint{0, 0, 0}),
            m_depth(0) {}

    ~Octree_node() {
      unsplit();
    }

    void unsplit() {
      if (m_children != NULL) {
        for (int i = 0; i < 8; i++)
          (*m_children)[i].unsplit();

        //delete m_children;
        m_children.reset();
      }
    }

    void split() {
      m_children = std::make_unique<ChildList>();
      for (int child_id = 0; child_id < 8; child_id++) {
        (*m_children)[child_id].set_parent(this);
        (*m_children)[child_id].depth() = m_depth + 1;

        for (int j = 0; j < 3; j++) {
          (*m_children)[child_id].location()[j] = 2 * m_location[j] + ((child_id >> j) & 1);
        }
      }
    }

    bool is_leaf() const { return (m_children == NULL); }

    Node *child(const unsigned int index) const {
      if (m_children == NULL || index > 7)
        return NULL;
      else
        return &((*m_children)[index]);
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

  }; // end class Octree_node
}

#endif //OCTREE_OCTREE_NODE_H
