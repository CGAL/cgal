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

  namespace Octree {

    /*!
     * \ingroup PkgOctreeClasses
     *
     * \brief A class representing a single node of the tree. Alternatively referred to as a cell, octant, or subtree
     *
     * \details A longer description ...
     * Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.
     *
     * \tparam Kernel
     * \tparam PointRange
     */
    template<class Kernel,
            class PointRange>
    class Octree_node {

    public: // types :

      /// \name Types
      /// @{

      /*!
       * \brief
       * \todo
       */
      typedef std::array<uint32_t, 3> IntPoint;

      /*!
       * \brief
       * \todo
       */
      typedef std::array<Node, 8> ChildList;

      /// @}

    private: // Private types

      typedef Octree_node<Kernel, PointRange> Node;
      typedef typename Kernel::FT FT;
      typedef typename Kernel::Point_3 Point;
      typedef typename PointRange::iterator Range_iterator;
      typedef typename std::iterator_traits<Range_iterator>::value_type Range_type;

    private: // members :

      //Node *m_children; /* pointer the the 8 possible child nodes. Leaf if NULL */
      std::unique_ptr <ChildList> m_children; /* pointer the the 8 possible child nodes. Leaf if NULL */

      Node *m_parent;    /* pointer the the single parent node. Root if NULL */
      IntPoint m_location;    /* integer location of current node (x,y,z) on the current depth grid */
      uint8_t m_depth;    /* current depth inside the octree */

      Range_iterator m_points_begin, m_points_end;

    public: // functions :

      /// \name Construction, Destruction
      /// @{

      /*!
       * \brief
       *
       * \todo
       */
      Octree_node() :
              m_children(),
              m_parent(NULL),
              m_location(IntPoint{0, 0, 0}),
              m_depth(0) {}

      /*!
       * \brief
       *
       * \todo
       */
      ~Octree_node() {
        unsplit();
      }

      /// @}

      /// \name Mutators
      /// @{

      /*!
       * \brief
       *
       * \todo
       */
      void unsplit() {
        if (m_children != NULL) {
          for (int i = 0; i < 8; i++)
            (*m_children)[i].unsplit();

          m_children.reset();
        }
      }

      /*!
       * \brief
       *
       * \todo
       */
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

      /// @}

      /// \name Accessors
      /// @{

      bool is_leaf() const { return (m_children == NULL); }

      Node &operator[](int index) {
        return (*m_children)[index];
      }

      const Node &operator[](int index) const {
        return (*m_children)[index];
      }

      Range_iterator &begin() { return m_points_begin; }

      const Range_iterator &begin() const { return m_points_begin; }

      Range_iterator &end() { return m_points_end; }

      const Range_iterator &end() const { return m_points_end; }

      Node *parent() { return m_parent; }

      const Node *parent() const { return m_parent; }

      void set_parent(Node *parent) { m_parent = parent; }

      size_t num_points() const {
        return std::distance(m_points_begin, m_points_end);
      }

      bool is_empty() const { return (num_points() == 0); }

      IntPoint &location() { return m_location; }

      const IntPoint &location() const { return m_location; }

      uint8_t &depth() { return m_depth; }

      const uint8_t &depth() const { return m_depth; }

      std::bitset<3> index() const {

        std::bitset<3> result{};

        result[0] = location()[0] & 1;
        result[1] = location()[1] & 1;
        result[2] = location()[2] & 1;

        return result;
      }

      bool is_sibling(Node *neighbor) const {
        return (m_parent == neighbor->parent());
      }

      /// @}

      /// \name Operators
      /// @{

      /*!
       * \brief
       *
       * \todo
       *
       * \param rhs
       * \return
       */
      bool operator==(Node &rhs) const {

        // Compare the points they contain
        // TODO

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

        // If both nodes are leaf nodes they must be identical (no check necessary)

        // If all other checks pass, the two trees are identical
        return true;
      }

      /// @}

    }; // end class Octree_node

  }
}

#endif //OCTREE_OCTREE_NODE_H
