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

#ifndef CGAL_OCTREE_WALKER_CRITERION_H
#define CGAL_OCTREE_WALKER_CRITERION_H

#include <iostream>
#include <boost/range/iterator_range.hpp>
#include "Node.h"
#include "Walker_iterator.h"

namespace CGAL {

namespace Octree {

namespace Walker {

template<class Value>
const Node::Node <Value> *next_sibling(const Node::Node <Value> *n) {

  // Passing null returns the first node
  if (nullptr == n)
    return nullptr;

  // If this node has no parent, it has no siblings
  if (nullptr == n->parent())
    return nullptr;

  // Find out which child this is
  std::size_t index = n->index().to_ulong();

  // Return null if this is the last child
  if (7 == index)
    return nullptr;

  // Otherwise, return the next child
  return &((*n->parent())[index + 1]);
}

template<class Value>
const Node::Node <Value> *next_sibling_up(const Node::Node <Value> *n) {

  if (!n)
    return nullptr;

  auto up = n->parent();

  while (nullptr != up) {

    if (nullptr != next_sibling(up))
      return next_sibling(up);

    up = up->parent();
  }

  return nullptr;
}

template<class Value>
const Node::Node <Value> *deepest_first_child(const Node::Node <Value> *n) {

  if (!n)
    return nullptr;

  // Find the deepest child on the left
  auto first = n;
  while (!first->is_leaf())
    first = &(*first)[0];
  return first;
}


/*!
 * \brief Walker for preorder traversal
 */
struct Preorder {

  /*!
   * \brief Retrieve the first node of a tree in a preorder traversal, given the root
   *
   * \tparam Value
   * \param root
   * \return
   */
  template<class Value>
  const Node::Node <Value> *first(const Node::Node <Value> *root) const {
    return root;
  }

  /*!
   * \brief Retrieve the next node of a tree in a preorder traversal, given the current one
   *
   * \tparam Value
   * \param n
   * \return
   */
  template<class Value>
  const Node::Node <Value> *next(const Node::Node <Value> *n) const {

    if (n->is_leaf()) {

      auto next = next_sibling(n);

      if (nullptr == next) {

        return next_sibling_up(n);
      }

      return next;

    } else {

      // Return the first child of this node
      return &(*n)[0];
    }

  }
};

/*!
 * \brief Walker for leaves-only traversal
 */
struct Postorder {

  /*!
   * \brief Retrieve the first node of a tree in a postorder traversal, given the root
   *
   * \tparam Value
   * \param root
   * \return
   */
  template<class Value>
  const Node::Node <Value> *first(const Node::Node <Value> *root) const {

    return deepest_first_child(root);
  }

  /*!
   * \brief Retrieve the next node of a tree in a postorder traversal, given the current one
   *
   * \tparam Value
   * \param n
   * \return
   */
  template<class Value>
  const Node::Node <Value> *next(const Node::Node <Value> *n) const {

    auto next = deepest_first_child(next_sibling(n));

    if (!next)
      next = n->parent();

    return next;
  }
};

/*!
 * \brief Walker for leaves-only traversal
 */
struct Leaves {

  /*!
   * \brief Retrieve the first node of a tree in a leaves-only traversal, given the root
   *
   * \tparam Value
   * \param root
   * \return
   */
  template<class Value>
  const Node::Node <Value> *first(const Node::Node <Value> *root) const {

    return deepest_first_child(root);
  }

  /*!
   * \brief Retrieve the next node of a tree in a leaves-only traversal, given the current one
   *
   * \tparam Value
   * \param n
   * \return
   */
  template<class Value>
  const Node::Node <Value> *next(const Node::Node <Value> *n) const {

    auto next = deepest_first_child(next_sibling(n));

    if (!next)
      next = deepest_first_child(next_sibling_up(n));

    return next;
  }
};

}

}

}

#endif //CGAL_OCTREE_WALKER_CRITERION_H
