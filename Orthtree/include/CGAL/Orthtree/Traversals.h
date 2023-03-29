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

#ifndef CGAL_ORTHTREE_TRAVERSALS_H
#define CGAL_ORTHTREE_TRAVERSALS_H

#include <CGAL/license/Orthtree.h>

#include <iostream>
#include <boost/range/iterator_range.hpp>
#include <CGAL/Orthtree/Traversal_iterator.h>

#include <stack>

namespace CGAL {

/// \cond SKIP_IN_MANUAL
// todo: is this necessary?
// Forward declaration
template <typename T, typename PR, typename PM>
class Orthtree;
/// \endcond

namespace Orthtrees {

/// \cond SKIP_IN_MANUAL

// todo: all of these could be members of Orthtree

template <typename Tree>
const typename Tree::Node* next_sibling(const Tree& orthtree, const typename Tree::Node* n) {

  // todo: maybe this should take a reference?
  if (nullptr == n)
    return nullptr;

  // If this node has no parent, it has no siblings
  if (n->is_root())
    return nullptr;

  // Find out which child this is
  std::size_t index = n->local_coordinates().to_ulong();

  constexpr static int degree = Tree::Node::Degree::value;
  // Return null if this is the last child
  if (int(index) == degree - 1)
    return nullptr;

  // Otherwise, return the next child
  return &(orthtree.children(orthtree.parent(*n))[index + 1]);
}

template <typename Tree>
const typename Tree::Node* next_sibling_up(const Tree& orthtree, const typename Tree::Node* n) {

  if (!n || n->is_root()) return nullptr;

  auto up = &orthtree.parent(*n);
  while (nullptr != up) {

    if (nullptr != next_sibling(orthtree, up))
      return next_sibling(orthtree, up);

    // todo: this could be cleaned up; it's probably not necessary to involve pointers here
    up = up->is_root() ? nullptr : &orthtree.parent(*up);
  }

  return nullptr;
}

template <typename Tree>
const typename Tree::Node* deepest_first_child(const Tree& orthtree, const typename Tree::Node* n) {

  if (n == nullptr)
    return nullptr;

  // Find the deepest child on the left
  auto first = n;
  while (!first->is_leaf())
    first = &orthtree.children(*first)[0];
  return first;
}


template <typename Tree>
const typename Tree::Node* first_child_at_depth(const Tree& orthtree, const typename Tree::Node* n, std::size_t depth) {

  if (!n)
    return nullptr;

  std::stack<const typename Tree::Node*> todo;
  todo.push(n);

  if (n->depth() == depth)
    return n;

  while (!todo.empty()) {
    const typename Tree::Node* node = todo.top();
    todo.pop();

    if (node->depth() == depth)
      return node;

    if (!node->is_leaf())
      for (int i = 0; i < Tree::Node::Degree::value; ++i)
        todo.push(&((*node)[std::size_t(Tree::Node::Degree::value - 1 - i)]));
  }

  return nullptr;
}

/// \endcond

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a preorder traversal.

  A preorder traversal starts from the root towards the leaves.

  \cgalModels `OrthtreeTraversal`
 */
template <typename Tree>
struct Preorder_traversal {
private:

  using Node = typename Tree::Node;

  const Tree& m_orthtree;

public:

  Preorder_traversal(const Tree& orthtree) : m_orthtree(orthtree) {}

  const Node* first() const {
    return &m_orthtree.root();
  }

  const Node* next(const Node* n) const {

    if (n->is_leaf()) {

      auto next = next_sibling(m_orthtree, n);

      if (nullptr == next) {

        return next_sibling_up(m_orthtree, n);
      }

      return next;

    } else {

      // Return the first child of this node
      return &m_orthtree.children(*n)[0];
    }

  }
};

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a postorder traversal.

  A postorder traversal starts from the leaves towards the root.

  \cgalModels `OrthtreeTraversal`
 */
template <typename Tree>
struct Postorder_traversal {
private:

  using Node = typename Tree::Node;

  const Tree& m_orthtree;

public:

  Postorder_traversal(const Tree& orthtree) : m_orthtree(orthtree) {}

  const Node* first() const {
    return deepest_first_child(m_orthtree, m_orthtree.root());
  }

  const Node* next(const Node* n) const {

    auto next = deepest_first_child(m_orthtree, next_sibling(m_orthtree, n));

    if (!next)
      next = n->parent();

    return next;
  }
};

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a traversal on leaves only.

  All non-leave nodes are ignored.

  \cgalModels `OrthtreeTraversal`
 */
template <typename Tree>
struct Leaves_traversal {
private:

  using Node = typename Tree::Node;

  const Tree& m_orthtree;

public:

  Leaves_traversal(const Tree& orthtree) : m_orthtree(orthtree) {}

  const Node* first() const {
    return deepest_first_child(m_orthtree, &m_orthtree.root());
  }

  const Node* next(const Node* n) const {

    auto next = deepest_first_child(m_orthtree, next_sibling(m_orthtree, n));

    if (!next)
      next = deepest_first_child(m_orthtree, next_sibling_up(m_orthtree, n));

    return next;
  }
};

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a traversal of a specific depth level.

  All trees at another depth are ignored. If the selected depth is
  higher than the maximum depth of the orthtree, no node will be traversed.

  \cgalModels `OrthtreeTraversal`
 */
template <typename Tree>
struct Level_traversal {
private:

  using Node = typename Tree::Node;

  const Tree& m_orthtree;
  const std::size_t m_depth;

public:

  /*!
    constructs a `depth`-level traversal.
  */
  Level_traversal(const Tree& orthtree, std::size_t depth) : m_orthtree(orthtree), m_depth(depth) {}

  template <typename Node>
  const Node* first() const {
    return first_child_at_depth(m_orthtree, m_orthtree.root(), m_depth);
  }

  template <typename Node>
  const Node* next(const Node* n) const {
    const Node* next = next_sibling(m_orthtree, n);

    if (!next) {
      const Node* up = n;
      do {
        up = next_sibling_up(m_orthtree, up);
        if (!up)
          return nullptr;

        next = first_child_at_depth(m_orthtree, up, m_depth);
      } while (!next);
    }

    return next;
  }
};

} // Orthtree
} // CGAL

#endif //CGAL_ORTHTREE_TRAVERSALS_H
