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

namespace Orthtrees {

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a preorder traversal.

  \tparam Tree an instance of `Orthtree`

  A preorder traversal starts from the root towards the leaves.

  \cgalModels{OrthtreeTraversal}
 */
template <typename Tree>
struct Preorder_traversal {
private:

  const Tree& m_orthtree;

public:

  using Node_index = typename Tree::Node_index;

  Preorder_traversal(const Tree& orthtree) : m_orthtree(orthtree) {}

  Node_index first_index() const {
    return m_orthtree.root();
  }

  std::optional<Node_index> next_index(Node_index n) const {

    if (m_orthtree.is_leaf(n)) {

      auto next = m_orthtree.next_sibling(n);

      if (!next)
        next = m_orthtree.next_sibling_up(n);

      return next;

    } else {
      return m_orthtree.child(n, 0);
    }
  }

};

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a postorder traversal.

  \tparam Tree an instance of `Orthtree`

  A postorder traversal starts from the leaves towards the root.

  \cgalModels{OrthtreeTraversal}
 */
template <typename Tree>
struct Postorder_traversal {
private:

  const Tree& m_orthtree;

public:

  using Node_index = typename Tree::Node_index;

  Postorder_traversal(const Tree& orthtree) : m_orthtree(orthtree) {}

  Node_index first_index() const {
    return m_orthtree.deepest_first_child(m_orthtree.root());
  }

  std::optional<Node_index> next_index(Node_index n) const {
    return m_orthtree.index(next(&m_orthtree[n]));
  }
};

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a traversal on leaves only.

  \tparam Tree an instance of `Orthtree`

  All non-leaf nodes are ignored.

  \cgalModels{OrthtreeTraversal}
 */
template <typename Tree>
struct Leaves_traversal {
private:

  const Tree& m_orthtree;

public:

  using Node_index = typename Tree::Node_index;

  Leaves_traversal(const Tree& orthtree) : m_orthtree(orthtree) {}

  Node_index first_index() const {
    return m_orthtree.deepest_first_child(m_orthtree.root());
  }

  std::optional<Node_index> next_index(Node_index n) const {

    if (m_orthtree.next_sibling(n))
      return m_orthtree.deepest_first_child(*m_orthtree.next_sibling(n));

    if (m_orthtree.next_sibling_up(n))
      return m_orthtree.deepest_first_child(*m_orthtree.next_sibling_up(n));

    return {};
  }
};

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a traversal of a specific depth level.

  \tparam Tree an instance of `Orthtree`

  All tree nodes at another depth are ignored. If the selected depth is
  higher than the maximum depth of the orthtree, no node will be traversed.

  \cgalModels{OrthtreeTraversal}
 */
template <typename Tree>
struct Level_traversal {
private:

  const Tree& m_orthtree;
  const std::size_t m_depth;

public:

  using Node_index = typename Tree::Node_index;

  /*!
    constructs a `depth`-level traversal.
  */
  Level_traversal(const Tree& orthtree, std::size_t depth) : m_orthtree(orthtree), m_depth(depth) {}

  Node_index first_index() const {
    // assumes the tree has at least one child at m_depth
    return *m_orthtree.first_child_at_depth(m_orthtree.root(), m_depth);
  }

  std::optional<Node_index> next_index(Node_index n) const {

    auto next = m_orthtree.next_sibling(n);

    if (!next) {

      auto up = n;
      do {

        if (!m_orthtree.next_sibling_up(up))
          return {};

        up = *m_orthtree.next_sibling_up(up);
        next = m_orthtree.first_child_at_depth(up, m_depth);

      } while (!next);
    }

    return next;
  }
};

} // Orthtrees
} // CGAL

#endif //CGAL_ORTHTREE_TRAVERSALS_H
