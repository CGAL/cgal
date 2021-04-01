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
// Forward declaration
template<typename T, typename PR, typename PM>
class Orthtree;
/// \endcond

namespace Orthtrees {

/// \cond SKIP_IN_MANUAL

template <typename Node>
Node next_sibling(Node n) {

  // Passing null returns the first node
  if (n.is_null())
    return Node();

  // If this node has no parent, it has no siblings
  if (n.parent().is_null())
    return Node();

  // Find out which child this is
  std::size_t index = n.local_coordinates().to_ulong();

  constexpr static int degree = Node::Degree::value;
  // Return null if this is the last child
  if (int(index) == degree - 1)
    return Node();

  // Otherwise, return the next child
  return n.parent()[index + 1];
}

template <typename Node>
Node next_sibling_up(Node n) {

  if (n.is_null())
    return Node();

  Node up = n.parent();

  while (!up.is_null()) {

    if (!next_sibling(up).is_null())
      return next_sibling(up);

    up = up.parent();
  }

  return Node();
}

template <typename Node>
Node deepest_first_child(Node n) {

  if (n.is_null())
    return Node();

  // Find the deepest child on the left
  Node first = n;
  while (!first.is_leaf())
    first = first[0];
  return first;
}

template <typename Node>
Node first_child_at_depth(Node n, std::size_t depth) {

  if (n.is_null())
    return Node();

  std::stack<Node> todo;
  todo.push(n);

  if (n.depth() == depth)
    return n;

  while (!todo.empty())
  {
    Node node = todo.top();
    todo.pop();

    if (node.depth() == depth)
      return node;

    if (!node.is_leaf())
      for (int i = 0; i < Node::Degree::value; ++ i)
        todo.push(node[std::size_t(Node::Degree::value - 1 - i)]);
  }

  return Node();
}

/// \endcond

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a preorder traversal.

  A preorder traversal starts from the root towards the leaves.

  \cgalModels `OrthtreeTraversal`
 */
struct Preorder_traversal {

  template <typename Node>
  Node first(Node root) const {
    return root;
  }

  template <typename Node>
  Node next(Node n) const {

    if (n.is_leaf()) {

      Node next = next_sibling(n);

      if (next.is_null())
        return next_sibling_up(n);

      return next;

    }
    else // Return the first child of this node
      return n[0];
  }
};

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a postorder traversal.

  A postorder traversal starts from the leaves towards the root.

  \cgalModels `OrthtreeTraversal`
 */
struct Postorder_traversal {

  template <typename Node>
  Node first(Node root) const {

    return deepest_first_child(root);
  }

  template <typename Node>
  Node next(Node n) const {

    Node next = deepest_first_child(next_sibling(n));

    if (next.is_null())
      next = n.parent();

    return next;
  }
};

/*!
  \ingroup PkgOrthtreeTraversal
  \brief A class used for performing a traversal on leaves only.

  All non-leave nodes are ignored.

  \cgalModels `OrthtreeTraversal`
 */
struct Leaves_traversal {

  template <typename Node>
  Node first(Node root) const {

    return deepest_first_child(root);
  }

  template <typename Node>
  Node next(Node n) const {

    Node next = deepest_first_child(next_sibling(n));

    if (next.is_null())
      next = deepest_first_child(next_sibling_up(n));

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
struct Level_traversal {

private:

  const std::size_t depth;

public:

  /*!
    constructs a `depth`-level traversal.
  */
  Level_traversal (std::size_t depth) : depth(depth) { }

  template <typename Node>
  Node first(Node root) const {
    return first_child_at_depth(root, depth);
  }

  template <typename Node>
  Node next(Node n) const {
    std::cerr << depth << " ";
    Node next = next_sibling(n);

    if (next.is_null())
    {
      Node up = n;
      do
      {
        up = next_sibling_up(up);
        if (up.is_null())
          return Node();
        next = first_child_at_depth(up, depth);
      }
      while (next.is_null());
    }

    return next;
  }
};

} // Orthtree
} // CGAL

#endif //CGAL_ORTHTREE_TRAVERSALS_H
