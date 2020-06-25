#ifndef OCTREE_TREE_WALKER_CRITERION_H
#define OCTREE_TREE_WALKER_CRITERION_H

namespace CGAL {

  template<class Node>
  Node *next_sibling(Node *n) {

    // Passing null returns the first node
    if (nullptr == n)
      return n;

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

  struct Depth_first {

    template<class Node>
    Node *first(Node *root) {
      return root;
    }

    template<class Node>
    Node *operator()(Node *n) {

      // Passing null returns the first node
      if (nullptr == n) {

      }

      if (n->is_leaf()) {

        // Check if this node is the last sibling
        if (7 != n->index().to_ulong()) {

          // Return the next sibling
          return (next_sibling(n));

        } else {

          // Return the next sibling of the parent
          // FIXME: this should be able to search upwards at the last sibling
          return (next_sibling(n->parent()));
        }


      } else {

        // Return the first child of this node
        return &(*n)[0];
      }

    }
  };
}

#endif //OCTREE_TREE_WALKER_CRITERION_H
