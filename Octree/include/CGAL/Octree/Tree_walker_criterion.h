#ifndef OCTREE_TREE_WALKER_CRITERION_H
#define OCTREE_TREE_WALKER_CRITERION_H

namespace CGAL {

  template<class Node>
  Node *next_sibling(Node *n) {

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

  struct Preorder {

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

        Node *next = next_sibling(n);

        if (nullptr == next) {

          Node *up = n->parent();

          while (nullptr != up) {

            if (nullptr != next_sibling(up))
              return next_sibling(up);

            up = up->parent();
          }

          return nullptr;
        }

        return next;

      } else {

        // Return the first child of this node
        return &(*n)[0];
      }

    }
  };
}

#endif //OCTREE_TREE_WALKER_CRITERION_H
