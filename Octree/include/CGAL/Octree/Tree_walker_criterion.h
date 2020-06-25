#ifndef OCTREE_TREE_WALKER_CRITERION_H
#define OCTREE_TREE_WALKER_CRITERION_H

namespace CGAL {

  struct Siblings {

    template<class Node>
    Node *operator()(Node *n) {

      // Null handler
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
  };

  struct Depth_first {

    template<class Node>
    Node *operator()(Node *n) {

      // Null handler
      if (nullptr == n)
        return n;

      if (n->isLeaf()) {

        // Check if this node is the last sibling
        if (false) {

          // Return the next sibling

        } else {

          // Return the next sibling of the parent

        }


      } else {

        // Return the first child of this node
      }


      // TODO: placeholder
      return n;

    }
  };
}

#endif //OCTREE_TREE_WALKER_CRITERION_H
