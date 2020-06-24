#ifndef OCTREE_TREE_WALKER_CRITERION_H
#define OCTREE_TREE_WALKER_CRITERION_H

namespace CGAL {

  struct Siblings {

    template<class Node>
    Node *operator()(const Node *n) {

      // TODO
    }
  };

  struct Depth_first {

    template<class Node>
    Node *operator()(const Node *n) {

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
