#ifndef OCTREE_WALKER_CRITERION_H
#define OCTREE_WALKER_CRITERION_H

#include <iostream>
#include "Node.h"

namespace CGAL {

  namespace Octree {

    namespace Walker {

      template<class Value>
      const Node::Node<Value> *next_sibling(const Node::Node<Value> *n) {

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
      const Node::Node<Value> *next_sibling_up(const Node::Node<Value> *n) {

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
      const Node::Node<Value> *deepest_first_child(const Node::Node<Value> *n) {

        if (!n)
          return nullptr;

        // Find the deepest child on the left
        auto first = n;
        while (!first->is_leaf())
          first = &(*first)[0];
        return first;
      }

      struct Preorder {

        template<class Value>
        const Node::Node<Value> *first(const Node::Node<Value> *root) {
          return root;
        }

        template<class Value>
        const Node::Node<Value> *operator()(const Node::Node<Value> *n) {

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

      struct Postorder {

        template<class Value>
        const Node::Node<Value> *first(const Node::Node<Value> *root) {

          return deepest_first_child(root);
        }

        template<class Value>
        const Node::Node<Value> *operator()(const Node::Node<Value> *n) {

          auto next = deepest_first_child(next_sibling(n));

          if (!next)
            next = n->parent();

          return next;
        }
      };

      struct Leaves {

        template<class Value>
        const Node::Node<Value> *first(const Node::Node<Value> *root) {

          return deepest_first_child(root);
        }

        template<class Value>
        const Node::Node<Value> *operator()(const Node::Node<Value> *n) {

          auto next = deepest_first_child(next_sibling(n));

          if (!next)
            next = deepest_first_child(next_sibling_up(n));

          return next;
        }
      };


//  class Tree_walker {
//
//    template<class Kernel, class PointRange>
//    static const Octree_node <Kernel, PointRange> *first(const Octree_node <Kernel, PointRange> *root) {
//      return root;
//    }
//
//    template<class Kernel, class PointRange>
//    static const Octree_node <Kernel, PointRange> *next(const Octree_node <Kernel, PointRange> *node) {
//      return node;
//    }
//  };

//  class Preorder_tree_walker {
//
//  public:
//
//    Preorder_tree_walker() {
//
//    }
//
//    template<class Kernel, class PointRange>
//    static const Octree_node<Kernel, PointRange> *first(const Octree_node <Kernel, PointRange> *root) {
//      return root;
//    }
//
//    template<class Kernel, class PointRange>
//    static const Octree_node <Kernel, PointRange> *next(const Octree_node <Kernel, PointRange> *n) {
//
//      if (n->is_leaf()) {
//
//        auto next = next_sibling(n);
//
//        if (nullptr == next) {
//
//          return next_sibling_up(n);
//        }
//
//        return next;
//
//      } else {
//
//        // Return the first child of this node
//        return &(*n)[0];
//      }
//
//    }
//  };


    }

  }

}

#endif //OCTREE_WALKER_CRITERION_H
