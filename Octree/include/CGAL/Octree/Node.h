#ifndef OCTREE_NODE_H
#define OCTREE_NODE_H

#include <array>
#include <memory>

namespace CGAL {
  namespace Octree {

    template <typename Value>
    class Node {

    public:

      typedef std::array<Node<Value>, 8> Child_list;

    private:

      Value m_value;

      const Node<Value> *m_parent;
      uint8_t m_depth;

      std::unique_ptr<Child_list> m_children;

    public:

      Node(Value value, Node<Value> *parent = nullptr) {

        m_value = value;
        m_parent = parent;
        if (parent)
          m_depth = parent->m_depth + 1;
      }

      void split() {

      }

      void unsplit() {

      }

    };
  }
}

#endif //OCTREE_NODE_H
