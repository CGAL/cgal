#ifndef OCTREE_NODE_H
#define OCTREE_NODE_H

#include <array>
#include <memory>
#include <bitset>

namespace CGAL {
  namespace Octree {

    template<typename Value>
    class Node {

    public:

      typedef std::array<Node<Value>, 8> Child_list;
      typedef std::array<uint32_t, 3> Int_location;
      typedef std::bitset<3> Index;

    private:

      Value m_value;

      const Node<Value> *m_parent;
      uint8_t m_depth;

      std::unique_ptr<Child_list> m_children;

    public:

      Node(Value value, Node<Value> *parent = nullptr, Index index = 0) : m_value(value), m_parent(parent) {

        if (parent) {

          m_depth = parent->m_depth + 1;
        }
      }

      Node(Node &&other) {

        m_value = other.m_value;
        m_parent = other.m_parent;
        m_children = std::move(other.m_children);
      }

      // The default constructor is enough

      void split() {

        m_children = std::make_unique<Child_list>();
        for (int child_id = 0; child_id < 8; child_id++) {

          (*m_children)[child_id] = {m_value, this};
        }
      }

      void unsplit() {

        // std::unique_ptr handles this nicely
        m_children.reset();
      }

      Node<Value> &operator[](int index) {
        return (*m_children)[index];
      }

      const Node<Value> &operator[](int index) const {
        return (*m_children)[index];
      }

      const uint8_t &depth() const { return m_depth; }


    };
  }
}

#endif //OCTREE_NODE_H
