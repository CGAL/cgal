
#ifndef OCTREE_WALKER_ITERATOR_H
#define OCTREE_WALKER_ITERATOR_H

#include <CGAL/Octree/Tree_walker_criterion.h>

#include <boost/iterator/iterator_facade.hpp>

namespace CGAL {

  template<class Value>
  class Walker_iterator :
          public boost::iterator_facade<Walker_iterator<Value>, Value, boost::forward_traversal_tag> {

  public:

    typedef std::function<Value *(Value *)> Walker_function;

  public:

    Walker_iterator() : m_value(nullptr), m_next(CGAL::Octree::Tree_walker::Preorder()) {}

    Walker_iterator(Value *first, const Walker_function &next) : m_value(first), m_next(next) {}

  private:
    friend class boost::iterator_core_access;

    bool equal(Walker_iterator<Value> const &other) const {
      return m_value == other.m_value;
    }

    void increment() {
      m_value = m_next(m_value);
    }

    Value &dereference() const {
      return *m_value;
    }

  private:

    Value *m_value;
    const Walker_function &m_next;
  };

}

#endif //OCTREE_WALKER_ITERATOR_H
