// Copyright (c) 2020  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_SMALL_LIST_H
#define CGAL_SMALL_LIST_H

#include <CGAL/assertions.h>
#include <CGAL/tags.h>

#include <boost/iterator/iterator_facade.hpp>

namespace CGAL
{

namespace internal
{

template <typename T>
class Small_list_node
{
private:
  using Self = Small_list_node<T>;

  T m_value;
  Self* m_prev;
  Self* m_next;

public:

  Small_list_node() : m_prev(nullptr), m_next(nullptr) { }

  Self* prev() const { return m_prev; }
  Self* next() const { return m_next; }
  const T& value() const { return m_value; }
        T& value()       { return m_value; }

  friend void connect (Self* a, Self* b)
  {
    CGAL_assertion (a->m_next == nullptr && b->m_prev == nullptr);
    a->m_next = b;
    b->m_prev = a;
  }

  friend void disconnect (Self* a, Self* b)
  {
    CGAL_assertion (a->m_next == b && b->m_prev == a);
    a->m_next = nullptr;
    b->m_prev = nullptr;
  }

};

template <typename T, std::size_t StackSize>
class Small_list_memory_pool
{
  T m_pool[StackSize];
  std::size_t m_next;

public:

  Small_list_memory_pool() : m_next(0) { }

  Small_list_memory_pool (const Small_list_memory_pool&) = delete;

  T* allocate()
  {
    if (m_next < StackSize)
      return m_pool + (m_next ++);
    return new T;
  }

  void deallocate (T* t)
  {
    if (m_pool <= t && t < m_pool + StackSize)
    {
      // If pointer was the last one allocated, we can use again this
      // chunk of memory
      if (t == m_pool + (StackSize - 1))
        m_next = StackSize - 1;
    }
    else
      delete t;
  }
};

template<typename T, typename ReverseTag>
class Small_list_iterator
  : public boost::iterator_facade<Small_list_iterator<T, ReverseTag>,
                                  T,
                                  std::bidirectional_iterator_tag>
{
  using Self = Small_list_iterator<T, ReverseTag>;
  using Node = Small_list_node<T>;

public:
  Small_list_iterator(Node* last = nullptr,
                      Node* node = nullptr)
    : m_node (node), m_last (last) { }

  Node* node() { return m_node; }

private:
  friend class boost::iterator_core_access;
  void increment() { increment(ReverseTag()); }
  void decrement()
  {
    if (m_node == nullptr)
      m_node = m_last;
    else
      decrement(ReverseTag());
  }

  // forward iterator
  void increment (const Tag_false&) { m_node = m_node->next(); }
  void decrement (const Tag_false&) { m_node = m_node->prev(); }

  // reverse iterator
  void increment (const Tag_true&) { m_node = m_node->prev(); }
  void decrement (const Tag_true&) { m_node = m_node->next(); }

  bool equal(const Self& other) const { return (this->m_node == other.m_node); }
  T& dereference() const { return const_cast<T&>(m_node->value()); }

  Node* m_node;
  Node* m_last;
};

} // namespace internal

template <typename T, std::size_t StackSize>
class Small_list
{
private:

  using Node = internal::Small_list_node<T>;
  using Pool = internal::Small_list_memory_pool<Node, StackSize>;

public:

  using iterator = internal::Small_list_iterator<T, Tag_false>;
  using reverse_iterator = internal::Small_list_iterator<T, Tag_true>;
  using const_iterator = iterator;
  using value_type = T;

private:

  Pool m_pool;
  Node* m_first;
  Node* m_last;
  std::size_t m_size;

public:

  Small_list()
    : m_first (nullptr)
    , m_last (nullptr)
    , m_size (0) { }

  Small_list(const Small_list& other)
    : m_first (nullptr)
    , m_last (nullptr)
    , m_size (0)
  {
    std::copy (other.begin(), other.end(),
               std::back_inserter (*this));
  }

  ~Small_list() { clear(); }

  bool empty() const { return (m_size == 0); }
  std::size_t size() const { return m_size; }

  void push_back (const T& t)
  {
    Node* n = m_pool.allocate();
    n->value() = t;
    if (empty())
    {
      m_first = n;
      m_last = n;
    }
    else
    {
      connect (m_last, n);
      m_last = n;
    }
    ++ m_size;
  }

  void push_front (const T& t)
  {
    Node* n = m_pool.allocate();
    n->value() = t;
    if (empty())
    {
      m_first = n;
      m_last = n;
    }
    else
    {
      connect (n, m_first);
      m_first = n;
    }
    ++ m_size;
  }

  const_iterator begin() const { return const_iterator (m_last, m_first); }
  const_iterator end() const { return const_iterator (m_last); }
  iterator begin() { return iterator (m_last, m_first); }
  iterator end() { return iterator (m_last); }
  reverse_iterator rbegin() { return reverse_iterator (m_first, m_last); }
  reverse_iterator rend() { return reverse_iterator (m_first); }

  iterator insert (iterator it, const T& t)
  {
    Node* next = it.node();

    if (next == m_first)
      push_front(t);
    else if (next == nullptr)
      push_back(t);
    else
    {
      Node* n = m_pool.allocate();
      n->value() = t;
      Node* prev = next->prev();
      disconnect (prev, next);
      connect (prev, n);
      connect (n, next);
      ++ m_size;
    }

    return next->prev();
  }

  iterator erase(iterator first, iterator last)
  {
    if (first == last)
      return last;

    Node* firstn = first.node();
    CGAL_assertion (firstn != nullptr);

    Node* beyondn = last.node();
    Node* lastn = (beyondn == nullptr ? m_last : beyondn->prev());

    Node* n = firstn;
    while (n != lastn)
    {
      Node* next = n->next();
      disconnect (n, n->next());
      if (n != firstn)
      {
        m_pool.deallocate(n);
        -- m_size;
      }
      n = next;
    }

    if (lastn != m_last)
    {
      disconnect (lastn, beyondn);
      connect (firstn, beyondn);
    }
    else
      m_last = firstn;

    m_pool.deallocate (lastn);
    -- m_size;

    return erase (first);
  }

  iterator erase(iterator it)
  {
    Node* n = it.node();
    CGAL_assertion (n != nullptr);

    Node* out = nullptr;
    if (m_size == 1)
    {
      m_first = nullptr;
      m_last = nullptr;
    }
    else if (n == m_first)
    {
      m_first = n->next();
      disconnect(n, m_first);
      out = m_first;
    }
    else if (n == m_last)
    {
      m_last = n->prev();
      disconnect (m_last, n);
    }
    else
    {
      Node* prev = n->prev();
      Node* next = n->next();
      disconnect (prev, n);
      disconnect (n, next);
      connect (prev, next);
      out = next;
    }

    m_pool.deallocate(n);
    -- m_size;
    return out;
  }

  void clear()
  {
    if (m_size == 0)
      return;

    Node* n = m_first;
    while (n != nullptr)
    {
      Node* next = n->next();
      m_pool.deallocate(n);
      n = next;
    }

    m_first = nullptr;
    m_last = nullptr;
    m_size = 0;
  }

  const T& front() const { return m_first->value(); }
        T& front()       { return m_first->value(); }
  const T& back() const { return m_last->value(); }
        T& back()       { return m_last->value(); }

};

} // namespace CGAL

#endif // CGAL_SMALL_LIST_H
