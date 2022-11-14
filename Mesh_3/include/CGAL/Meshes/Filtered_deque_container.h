// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Clement JAMIN

#ifndef CGAL_MESHES_FILTERED_DEQUE_CONTAINER_H
#define CGAL_MESHES_FILTERED_DEQUE_CONTAINER_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/algorithm.h>
#include <algorithm>
#include <utility>
#include <deque>
#ifdef CGAL_LINKED_WITH_TBB
  #include <tbb/enumerable_thread_specific.h>
#endif

namespace CGAL {

namespace Meshes {

  /************************************************
  // Class Filtered_deque_container_base
  // Two versions: sequential / parallel
  ************************************************/

  // Sequential
  template <typename Element, typename Quality,
            typename Concurrency_tag>
  class Filtered_deque_container_base
  {
  public:
    typedef std::deque<std::pair<Quality, Element> > Container;
    typedef typename Container::size_type size_type;
    typedef typename Container::value_type value_type;

    void add_to_TLS_lists_impl(bool) {}
    Element get_next_local_element_impl()
    { return Element(); }
    value_type get_next_local_raw_element_impl()
    { return value_type(); }
    void pop_next_local_element_impl() {}

  protected:
    Filtered_deque_container_base() {}
    Filtered_deque_container_base(bool) {}

    template<typename Container>
    void splice_local_lists_impl(Container &)
    {}

    template <typename Predicate>
    bool no_longer_local_element_to_refine_impl(const Predicate &)
    {
      return true;
    }

    template<typename Container>
    void insert_raw_element(const value_type &re, Container &container)
    {
      container.push_back(re);
    }
  };

#ifdef CGAL_LINKED_WITH_TBB
  // Parallel
  template <typename Element, typename Quality>
  class Filtered_deque_container_base<Element, Quality, Parallel_tag>
  {
  public:
    typedef std::deque<std::pair<Quality, Element> > Container;
    typedef typename Container::size_type size_type;
    typedef typename Container::value_type value_type;

    void add_to_TLS_lists_impl(bool add)
    {
      m_add_to_TLS_lists = add;
    }

    // Warning: no_longer_local_element_to_refine_impl must have been called
    // just before calling get_next_local_element_impl
    // (successive calls to "get_next_local_element_impl" are not allowed)
    Element get_next_local_element_impl()
    {
      CGAL_assertion(!m_local_lists.local().empty());
      // Add this? It shouldn't be necessary as user
      // is supposed to call "no_longer_element_to_refine_impl" first
      /*while( !test(container.front()) )
      {
        container.pop_front();
      }*/
      return m_local_lists.local().front().second;
    }

    // Warning: no_longer_local_element_to_refine_impl must have been called
    // just before calling get_next_local_raw_element_impl
    // (successive calls to "get_next_local_raw_element_impl" are not allowed)
    value_type get_next_local_raw_element_impl()
    {
      CGAL_assertion(!m_local_lists.local().empty());
      return m_local_lists.local().front();
    }

    void pop_next_local_element_impl()
    {
      // Erase last element
      m_local_lists.local().pop_front();
    }

  protected:
    Filtered_deque_container_base(bool add_to_TLS_lists = false)
      : m_add_to_TLS_lists(add_to_TLS_lists) {}


    template<typename Container>
    void splice_local_lists_impl(Container &container)
    {
      for(typename LocalList::iterator it_list = m_local_lists.begin() ;
          it_list != m_local_lists.end() ;
          ++it_list )
      {
#ifdef _DEBUG
        size_t deque_size = container.size();
        size_t local_list_size = it_list->size();
#endif
        container.insert(container.end(), it_list->begin(), it_list->end());
        it_list->clear();
      }
    }

    template <typename Predicate>
    bool no_longer_local_element_to_refine_impl(const Predicate &test)
    {
      bool is_empty = m_local_lists.local().empty();
      while( !is_empty && !test(m_local_lists.local().front().second) )
      {
        pop_next_local_element_impl();
        is_empty = m_local_lists.local().empty();
      }
      return is_empty;
    }

    template<typename Container>
    void insert_raw_element(const value_type &re, Container &container)
    {
      if (m_add_to_TLS_lists)
        m_local_lists.local().push_back(re);
      else
        container.push_back(re);
    }


    // === Member variables ===

    typedef tbb::enumerable_thread_specific<
      std::deque<std::pair<Quality, Element> > > LocalList;
    LocalList m_local_lists;
    bool m_add_to_TLS_lists;
  };
#endif // CGAL_LINKED_WITH_TBB

  /************************************************
  // Class Filtered_deque_container
  //
  // This container is a filtered deque:
  // front() and empty() use an object predicate
  // to test if the element is ok.
  ************************************************/

  template <typename Element_, typename Quality_,
            typename Predicate, typename Concurrency_tag>
  class Filtered_deque_container
    : public Filtered_deque_container_base<Element_, Quality_, Concurrency_tag>
  {
  public:
    typedef Filtered_deque_container_base<Element_, Quality_, Concurrency_tag> Base;
    typedef typename Base::Container Container;
    typedef typename Base::value_type value_type;
    typedef typename Base::size_type size_type;
    typedef Quality_ Quality;
    typedef Element_ Element;

  protected:
    // --- protected datas ---
    Container container;
    Predicate test;

    static bool CompareTwoElements(std::pair<Quality, Element> e1,
                                    std::pair<Quality, Element> e2)
    {
      return (e1.first < e2.first);
    }

  public:

    // Constructors - For sequential
    Filtered_deque_container() {}
    explicit Filtered_deque_container(const Predicate &p)
      : test(p) {}

    // Constructors - For parallel
    explicit Filtered_deque_container(bool add_to_TLS_lists)
      : Base(add_to_TLS_lists) {}
    explicit Filtered_deque_container(const Predicate &p, bool add_to_TLS_lists)
      : test(p), Base(add_to_TLS_lists) {}

    void splice_local_lists_impl()
    {
      Base::splice_local_lists_impl(container);
    }

    bool no_longer_local_element_to_refine_impl()
    {
      return Base::no_longer_local_element_to_refine_impl(test);
    }

    void insert_raw_element(const value_type &re)
    {
      Base::insert_raw_element(re, container);
    }

    bool no_longer_element_to_refine_impl()
    {
#ifdef _DEBUG
      size_t deque_size = container.size();
#endif
      bool is_empty = container.empty();
      while( !is_empty && !test(container.front().second) )
      {
        pop_next_element_impl();
        is_empty = container.empty();
      }
      return is_empty;
    }

    // Warning: no_longer_element_to_refine_impl must have been called
    // just before calling get_next_element_impl
    // (successive calls to "get_next_element_impl" are not allowed)
    Element get_next_element_impl() const
    {
      CGAL_assertion(!container.empty());
      // Add this? It shouldn't be necessary as user
      // is supposed to call "no_longer_element_to_refine_impl" first
      /*while( !test(container.front()) )
      {
        container.pop_front();
      }*/
      return container.front().second;
    }

    void add_bad_element(const Element& e, const Quality& q)
    {
      insert_raw_element(std::make_pair(q, e));
    }

    void pop_next_element_impl()
    {
      // Erase last element
      container.pop_front();
    }

    // Sort
    // Worst (smallest) quality first
    void sort ()
    {
      std::sort(container.begin(), container.end(), CompareTwoElements);
    }

    // Clear
    void clear ()
    {
      container.clear();
    }

    // Random shuffle
    void random_shuffle ()
    {
      CGAL::cpp98::random_shuffle(container.begin(), container.end());
    }

    size_type size() const
    {
            return container.size();
    }

    // Warning: no_longer_element_to_refine_impl must have been called
    // just before calling get_next_raw_element_impl
    // (successive calls to "get_next_raw_element_impl" are not allowed)
    value_type get_next_raw_element_impl()
    {
      CGAL_assertion(!container.empty());
      return container.front();
    }

    bool is_zombie(const Element &e) const
    {
      return !test(e);
    }

  }; // end Filtered_deque_container

} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESHES_FILTERED_DEQUE_CONTAINER_H
