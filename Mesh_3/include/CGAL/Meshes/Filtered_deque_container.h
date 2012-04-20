// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Clement JAMIN

#ifndef CGAL_MESHES_FILTERED_DEQUE_CONTAINER_H
#define CGAL_MESHES_FILTERED_DEQUE_CONTAINER_H

#include <utility>
#include <deque>
#ifdef CONCURRENT_MESH_3
  #include <tbb/enumerable_thread_specific.h>
#endif

namespace CGAL {

  namespace Meshes {
    
    /** This container is a filtered deque: front() and empty() use an
        object predicate to test if the element is ok. */

    template <typename Element_, typename Quality_, 
              typename Predicate>
    class Filtered_deque_container 
    {
    public:
      typedef Quality_ Quality;
      typedef Element_ Element;
      typedef std::deque<std::pair<Quality, Element> > Container;
      typedef typename Container::size_type size_type;
      typedef typename Container::value_type value_type;

    protected:
      // --- protected datas ---
      Container container;
      Predicate test;
#ifdef CONCURRENT_MESH_3
      typedef tbb::enumerable_thread_specific< 
        std::deque<std::pair<Quality, Element> > > LocalList;
      LocalList m_local_lists;
      bool m_add_to_TLS_lists;
#endif

    static bool CompareTwoElements(std::pair<Quality, Element> e1, 
                                   std::pair<Quality, Element> e2) 
    {
      return (e1.first < e2.first);
    }

    public:
		
#ifdef CONCURRENT_MESH_3
      explicit Filtered_deque_container(bool add_to_TLS_lists = false) 
        : m_add_to_TLS_lists(add_to_TLS_lists) {}
      explicit Filtered_deque_container(const Predicate &p, bool add_to_TLS_lists=false)
        : test(p), m_add_to_TLS_lists(add_to_TLS_lists) {}
#else
      Filtered_deque_container() {}
      explicit Filtered_deque_container(const Predicate &p)
        : test(p) {}
#endif

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
        // CJTODO BUG: add this? It shouldn't be necessary as user 
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

#ifdef CONCURRENT_MESH_3
      void add_to_TLS_lists_impl(bool add = true)
      {
        m_add_to_TLS_lists = add;
      }

      void splice_local_lists_impl()
      {
        for( LocalList::iterator it_list = m_local_lists.begin() ; 
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

      
      bool no_longer_local_element_to_refine_impl()
      {
        bool is_empty = m_local_lists.local().empty();
        while( !is_empty && !test(m_local_lists.local().front().second) )
        {
          pop_next_local_element_impl();
          is_empty = m_local_lists.local().empty();
        }
        return is_empty;
      }

      // Warning: no_longer_local_element_to_refine_impl must have been called
      // just before calling get_next_local_element_impl
      // (successive calls to "get_next_local_element_impl" are not allowed)
      Element get_next_local_element_impl()
      {
        CGAL_assertion(!m_local_lists.local().empty());
        // CJTODO BUG: add this? It shouldn't be necessary as user 
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
#endif

      void pop_next_element_impl()
      {
        // Erase last element
        container.pop_front();
      }

      void sort ()
      {
        std::sort(container.begin(), container.end(), CompareTwoElements);
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

      void insert_raw_element(const value_type &re)
      {
#ifdef CONCURRENT_MESH_3
        if (m_add_to_TLS_lists)
          m_local_lists.local().push_back(re);
        else
          container.push_back(re);
#else
        container.push_back(re);
#endif
      }

      bool is_zombie(const Element &e) const
      {
        return !test(e);
      }

    }; // end Filtered_deque_container
    
  } // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESHES_FILTERED_DEQUE_CONTAINER_H
