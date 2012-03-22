// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_MESHES_FILTERED_MULTIMAP_CONTAINER_H
#define CGAL_MESHES_FILTERED_MULTIMAP_CONTAINER_H

#include <map>
#ifdef CONCURRENT_MESH_3
  #include <tbb/enumerable_thread_specific.h>
#endif

namespace CGAL {

  namespace Meshes {
    
    /** This container is a filtered multimap: front() and empty() use an
        object predicate to test if the element is ok. */

    template <typename Element_, typename Quality_, 
              typename Predicate>
    class Filtered_multimap_container 
    {
    public:
      typedef Quality_ Quality;
      typedef Element_ Element;
      typedef std::multimap<Quality, Element> Map;
      typedef typename Map::size_type size_type;
      typedef typename Map::value_type value_type;

    protected:
      // --- protected datas ---
      Map multimap;
      Predicate test;
#ifdef CONCURRENT_MESH_3
      typedef tbb::enumerable_thread_specific< std::vector<std::pair<Quality, Element> > > LocalList;
      LocalList localList;
      bool m_addToTLSLists;
#endif

    public:
		
#ifdef CONCURRENT_MESH_3
      explicit Filtered_multimap_container(bool addToTLSLists = false) 
        : m_addToTLSLists(addToTLSLists) {}
      explicit Filtered_multimap_container(const Predicate &p, bool addToTLSLists=false)
        : test(p), m_addToTLSLists(addToTLSLists) {}
#else
      Filtered_multimap_container() {}
      explicit Filtered_multimap_container(const Predicate &p)
        : test(p) {}
#endif

      bool no_longer_element_to_refine_impl()
      {
        bool is_empty = multimap.empty();
        while( !is_empty && !test(multimap.begin()->second) )
        {
          pop_next_element_impl();
          is_empty = multimap.empty();
        }
        return is_empty;
      }

      // Warning: no_longer_element_to_refine_impl must have been called
      // just before calling get_next_element_impl
      // (successive calls to "get_next_element_impl" are not allowed)
      Element get_next_element_impl() const
      {
        CGAL_assertion(!multimap.empty());
        // CJTODO BUG: add this? It shouldn't be necessary as user 
        // is supposed to call "no_longer_element_to_refine_impl" first
        /*while( !test(multimap.front()) )
        {
          multimap.pop_front();
        }*/
        return multimap.begin()->second;
      }

      void add_bad_element(const Element& e, const Quality& q)
      {
        insert_raw_element(std::make_pair(q, e));
      }

#ifdef CONCURRENT_MESH_3
      void addToTLSLists(bool add = true)
      {
        m_addToTLSLists = add;
      }

      void spliceLocalLists()
      {
        for( LocalList::iterator it_list = localList.begin() ; 
             it_list != localList.end() ; 
             ++it_list )
        {
          multimap.insert(it_list->begin(), it_list->end());
          it_list->clear();
        }
      }
#endif

      void pop_next_element_impl()
      {
        // Erase last element
        multimap.erase( multimap.begin() );
      }

      /*void remove_element(const Element& e)
      {
        multimap.erase(multimap.find(e));
      }

      const Quality& quality(const Element& e)
      {
        return multimap[e];
      }*/

      size_type size() const
      {
	      return multimap.size();
      }

      // Warning: no_longer_element_to_refine_impl must have been called
      // just before calling get_next_element_impl
      // (successive calls to "get_next_element_impl" are not allowed)
      value_type get_next_raw_element_impl()
      {
        CGAL_assertion(!multimap.empty());
        return *multimap.begin();
      }

      void insert_raw_element(const value_type &re)
      {
#ifdef CONCURRENT_MESH_3
        if (m_addToTLSLists)
          localList.local().push_back(re);
        else
          multimap.insert(re);
#else
        multimap.insert(re);
#endif
      }

      bool is_zombie(const Element &e) const
      {
        return !test(e);
      }

    }; // end Filtered_multimap_container
    
  } // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESHES_FILTERED_MULTIMAP_CONTAINER_H
