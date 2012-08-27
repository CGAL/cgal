// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_FILTRED_QUEUE_CONTAINER_H
#define CGAL_FILTRED_QUEUE_CONTAINER_H

#include <deque>
#include <queue>

namespace CGAL {

  namespace Meshes {

    /** This container is a filtered queue: front() and empty() use an
        object predicate to test if the element is ok. */
    template <typename Elt, typename Pred>
    class Filtered_queue_container 
    {
    public:
      typedef Elt Element;
      typedef Pred Predicate;

      typedef typename std::deque<Element>::const_iterator const_iterator;
    private:
      // --- private datas ---
      std::deque<Element> d;
      Predicate test;

    public:
      Filtered_queue_container(Predicate p) : d(), test(p) {}
      Filtered_queue_container() : test() {}

      void clear()
      {
        d.clear();
      }

      // backward compatibility with a previous API
      bool empty()
      {
        return no_longer_element_to_refine_impl();
      }

      bool no_longer_element_to_refine_impl()
      {
        if(d.empty())
          return true;

        while( !test(d.front()) )
          {
            d.pop_front();
            if( d.empty() )
              return true;
          }

        return false;
      }

      typename Pred::Result_type get_next_element_impl()
      {
        while( !test(d.front()) )
          d.pop_front();
        return test.result();
      }

      void add_bad_element(const Element& e)
      {
        d.push_back(e);
      }

      // backward compatibility with a previous API
      void remove_next_element()
      {
        return pop_next_element_impl();
      }

      void pop_next_element_impl()
      {
        d.pop_front();
      }

      const_iterator begin() const
      {
        return d.begin();
      }
      
      const_iterator end() const
      {
        return d.end();
      }
    }; // end Simple_queue_container

  } // end namespace Meshes

} // end namespace CGAL

#endif // CGAL_FILTRED_QUEUE_CONTAINER_H
