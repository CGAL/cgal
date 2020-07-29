// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
