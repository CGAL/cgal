// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_FILTRED_QUEUE_CONTAINER_H
#define CGAL_FILTRED_QUEUE_CONTAINER_H

#include <deque>
#include <queue>

namespace CGAL {

  namespace Mesh_2 {

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

      bool empty()
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

      Element& get_next_element()
      {
        while( !test(d.front()) )
          d.pop_front();
        return d.front();
      }

      void add_element(const Element& e)
      {
        d.push_back(e);
      }

      void remove_next_element()
      {
        d.pop_front();
      }

      const_iterator begin() const
      {
        std::cerr << "edges_to_be_conformed.size()= "
                  << d.size() << std::endl;
        return d.begin();
      }
      
      const_iterator end() const
      {
        return d.end();
      }
    }; // end Simple_queue_container

  } // end namespace Mesh_2

} // end namespace CGAL

#endif // CGAL_FILTRED_QUEUE_CONTAINER_H
