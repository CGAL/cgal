// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_MESH_3_SIMPLE_QUEUE_CONTAINER_H
#define CGAL_MESH_3_SIMPLE_QUEUE_CONTAINER_H

#include <queue>

namespace CGAL {

  namespace Mesh_3 {

    template <typename Elt>
    class Simple_queue_container
    {
    public:
      typedef Elt Element;
      typedef std::queue<Element> Queue;
      typedef typename Queue::size_type size_type;

    protected:
      // --- protected datas ---
      Queue q;

    public:
      bool no_longer_element_to_refine_impl() const
      {
        return q.empty();
      }

      Element& get_next_element_impl()
      {
        return q.front();
      }

      void add_bad_element(const Element& e)
      {
        q.push(e);
      }

      void pop_next_element_impl()
      {
        q.pop();
      }

      size_type size() const
      {
        return q.size();
      }
    }; // end Simple_queue_container

  } // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_SIMPLE_QUEUE_CONTAINER_H
