// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_MESH_3_SIMPLE_QUEUE_CONTAINER_H
#define CGAL_MESH_3_SIMPLE_QUEUE_CONTAINER_H

#include <queue>

namespace CGAL {

  namespace Mesh_3 {

    /** Example of a \c NON-REMOVABLE container:
        one cannot remove elements from it, but the front. */
    template <typename Elt>
    class Simple_queue_container 
    {
    public:
      typedef Elt Element;

    private:
      // --- private datas ---
      std::queue<Element> q;

    public:
      bool no_longer_element_to_refine_impl() const
      {
        return q.empty();
      }

      Element& get_next_element_impl()
      {
        return q.front();
      }

      void add_element(const Element& e)
      {
        q.push(e);
      }

      void remove_next_element_impl()
      {
        q.pop();
      }
    }; // end Simple_queue_container
    
  }; // end namespace Mesh_3
}; // end namespace CGAL

#endif // CGAL_MESH_3_SIMPLE_QUEUE_CONTAINER_H
