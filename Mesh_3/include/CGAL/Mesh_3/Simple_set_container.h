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

#ifndef CGAL_MESH_3_SIMPLE_SET_CONTAINER_H
#define CGAL_MESH_3_SIMPLE_SET_CONTAINER_H

#include <set>

namespace CGAL {

  namespace Mesh_3 {

    /** Example of a \c REMOVABLE container:
        one can remove any elements from it. */
    template <typename Element>
    class Simple_set_container 
    {
      // --- private datas ---
      std::set<Element> s;

    public:
      bool empty() const
      {
        return s.empty();
      }

      Element get_next_element()
      {
        CGAL_assertion(!s.empty());
        
        return *(s.begin());

      }

      void add_element(const Element& e)
      {
        s.insert(e);
      }

      void remove_next_element()
      {
        s.erase(s.begin());
      }

      void remove_element(const Element& e)
      {
        s.erase(e);
      }
    }; // end Simple_set_container
    
  }; // end namespace Mesh_3
}; // end namespace CGAL

#endif // CGAL_MESH_3_SIMPLE_SET_CONTAINER_H
