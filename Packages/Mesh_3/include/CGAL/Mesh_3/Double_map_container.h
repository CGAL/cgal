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

#ifndef CGAL_MESH_3_DOUBLE_MAP_CONTAINER_H
#define CGAL_MESH_3_DOUBLE_MAP_CONTAINER_H

#include <set>
#include <CGAL/Double_map.h>

namespace CGAL {

  namespace Mesh_3 {

    /** Example of a \c NON-REMOVABLE container:
        one cannot remove elements from it, but the front. */
    template <typename Elt, class Quality>
    class Double_map_container 
    {
    public:
      typedef Elt Element;

    private:
      // --- private datas ---
      Double_map<Element, Quality> m;

    public:
      bool no_longer_element_to_refine_impl() const
      {
        return m.empty();
      }

      Element get_next_element_impl()
      {
        CGAL_assertion(!m.empty());
#if CGAL_MESH_3_DEBUG_DOUBLE_MAP
	std::cerr << "get_next_element_impl(" << &*(m.front()->second)
		  << ")\n";
#endif
        return m.front()->second;

      }

      void add_element(const Element& e, const Quality& q)
      {
#if CGAL_MESH_3_DEBUG_DOUBLE_MAP
	std::cerr << "add_element(" << &*e << ")\n";
#endif
        m.insert(e, q);
      }

      void pop_next_element_impl()
      {
#if CGAL_MESH_3_DEBUG_DOUBLE_MAP
	std::cerr << "pop_front(" << &*(m.front()->second) << ")\n";
#endif
        m.pop_front();
      }

      void remove_element(Element& e)
      {
#if CGAL_MESH_3_DEBUG_DOUBLE_MAP
	std::cerr << "remove_element(" << &*e << ")\n";
#endif
        m.erase(e);
      }

      int queue_size()
      {
	return m.size();
      }

    }; // end Double_map_container
    
  }; // end namespace Mesh_3
}; // end namespace CGAL

#endif // CGAL_MESH_3_DOUBLE_MAP_CONTAINER_H
