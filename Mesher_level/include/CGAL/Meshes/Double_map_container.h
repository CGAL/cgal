// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifndef CGAL_MESHES_DOUBLE_MAP_CONTAINER_H
#define CGAL_MESHES_DOUBLE_MAP_CONTAINER_H

#include <set>
#include <iostream>
#include <CGAL/Double_map.h>

// backward compatibility
#ifdef CGAL_MESH_3_DEBUG_DOUBLE_MAP
#  define CGAL_MESHES_DEBUG_DOUBLE_MAP CGAL_MESH_3_DEBUG_DOUBLE_MAP
#endif

namespace CGAL {

  namespace Meshes {

    template <typename Elt, class Quality>
    class Double_map_container 
    {
    public:
      typedef Elt Element;

    protected:
      // --- protected datas ---
      Double_map<Element, Quality> m;

    public:
      bool no_longer_element_to_refine_impl() const
      {
        return m.empty();
      }

      Element get_next_element_impl()
      {
        CGAL_assertion(!m.empty());
#if CGAL_MESHES_DEBUG_DOUBLE_MAP
	std::cerr << "get_next_element_impl(" << &*(m.front()->second)
		  << ")\n";
#endif
        return m.front()->second;

      }

      void add_bad_element(const Element& e, const Quality& q)
      {
#if CGAL_MESHES_DEBUG_DOUBLE_MAP
	std::cerr << "add_bad_element(" << &*e << ")\n";
#endif
        m.insert(e, q);
      }

      void pop_next_element_impl()
      {
        m.pop_front();
      }

      void remove_element(const Element& e)
      {
#if CGAL_MESHES_DEBUG_DOUBLE_MAP
	std::cerr << "remove_element(" << &*e << ")\n";
#endif
        m.erase(e);
      }

      typename Double_map<Element, Quality>::size_type
      size() const
      {
	return m.size();
      }
    }; // end Double_map_container
    
  } // end namespace Meshes
} // end namespace CGAL

#endif // CGAL_MESHES_DOUBLE_MAP_CONTAINER_H
