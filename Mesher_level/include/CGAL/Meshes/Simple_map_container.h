// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_MESHES_SIMPLE_MAP_CONTAINER_H
#define CGAL_MESHES_SIMPLE_MAP_CONTAINER_H

#include <map>

namespace CGAL {

  namespace Meshes {

    template <typename Element, typename Info>
    class Simple_map_container
    {
    public:
      typedef std::map<Element, Info> Map;
      typedef typename Map::size_type size_type;
      typedef typename Map::value_type value_type;

    protected:
      // --- protected datas ---
      Map map;

    public:
      bool no_longer_element_to_refine_impl() const
      {
        return map.empty();
      }

      const value_type& get_next_element_impl()
      {
        CGAL_assertion(!map.empty());

        return *(map.begin());
      }

      void add_bad_element(const Element& e, const Info& i)
      {
        map.insert(std::make_pair(e, i));
      }

      void pop_next_element_impl()
      {
        map.erase(map.begin());
      }

      void remove_element(const Element& e)
      {
        map.erase(e);
      }

      const Info& info(const Element& e)
      {
        return map[e];
      }

      size_type size() const
      {
        return map.size();
      }
    }; // end Simple_map_container

  } // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESHES_SIMPLE_MAP_CONTAINER_H
