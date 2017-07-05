// Copyright (c) 2016  GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_BOOL_PROPERTY_MAP_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_BOOL_PROPERTY_MAP_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/Kernel_traits.h>

namespace CGAL {

namespace Surface_mesh_parameterization {

namespace internal {

template <typename Set>
class Bool_property_map
{
  typedef Set S;
  typedef Bool_property_map<S> Self;

public:
  typedef typename Set::key_type key_type;
  typedef bool value_type;
  typedef bool reference;
  typedef boost::read_write_property_map_tag category;

  friend value_type get(const Self& pm, const key_type& k)
  {
    return pm.m_s->find(k) != pm.m_s->end();
  }

  friend void put(const Self& pm, key_type& k, const value_type& v)
  {
    if(v) {
      pm.m_s->insert(k);
    } else {
      pm.m_s->erase(k);
    }
  }

  Bool_property_map() : m_s(0) { }
  Bool_property_map(S& s) : m_s(&s) { }

private:
  S* m_s;
};

} // namespace internal

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_INTERNAL_BOOL_PROPERTY_MAP_H
