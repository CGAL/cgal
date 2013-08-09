// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_MCFSKEL_FIXED_EDGE_MAP_H
#define CGAL_MCFSKEL_FIXED_EDGE_MAP_H

/// @cond CGAL_DOCUMENT_INTERNAL

namespace CGAL {
namespace internal {

// Map used to mark edges as fixed
#include <CGAL/Unique_hash_map.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <boost/graph/graph_traits.hpp>

//
// BGL property map which indicates whether an edge is border OR is marked as non-removable
//
template <class HalfedgeGraph>
class Fixed_edge_map : public boost::put_get_helper<bool, Fixed_edge_map<HalfedgeGraph> >
{
public:

  typedef boost::readable_property_map_tag                                   category;
  typedef bool                                                               value_type;
  typedef bool                                                               reference;
  typedef typename boost::graph_traits<HalfedgeGraph const>::edge_descriptor key_type;

  Fixed_edge_map() : mFixed(false) {}

  reference operator[](key_type const& e) const
  {
    return e->is_border() || is_fixed(e);
  }

  void set_is_fixed (key_type const& e, bool is)
  {
    mFixed[e] = is;
  }

  bool is_fixed(key_type const& e) const
  {
    return mFixed.is_defined(e) ? mFixed[e] : false;
  }

private:

  CGAL::Unique_hash_map<key_type, bool> mFixed;
};

} //namespace internal
} //namespace CGAL

/// @endcond

#endif //CGAL_MCFSKEL_FIXED_EDGE_MAP_H
