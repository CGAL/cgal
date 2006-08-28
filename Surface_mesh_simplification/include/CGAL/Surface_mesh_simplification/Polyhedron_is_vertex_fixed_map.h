// Copyright (c) 2006 Geometry Factory (France).
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
// $URL$
// $Id$
// 
//
// Author(s): Fernando Cacciola <fernando.cacciola@gmail.com>

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLYHEDRON_IS_VERTEX_FIXED_MAP_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLYHEDRON_IS_VERTEX_FIXED_MAP_H

#include <CGAL/boost/graph/Polyhedron_BGL_properties.h>
#include <CGAL/Surface_mesh_simplification/Detail/TSMS_common.h>

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
#  define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS
#else
#  define CGAL_HDS_PARAM_ class HDS
#endif

namespace boost
{

struct vertex_is_fixed_t {} ;

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class Polyhedron_is_vertex_fixed_map : public put_get_helper<bool, Polyhedron_is_vertex_fixed_map<Gt, I, HDS, A> >
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;
  
public:

  typedef readable_property_map_tag category;
  typedef bool value_type;
  typedef bool reference;
  typedef typename graph_traits<Polyhedron const>::vertex_descriptor  key_type;

  reference operator[](key_type const& e) const { e->is_fixed() ; }
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline Polyhedron_is_vertex_fixed_map<Gt, I, HDS, A>
get(vertex_is_fixed_t, const CGAL::Polyhedron_3<Gt,I,HDS,A>& p) 
{
  Polyhedron_is_vertex_fixed_map<Gt, I, HDS, A> m(p);
  return m;
}

template <>
struct Polyhedron_property_map<vertex_is_fixed_t> 
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_ {
    typedef Polyhedron_is_vertex_fixed_map<Gt, I, HDS, A> type;
    typedef Polyhedron_is_vertex_fixed_map<Gt, I, HDS, A> const_type;
  };
};
        
} // namespace boost

#undef CGAL_HDS_PARAM_

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLYHEDRON_IS_VERTEX_FIXED_MAP_H
