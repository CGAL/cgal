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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Surface_mesh_simplification/Polyhedron_edge_cached_pointer_map.h $
// $Id: Polyhedron_edge_cached_pointer_map.h 32048 2006-06-23 13:59:36Z lsaboret $
// 
//
// Author(s): Fernando Cacciola <fernando.cacciola@gmail.com>

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLYHEDRON_EDGE_CACHED_POINTER_MAP_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLYHEDRON_EDGE_CACHED_POINTER_MAP_H

#include <CGAL/boost/graph/Polyhedron_BGL_properties.h>
#include <CGAL/Surface_mesh_simplification/Detail/TSMS_common.h>

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
#  define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS
#else
#  define CGAL_HDS_PARAM_ class HDS
#endif

namespace boost
{

struct edge_cached_pointer_t {} ;

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class Polyhedron_edge_cached_pointer_map : public put_get_helper<void*&, Polyhedron_edge_cached_pointer_map<Gt, I, HDS, A> >
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;
  
public:

  typedef lvalue_property_map_tag category;
  typedef void*  value_type;
  typedef void*& reference;
  typedef typename graph_traits<Polyhedron>::edge_descriptor  key_type;

  reference operator[](key_type const& e) const { return e->cached_pointer() ; }
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline Polyhedron_edge_cached_pointer_map<Gt, I, HDS, A>
get(edge_cached_pointer_t, CGAL::Polyhedron_3<Gt,I,HDS,A>& p) 
{
  Polyhedron_edge_cached_pointer_map<Gt, I, HDS, A> m(p);
  return m;
}

template <>
struct Polyhedron_property_map<edge_cached_pointer_t> 
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_ {
    typedef Polyhedron_edge_cached_pointer_map<Gt, I, HDS, A> type;
    typedef Polyhedron_edge_cached_pointer_map<Gt, I, HDS, A> const_type;
  };
};
        
} // namespace boost

#undef CGAL_HDS_PARAM_

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLYHEDRON_EDGE_CACHED_POINTER_MAP_H
