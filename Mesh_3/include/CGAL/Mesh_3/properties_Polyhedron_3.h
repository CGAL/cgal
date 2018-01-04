// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_MESH_3_PROPERTIES_POLYHEDRON_3_H
#define CGAL_MESH_3_PROPERTIES_POLYHEDRON_3_H

#include <CGAL/license/Mesh_3.h>
#include <CGAL/Mesh_3/properties.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS

namespace CGAL {

struct Polyhedron_face_time_stamp_pmap
{
  typedef void                               key_type;
  typedef std::size_t                        value_type;
  typedef std::size_t                        reference;
  typedef boost::read_write_property_map_tag category;
};

template <typename Handle_type>
std::size_t get(Polyhedron_face_time_stamp_pmap, Handle_type h)
{
  return h->time_stamp();
}

template <typename Handle_type>
void put(Polyhedron_face_time_stamp_pmap, Handle_type h,
         std::size_t ts)
{
  h->set_time_stamp(ts);
}

template <>
struct Polyhedron_property_map<CGAL::vertex_time_stamp_t>
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_
  {
    typedef Polyhedron_face_time_stamp_pmap type;
    typedef type const_type;
  };
};

template <>
struct Polyhedron_property_map<CGAL::halfedge_time_stamp_t>
  : public Polyhedron_property_map<CGAL::vertex_time_stamp_t>
{};

template <>
struct Polyhedron_property_map<CGAL::face_time_stamp_t>
  : public Polyhedron_property_map<CGAL::vertex_time_stamp_t>
{};


template <typename Patch_id>
struct Polyhedron_face_patch_id_pmap {
  typedef void                               key_type;
  typedef Patch_id                           value_type;
  typedef Patch_id                           reference;
  typedef boost::read_write_property_map_tag category;
};

template <typename Patch_id, typename Handle_type>
Patch_id get(Polyhedron_face_patch_id_pmap<Patch_id>, Handle_type h)
{
  return h->patch_id();
}

template <typename Patch_id, typename Handle_type>
void put(Polyhedron_face_patch_id_pmap<Patch_id>, Handle_type h,
         Patch_id pid)
{
  h->set_patch_id(pid);
}

template <typename Patch_id>
struct Polyhedron_property_map<CGAL::face_patch_id_t<Patch_id> >
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_
  {
    typedef Polyhedron_face_patch_id_pmap<Patch_id> type;
    typedef type const_type;
  };
};



namespace internal{
BOOST_MPL_HAS_XXX_TRAIT_DEF(Plane_3)

template <class Gt, class I, CGAL_HDS_PARAM_, class A>
struct Get_static_property_map {
  typedef boost::graph_traits<CGAL::Polyhedron_3<Gt,I,HDS,A> > Graph_traits;
  typedef CGAL::Static_property_map<typename Graph_traits::face_descriptor,
                                    std::pair<int,int> > type;
};

} // end namespace internal

template <class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::lazy_enable_if<
  internal::has_Plane_3<Gt>,
  internal::Get_static_property_map<Gt, I, HDS, A>
  >::type
get(CGAL::face_patch_id_t<void>, const Polyhedron_3<Gt,I,HDS,A>&)
{
  typedef typename internal::Get_static_property_map<Gt, I, HDS, A>::type Pmap;
  return Pmap( std::pair<int,int>(0,1) );
}

template <>
struct Polyhedron_property_map<CGAL::face_patch_id_t<void> >
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_
  {
    typedef typename internal::Get_static_property_map<Gt,I,HDS,A>::type type;
    typedef type const_type;
  };
};



struct Polyhedron_num_feature_edges_pmap {
  typedef void                               key_type;
  typedef int                                value_type;
  typedef int                                reference;
  typedef boost::read_write_property_map_tag category;
};

template <typename Handle_type>
int get(Polyhedron_num_feature_edges_pmap, Handle_type h)
{
  return h->nb_of_feature_edges;
}

template <typename Handle_type>
void put(Polyhedron_num_feature_edges_pmap, Handle_type h, int n)
{
  h->nb_of_feature_edges = n;
}


template <>
struct Polyhedron_property_map<CGAL::vertex_num_feature_edges_t>
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_
  {
    typedef Polyhedron_num_feature_edges_pmap type;
    typedef type const_type;
  };
};

struct Polyhedron_is_feature_edge_pmap {
  typedef void                               key_type;
  typedef bool                               value_type;
  typedef bool                               reference;
  typedef boost::read_write_property_map_tag category;
};

template <typename Handle_type>
bool get(Polyhedron_is_feature_edge_pmap, Handle_type h)
{
  return h->is_feature_edge();
}

template <typename Handle_type>
void put(Polyhedron_is_feature_edge_pmap, Handle_type h, bool b)
{
  h->set_feature_edge(b);
}

template <>
struct Polyhedron_property_map<CGAL::halfedge_is_feature_t>
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_
  {
    typedef Polyhedron_is_feature_edge_pmap type;
    typedef type const_type;
  };
};



template <typename Patch_id>
struct Polyhedron_incident_patches_pmap {
  typedef void                               key_type;
  typedef std::set<Patch_id>                 value_type;
  typedef std::set<Patch_id>&                reference;
  typedef boost::read_write_property_map_tag category;
};

template <typename Patch_id, typename Handle_type>
std::set<Patch_id>& get(Polyhedron_incident_patches_pmap<Patch_id>,
                        Handle_type h)
{
  return h->incident_patches_ids_set();
}

template <typename Patch_id, typename Handle_type>
void put(Polyhedron_incident_patches_pmap<Patch_id>,
         Handle_type h, const std::set<Patch_id>& v)
{
  BOOST_FOREACH(Patch_id n, v)
    h->add_incident_patch(n);
}

template <typename Patch_id>
struct Polyhedron_property_map<CGAL::vertex_incident_patches_t<Patch_id> >
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_
  {
    typedef Polyhedron_incident_patches_pmap<Patch_id> type;
    typedef type const_type;
  };
};


} // end namespace CGAL

#undef CGAL_HDS_PARAM_

#endif // CGAL_MESH_3_PROPERTIES_POLYHEDRON_3_H
