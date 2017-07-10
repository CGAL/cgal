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

#ifndef CGAL_MESH_3_PROPERTIES_SURFACE_MESH_H
#define CGAL_MESH_3_PROPERTIES_SURFACE_MESH_H

#include <CGAL/license/Mesh_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Mesh_3/properties.h>
#include <set>

namespace boost {

template <typename P, typename I>
struct property_map<CGAL::Surface_mesh<P>, CGAL::face_patch_id_t<I> >
{

  typedef typename boost::graph_traits<CGAL::Surface_mesh<P> >::face_descriptor face_descriptor;

  typedef typename CGAL::Surface_mesh<P>::template Property_map<face_descriptor, I> type;
  typedef type const_type;
};


template <typename P>
struct property_map<CGAL::Surface_mesh<P>, CGAL::face_patch_id_t<void> >
{

  typedef typename boost::graph_traits<CGAL::Surface_mesh<P> >::face_descriptor face_descriptor;

  typedef CGAL::Static_property_map<typename boost::graph_traits<CGAL::Surface_mesh<P> >::face_descriptor,std::pair<int,int> > type;
  typedef type const_type;
};


template<typename P>
struct property_map<CGAL::Surface_mesh<P>, CGAL::halfedge_is_feature_t>
{
  typedef typename boost::graph_traits<CGAL::Surface_mesh<P> >::halfedge_descriptor halfedge_descriptor;

  typedef typename CGAL::Surface_mesh<P>::template Property_map<halfedge_descriptor, bool> type;
  typedef type const_type;
};


template <typename P>
struct property_map<CGAL::Surface_mesh<P>, CGAL::vertex_num_feature_edges_t>
{

  typedef typename boost::graph_traits<CGAL::Surface_mesh<P> >::vertex_descriptor vertex_descriptor;

  typedef typename CGAL::Surface_mesh<P>::template Property_map<vertex_descriptor, int> type;
  typedef type const_type;
};


template <typename P, typename I>
struct property_map<CGAL::Surface_mesh<P>, CGAL::vertex_incident_patches_t<I> >
{

  typedef typename boost::graph_traits<CGAL::Surface_mesh<P> >::vertex_descriptor vertex_descriptor;

  typedef typename CGAL::Surface_mesh<P>::template Property_map<vertex_descriptor, std::set<I> > type;
  typedef type const_type;
};


template <typename P>
struct property_map<CGAL::Surface_mesh<P>, CGAL::vertex_time_stamp_t>
{

  typedef typename boost::graph_traits<CGAL::Surface_mesh<P> >::vertex_descriptor vertex_descriptor;

  typedef typename CGAL::Surface_mesh<P>::template Property_map<vertex_descriptor, std::size_t> type;
  typedef type const_type;
};


template <typename P>
struct property_map<CGAL::Surface_mesh<P>, CGAL::halfedge_time_stamp_t>
{

  typedef typename boost::graph_traits<CGAL::Surface_mesh<P> >::halfedge_descriptor halfedge_descriptor;

  typedef typename CGAL::Surface_mesh<P>::template Property_map<halfedge_descriptor, std::size_t> type;
  typedef type const_type;
};


template <typename P>
struct property_map<CGAL::Surface_mesh<P>, CGAL::face_time_stamp_t>
{

  typedef typename boost::graph_traits<CGAL::Surface_mesh<P> >::face_descriptor face_descriptor;

  typedef typename CGAL::Surface_mesh<P>::template Property_map<face_descriptor, std::size_t> type;
  typedef type const_type;
};

} // namespace boost

namespace CGAL {

template <typename P, typename Property_tag>
struct Get_pmap_of_surface_mesh {
  typedef typename boost::property_map<Surface_mesh<P>, Property_tag >::type type;
};

template <typename P, typename I>
typename boost::lazy_disable_if<
   boost::is_const<P>,
   Get_pmap_of_surface_mesh<P, CGAL::face_patch_id_t<I> >
 >::type
inline get(CGAL::face_patch_id_t<I>, Surface_mesh<P> & smesh)
{
 typedef typename boost::graph_traits<Surface_mesh<P> >::face_descriptor face_descriptor;
  return smesh. template add_property_map<face_descriptor,I>("f:patch_id", 1).first;
}


#define CGAL_MESH_3_PROPERTY_SURFACE_MESH_RETURN_TYPE(Tag) \
  typename boost::lazy_disable_if<                      \
     boost::is_const<P>,                                \
     Get_pmap_of_surface_mesh<P, Tag >                  \
   >::type

template <typename P>
CGAL_MESH_3_PROPERTY_SURFACE_MESH_RETURN_TYPE(CGAL::face_patch_id_t<void>)
inline get(CGAL::face_patch_id_t<void>, Surface_mesh<P> &)
{
  typedef CGAL::Static_property_map<typename boost::graph_traits<Surface_mesh<P> >::face_descriptor,std::pair<int,int> > Pmap;

  return Pmap(std::make_pair(0,1));
}


template <typename P>
CGAL_MESH_3_PROPERTY_SURFACE_MESH_RETURN_TYPE(CGAL::halfedge_is_feature_t)
inline get(CGAL::halfedge_is_feature_t, Surface_mesh<P>& smesh)
{
  typedef typename boost::graph_traits<Surface_mesh<P> >::halfedge_descriptor halfedge_descriptor;
  return smesh. template add_property_map<halfedge_descriptor,bool>("h:is_feature").first;
}


 template <typename P>
CGAL_MESH_3_PROPERTY_SURFACE_MESH_RETURN_TYPE(CGAL::vertex_num_feature_edges_t)
 inline get(CGAL::vertex_num_feature_edges_t, Surface_mesh<P> & smesh)
{
  typedef typename boost::graph_traits<Surface_mesh<P> >::vertex_descriptor vertex_descriptor;
  return smesh. template add_property_map<vertex_descriptor,int>("v:nfe").first;
}

  template <typename P, typename I>
typename boost::lazy_disable_if<
   boost::is_const<P>,
   Get_pmap_of_surface_mesh<P, CGAL::vertex_incident_patches_t<I> >
 >::type
  inline get(CGAL::vertex_incident_patches_t<I>, Surface_mesh<P> & smesh)
{
  typedef typename boost::graph_traits<Surface_mesh<P> >::vertex_descriptor vertex_descriptor;
  return smesh. template add_property_map<vertex_descriptor,std::set<I> >("v:ip").first;
}


 template <typename P>
 CGAL_MESH_3_PROPERTY_SURFACE_MESH_RETURN_TYPE(CGAL::vertex_time_stamp_t)
inline get(CGAL::vertex_time_stamp_t, Surface_mesh<P> & smesh)
{
  typedef typename boost::graph_traits<Surface_mesh<P> >::vertex_descriptor vertex_descriptor;
  return smesh. template add_property_map<vertex_descriptor,std::size_t>("v:time_stamp").first;
}


 template <typename P>
 CGAL_MESH_3_PROPERTY_SURFACE_MESH_RETURN_TYPE(CGAL::halfedge_time_stamp_t)
inline get(CGAL::halfedge_time_stamp_t, Surface_mesh<P> & smesh)
{
  typedef typename boost::graph_traits<Surface_mesh<P> >::halfedge_descriptor halfedge_descriptor;
  return smesh. template add_property_map<halfedge_descriptor,std::size_t>("h:time_stamp").first;
}


 template <typename P>
 CGAL_MESH_3_PROPERTY_SURFACE_MESH_RETURN_TYPE(CGAL::face_time_stamp_t)
inline get(CGAL::face_time_stamp_t, Surface_mesh<P> & smesh)
{
  typedef typename boost::graph_traits<Surface_mesh<P> >::face_descriptor face_descriptor;
  return smesh. template add_property_map<face_descriptor,std::size_t>("v:time_stamp").first;
}

} // namespace CGAL

#undef CGAL_MESH_3_PROPERTY_SURFACE_MESH_RETURN_TYPE

#endif //CGAL_MESH_3_PROPERTIES_SURFACE_MESH_H
