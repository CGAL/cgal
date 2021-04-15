// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri


#ifndef CGAL_PROPERTIES_SEAM_MESH_H
#define CGAL_PROPERTIES_SEAM_MESH_H

#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/boost/graph/properties.h>

namespace CGAL {

template <class TM,class SEM, class SVM>
class Seam_mesh;

template<class TM,class SEM, class SVM>
class Seam_mesh_point_map
{
  typedef typename boost::property_map<TM,vertex_point_t>::const_type Graph_pmap;
public:
  typedef boost::readable_property_map_tag                          category;
  typedef typename boost::property_traits<Graph_pmap>::value_type   value_type;
  typedef typename boost::property_traits<Graph_pmap>::reference    reference;
  typedef typename boost::graph_traits<Seam_mesh<TM, SEM, SVM> >::vertex_descriptor
                                                                    key_type;

public:
  Seam_mesh_point_map(const Seam_mesh<TM, SEM, SVM>& mesh, Graph_pmap map)
    : mesh(mesh), map(map)
  { }

  friend reference get(const Seam_mesh_point_map& pmap, key_type vd)
  {
    typename boost::graph_traits<TM>::halfedge_descriptor hd = vd;
    return get(pmap.map, target(hd, pmap.mesh.mesh()));
  }

private:
  const Seam_mesh<TM, SEM, SVM>& mesh;
  Graph_pmap map;
};

template<class TM, class SEM, class SVM, class Map>
class Seam_mesh_uv_map
{
  typedef Seam_mesh_uv_map<TM, SEM, SVM, Map>                   Self;
public:
  typedef boost::read_write_property_map_tag                    category;
  typedef typename boost::property_traits<Map>::value_type      value_type;
  typedef typename boost::property_traits<Map>::reference       reference;
  typedef typename boost::graph_traits<Seam_mesh<TM, SEM, SVM> >::vertex_descriptor
                                                                key_type;
  // assert that key_type equals boost::property_traits<Map>::key_type

    typedef Seam_mesh<TM, SEM, SVM> Mesh;
public:

  Seam_mesh_uv_map(const Seam_mesh<TM, SEM, SVM>& mesh, Map map)
    : mesh(mesh), map(map)
  { }

    Seam_mesh_uv_map(const Self& other)
    : mesh(other.mesh), map(other.map)
  { }

    //reference operator[](const key_type& vd) const
    //{
    //  typename boost::graph_traits<TM, SEM, SVM>::halfedge_descriptor hd = vd;
    //  return map[target(hd, mesh.mesh())];
    //}

  inline friend reference get(const Self& pm, key_type vd)
  {
    typename boost::graph_traits<TM>::halfedge_descriptor hd = vd;
    return get(pm.map, target(hd, pm.mesh.mesh()));
  }

  inline friend void put(const Self& pm, key_type vd, const value_type& uv)
  {
    typename boost::graph_traits<TM>::halfedge_descriptor hd = vd;
    put(pm.map, target(hd, pm.mesh.mesh()), uv);
  }

private:
  const Seam_mesh<TM, SEM, SVM>& mesh;
  Map map;
};

} // namespace CGAL

// overloads and specializations in the boost namespace
namespace boost {

template<class TM, class SEM, class SVM, typename T>
struct property_map<CGAL::Seam_mesh<TM, SEM, SVM>, T>
  : public cgal_no_property
{ };

template<class TM, class SEM, class SVM>
struct property_map<CGAL::Seam_mesh<TM, SEM, SVM>, CGAL::vertex_point_t>
{
  typedef CGAL::Seam_mesh<TM, SEM, SVM>                  SM;
  typedef CGAL::Seam_mesh_point_map<TM, SEM, SVM>        type;
  typedef type                                           const_type;
};

template <class TM, class SEM, class SVM, typename T>
struct property_map<CGAL::Seam_mesh<TM, SEM, SVM>, CGAL::dynamic_vertex_property_t<T> >
{
  typedef typename boost::graph_traits<CGAL::Seam_mesh<TM, SEM, SVM> >::vertex_descriptor vertex_descriptor;
  typedef CGAL::internal::Dynamic_property_map<vertex_descriptor,T> type;
  typedef type const_type;
};

template <class TM, class SEM, class SVM, typename T>
struct property_map<CGAL::Seam_mesh<TM, SEM, SVM>, CGAL::dynamic_halfedge_property_t<T> >
{
  typedef typename boost::graph_traits<CGAL::Seam_mesh<TM, SEM, SVM> >::halfedge_descriptor halfedge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<halfedge_descriptor,T> type;
  typedef type const_type;
};


template <class TM, class SEM, class SVM, typename T>
struct property_map<CGAL::Seam_mesh<TM, SEM, SVM>, CGAL::dynamic_edge_property_t<T> >
{
  typedef typename boost::graph_traits<CGAL::Seam_mesh<TM, SEM, SVM> >::edge_descriptor edge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<edge_descriptor,T> type;
  typedef type const_type;
};

template <class TM, class SEM, class SVM, typename T>
struct property_map<CGAL::Seam_mesh<TM, SEM, SVM>, CGAL::dynamic_face_property_t<T> >
{
  typedef typename boost::graph_traits<CGAL::Seam_mesh<TM, SEM, SVM> >::face_descriptor face_descriptor;
  typedef CGAL::internal::Dynamic_property_map<face_descriptor,T> type;
  typedef type const_type;
};
} // namespace boost

namespace CGAL {

template<class TM, class SEM, class SVM>
typename boost:: property_map<CGAL::Seam_mesh<TM, SEM, SVM>, vertex_point_t>::const_type
get(vertex_point_t, const Seam_mesh<TM, SEM, SVM>& sm)
{
  return Seam_mesh_point_map<TM, SEM, SVM>(sm, get(vertex_point, sm.mesh()));
}

template<class TM, class SEM, class SVM>
struct graph_has_property<CGAL::Seam_mesh<TM, SEM, SVM>, CGAL::vertex_point_t>
  : CGAL::Tag_true {};
} // namespace CGAL




#endif // CGAL_PROPERTIES_SEAM_MESH_H
