
// Copyright (c) 2019  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
// 
//
// Author(s)     : Andreas Fabri


#ifndef CGAL_PROPERTIES_SURFACE_MESH_H
#define CGAL_PROPERTIES_SURFACE_MESH_H

#ifndef DOXYGEN_RUNNING

#include <pmp/SurfaceMesh.h>

#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/cstdint.hpp>

namespace CGAL {


class pmp_edge_weight_pmap 
 : public boost::put_get_helper<double, pmp_edge_weight_pmap >
{
  typedef pmp::SurfaceMesh SM;
public:
  typedef boost::readable_property_map_tag                category;
  typedef double                                          value_type;
  typedef value_type                                      reference;
  typedef typename boost::graph_traits<pmp::SurfaceMesh>::edge_descriptor key_type;

  pmp_edge_weight_pmap(const pmp::SurfaceMesh& sm)
    : sm(sm)
    {}

  value_type operator[](const key_type& e) const
  {
    return pmp::distance(sm.position(source(e,sm)),
                         sm.position(target(e,sm)));
  }

private:
  
  const SM& sm;
};

  

template <typename VEF>
class pmp_index_pmap : public boost::put_get_helper<boost::uint32_t, pmp_index_pmap<VEF> >
{
public:
  typedef boost::readable_property_map_tag category;
  typedef boost::uint32_t                  value_type;
  typedef boost::uint32_t                  reference;
  typedef VEF                              key_type;

  value_type operator[](const key_type& vd) const
  {
    return vd.idx();
  }
};

} // CGAL


namespace boost {

//
// edge_weight
//

template <>
struct property_map<pmp::SurfaceMesh, boost::edge_weight_t >
{
  typedef CGAL::pmp_edge_weight_pmap type;
  typedef CGAL::pmp_edge_weight_pmap const_type;
};
} // namespace boost

namespace pmp {

typename boost::property_map<SurfaceMesh, boost::edge_weight_t>::const_type
get(boost::edge_weight_t, const SurfaceMesh& sm)
{
  return CGAL::pmp_edge_weight_pmap(sm);
}


Scalar
get(boost::edge_weight_t, const SurfaceMesh& sm,
    const boost::graph_traits<SurfaceMesh>::edge_descriptor& e)
{
  return CGAL::pmp_edge_weight_pmap(sm)[e];
}
  
} // namespace pmp

//
// vertex_index
//

namespace boost {
template <>  
struct property_map<pmp::SurfaceMesh, boost::vertex_index_t >
{
  typedef CGAL::pmp_index_pmap<pmp::Vertex> type;
  typedef CGAL::pmp_index_pmap<pmp::Vertex> const_type;
};
}


namespace pmp{

CGAL::pmp_index_pmap<Vertex>
get(const boost::vertex_index_t&, const SurfaceMesh&)
{
  return CGAL::pmp_index_pmap<Vertex>();
}
  
}

//
// face_index
//
namespace boost{
 
template<>
struct property_map<pmp::SurfaceMesh, boost::face_index_t >
{
  typedef CGAL::pmp_index_pmap<pmp::Face> type;
  typedef CGAL::pmp_index_pmap<pmp::Face> const_type;
};
  
}

namespace pmp{

CGAL::pmp_index_pmap<Face>
get(const boost::face_index_t&, const SurfaceMesh&)
{
  return CGAL::pmp_index_pmap<pmp::Face>();
}
}

//
// edge_index
//
namespace boost{
  
template <>
struct property_map<pmp::SurfaceMesh, boost::edge_index_t >
{
  typedef CGAL::pmp_index_pmap<CGAL::internal::PMP_edge> type;
  typedef CGAL::pmp_index_pmap<CGAL::internal::PMP_edge> const_type;
};
  
}

namespace boost{

template <typename T>
struct property_traits<pmp::VertexProperty<T>> {
  typedef pmp::Vertex key_type;
  typedef T value_type;
  typedef boost::lvalue_property_map_tag category;
  typedef typename pmp::VertexProperty<T>::reference reference;
};
 
template <typename T>
struct property_traits<pmp::HalfedgeProperty<T>> {
  typedef pmp::Halfedge key_type;
  typedef T value_type;
  typedef boost::lvalue_property_map_tag category;
  typedef typename pmp::HalfedgeProperty<T>::reference reference;
};
  
template <typename T>
struct property_traits<pmp::EdgeProperty<T>> {
  typedef pmp::Edge key_type;
  typedef T value_type;
  typedef boost::lvalue_property_map_tag category;
  typedef typename pmp::EdgeProperty<T>::reference reference;
};
  
template <typename T>
struct property_traits<pmp::FaceProperty<T>> {
  typedef pmp::Face key_type;
  typedef T value_type;
  typedef boost::lvalue_property_map_tag category;
  typedef typename pmp::FaceProperty<T>::reference reference;
}; 
}


namespace pmp{

CGAL::pmp_index_pmap<CGAL::internal::PMP_edge>
get(const boost::edge_index_t&, const SurfaceMesh&)
{
  return CGAL::pmp_index_pmap<CGAL::internal::PMP_edge>();
}
}

//
// halfedge_index
//
namespace boost{
  
template <>
struct property_map<pmp::SurfaceMesh, boost::halfedge_index_t >
{
  typedef CGAL::pmp_index_pmap<pmp::Halfedge> type;
  typedef CGAL::pmp_index_pmap<pmp::Halfedge> const_type;
};
  
}

namespace pmp{

CGAL::pmp_index_pmap<Halfedge>
get(const boost::halfedge_index_t&, const SurfaceMesh&)
{
  return CGAL::pmp_index_pmap<Halfedge>();
}
  
}

//
// vertex_point
//
namespace pmp {

template <typename P>
struct VertexPointMap {
  typedef Vertex key_type;
  typedef P value_type;
  typedef P reference;
  typedef boost::read_write_property_map_tag category;

  VertexPointMap(const SurfaceMesh& sm)
    : sm(&const_cast<SurfaceMesh&>(sm))
  {}
  
  inline friend value_type get(const VertexPointMap& pm, key_type v)
  {
    const Point& p = pm.sm->position(v);
    return P(p[0],p[1],p[2]);   
  }
  
  inline friend void put(const VertexPointMap& pm, key_type v, const value_type& p)
  {
    pm.sm->position(v) = Point(p.x(),p.y(),p.z());
  }
  
  SurfaceMesh* sm;
};
  
}


namespace boost{
  
template <>
struct property_map<pmp::SurfaceMesh, CGAL::vertex_point_t >
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 P;
  typedef pmp::SurfaceMesh SM;
  typedef pmp::VertexPointMap<P> type;
  typedef type const_type;
};
  
}


namespace pmp{

boost::property_map<SurfaceMesh, CGAL::vertex_point_t >::type
get(CGAL::vertex_point_t, const SurfaceMesh& g) {
  return boost::property_map<SurfaceMesh, CGAL::vertex_point_t >::type(g);
}


// get for intrinsic properties
#define CGAL_SM_INTRINSIC_PROPERTY(RET, PROP, TYPE)                     \
 RET                                                                   \
 get(PROP p, const SurfaceMesh& sm,                      \
     const TYPE& x)                                                \
 { return get(get(p, sm), x); }                                        \

CGAL_SM_INTRINSIC_PROPERTY(boost::uint32_t, boost::vertex_index_t, Vertex)

CGAL_SM_INTRINSIC_PROPERTY(boost::uint32_t, boost::edge_index_t, CGAL::internal::PMP_edge)

CGAL_SM_INTRINSIC_PROPERTY(boost::uint32_t, boost::halfedge_index_t, Halfedge)

CGAL_SM_INTRINSIC_PROPERTY(boost::uint32_t, boost::face_index_t, Face)

CGAL_SM_INTRINSIC_PROPERTY(CGAL::Exact_predicates_inexact_constructions_kernel::Point_3, CGAL::vertex_point_t, Vertex)

#undef CGAL_SM_INTRINSIC_PROPERTY

// put for intrinsic properties
// only available for vertex_point


void
put(CGAL::vertex_point_t p, SurfaceMesh& sm,
    typename boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
    const Point& point)
{
  sm.position(v) = point;
}

} // namespace pmp

#if 0
namespace boost {

template <>
struct edge_property_type<pmp::SurfaceMesh>
{
  typedef boost::edge_weight_t type;
};
 
}
#endif

namespace CGAL {

template <>
struct graph_has_property<pmp::SurfaceMesh, boost::vertex_index_t>
  : CGAL::Tag_true {};
  
template <>
struct graph_has_property<pmp::SurfaceMesh, boost::edge_index_t>
  : CGAL::Tag_true {};
  
template <>
struct graph_has_property<pmp::SurfaceMesh, boost::halfedge_index_t>
  : CGAL::Tag_true {};

  template <>
struct graph_has_property<pmp::SurfaceMesh, boost::face_index_t>
  : CGAL::Tag_true {};

template <>
struct graph_has_property<pmp::SurfaceMesh, CGAL::vertex_point_t>
  : CGAL::Tag_true {};
  
template <>
struct graph_has_property<pmp::SurfaceMesh, boost::edge_weight_t>
  : CGAL::Tag_true {};

} // namespace CGAL


// dynamic properties
namespace boost {

template <typename T>
struct property_map<pmp::SurfaceMesh, CGAL::dynamic_vertex_property_t<T> >
{
  typedef pmp::SurfaceMesh SM;
  typedef pmp::VertexProperty<T> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef CGAL::internal::Dynamic_with_index<typename pmp::Vertex, T> const_type;
};

template <typename T>
struct property_map<pmp::SurfaceMesh, CGAL::dynamic_face_property_t<T> >
{
  typedef pmp::SurfaceMesh SM;
  typedef pmp::FaceProperty<T> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef CGAL::internal::Dynamic_with_index<typename pmp::Face, T> const_type;
};

template <typename T>
struct property_map<pmp::SurfaceMesh, CGAL::dynamic_halfedge_property_t<T> >
{
  typedef pmp::SurfaceMesh SM;
  typedef pmp::HalfedgeProperty<T> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef CGAL::internal::Dynamic_with_index<typename pmp::Halfedge, T> const_type;
};

template <typename T>
struct property_map<pmp::SurfaceMesh, CGAL::dynamic_edge_property_t<T> >
{
  typedef pmp::SurfaceMesh SM;
  typedef pmp::EdgeProperty<T> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef CGAL::internal::Dynamic_with_index<typename pmp::Edge, T> const_type;
};

} // namespace boost

namespace CGAL {

// get functions for dynamic properties of mutable SurfaceMesh
template <typename T>
typename boost::property_map<pmp::SurfaceMesh, dynamic_vertex_property_t<T> >::type
get(dynamic_vertex_property_t<T>, pmp::SurfaceMesh& sm)
{
  typedef typename boost::property_map<pmp::SurfaceMesh, dynamic_vertex_property_t<T> >::SMPM SMPM;
  typedef typename boost::property_map<pmp::SurfaceMesh, dynamic_vertex_property_t<T> >::type DPM;
  return DPM(sm, new SMPM(sm.template add_vertex_property<T>(std::string())));
}

template <typename T>
typename boost::property_map<pmp::SurfaceMesh, dynamic_face_property_t<T> >::type
get(dynamic_face_property_t<T>, pmp::SurfaceMesh& sm)
{
  typedef typename boost::property_map<pmp::SurfaceMesh, dynamic_face_property_t<T> >::SMPM SMPM;
  typedef typename boost::property_map<pmp::SurfaceMesh, dynamic_face_property_t<T> >::type DPM;
  return DPM(sm, new SMPM(sm.template add_face_property<T>(std::string())));
}

template <typename T>
typename boost::property_map<pmp::SurfaceMesh, dynamic_edge_property_t<T> >::type
get(dynamic_edge_property_t<T>, pmp::SurfaceMesh& sm)
{
  typedef typename boost::property_map<pmp::SurfaceMesh, dynamic_edge_property_t<T> >::SMPM SMPM;
  typedef typename boost::property_map<pmp::SurfaceMesh, dynamic_edge_property_t<T> >::type DPM;
  return DPM(sm, new SMPM(sm.template add_edge_property<T>(std::string())));
}

template <typename T>
typename boost::property_map<pmp::SurfaceMesh, dynamic_halfedge_property_t<T> >::type
get(dynamic_halfedge_property_t<T>, pmp::SurfaceMesh& sm)
{
  typedef typename boost::property_map<pmp::SurfaceMesh, dynamic_halfedge_property_t<T> >::SMPM SMPM;
  typedef typename boost::property_map<pmp::SurfaceMesh, dynamic_halfedge_property_t<T> >::type DPM;
  return DPM(sm, new SMPM(sm.template add_halfedge_property<T>(std::string())));
}

// get functions for dynamic properties of const SurfaceMesh
template <typename T>
typename boost::property_map<pmp::SurfaceMesh, dynamic_vertex_property_t<T> >::const_type
get(dynamic_vertex_property_t<T>, const pmp::SurfaceMesh& sm)
{
  return CGAL::internal::Dynamic_with_index<typename pmp::Vertex, T>(num_vertices(sm));
}

template <typename T>
typename boost::property_map<pmp::SurfaceMesh, dynamic_face_property_t<T> >::const_type
get(dynamic_face_property_t<T>, const pmp::SurfaceMesh& sm)
{
  return CGAL::internal::Dynamic_with_index<typename pmp::Face, T>(num_faces(sm));
}

template <typename T>
typename boost::property_map<pmp::SurfaceMesh, dynamic_halfedge_property_t<T> >::const_type
get(dynamic_halfedge_property_t<T>, const pmp::SurfaceMesh& sm)
{
  return CGAL::internal::Dynamic_with_index<typename pmp::Halfedge, T>(num_halfedges(sm));
}

template <typename Point, typename T>
typename boost::property_map<pmp::SurfaceMesh, dynamic_edge_property_t<T> >::const_type
get(dynamic_edge_property_t<T>, const pmp::SurfaceMesh& sm)
{
  return CGAL::internal::Dynamic_with_index<typename pmp::Edge, T>(num_edges(sm));
}

// implementation detail: required by Dynamic_property_map_deleter
template <typename T>
void
remove_property(pmp::VertexProperty<T> pm, pmp::SurfaceMesh& sm)
{
  return sm.remove_vertex_property(pm);
}

template <typename T>
void
remove_property(pmp::EdgeProperty<T> pm, pmp::SurfaceMesh& sm)
{
  return sm.remove_edge_property(pm);
}
  
template <typename T>
void
remove_property(pmp::HalfedgeProperty<T> pm, pmp::SurfaceMesh& sm)
{
  return sm.remove_halfedge_property(pm);
}

template <typename T>
void
remove_property(pmp::FaceProperty<T> pm, pmp::SurfaceMesh& sm)
{
  return sm.remove_face_property(pm);
}

template <typename Property_tag>
struct Get_pmap_of_surface_mesh {
  typedef typename boost::property_map<pmp::SurfaceMesh, Property_tag >::type type;
};


} // namespace CGAL


#endif // DOXYGEN_RUNNING

#endif /* CGAL_PROPERTIES_SURFACE_MESH_H */
