// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Philipp Möller


#include <CGAL/assertions.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <boost/mpl/if.hpp>

#ifndef OPEN_MESH_CLASS
  #error OPEN_MESH_CLASS is not defined
#endif

// note only the properties below are protected by the macro,
// the rest of the file is the shared implementation of properties for
// OpenMesh::PolyMesh_ArrayKernelT<K> and OpenMesh::TriMesh_ArrayKernelT<K>
#ifndef CGAL_BOOST_GRAPH_PROPERTIES_OPENMESH_H
#define CGAL_BOOST_GRAPH_PROPERTIES_OPENMESH_H
namespace CGAL {

template <typename Mesh, typename Descriptor, typename Value>
class OM_pmap {
public:
  typedef typename boost::mpl::if_<boost::is_same<Descriptor, typename boost::graph_traits<Mesh>::vertex_descriptor>,
                                   OpenMesh::VPropHandleT<Value>,
                                   typename boost::mpl::if_<boost::is_same<Descriptor, typename boost::graph_traits<Mesh>::face_descriptor>,
                                                            OpenMesh::FPropHandleT<Value>,
                                                            typename boost::mpl::if_<boost::is_same<Descriptor, typename boost::graph_traits<Mesh>::halfedge_descriptor>,
                                                                                     OpenMesh::HPropHandleT<Value>,
                                                                                     OpenMesh::EPropHandleT<Value> >::type>::type>::type H;

  typedef boost::read_write_property_map_tag category;

  typedef Descriptor key_type;
  typedef Value value_type;

  typedef value_type& reference;

  OM_pmap()
  {}

  OM_pmap(Mesh& m)
    : mesh(&m)
  {
    mesh->add_property(h);
  }

  OM_pmap(Mesh& m, H h)
    : mesh(&m), h(h)
  {
    mesh->add_property(h);
  }

  inline friend reference get(const OM_pmap<Mesh,Descriptor,Value>& pm,
                              typename boost::graph_traits<Mesh>::vertex_descriptor k)
  {
    return pm.mesh->property(pm.h,k);
  }

  inline friend reference get(const OM_pmap<Mesh,Descriptor,Value>& pm,
                         typename boost::graph_traits<Mesh>::face_descriptor k)
  {
    return pm.mesh->property(pm.h,k);
  }

  inline friend reference get(const OM_pmap<Mesh,Descriptor,Value>& pm,
                         typename boost::graph_traits<Mesh>::halfedge_descriptor k)
  {
    return pm.mesh->property(pm.h,k);
  }

  inline friend reference get(const OM_pmap<Mesh,Descriptor,Value>& pm,
                         typename boost::graph_traits<Mesh>::edge_descriptor k)
  {
    typename Mesh::EdgeHandle eh(k.idx());
    return pm.mesh->property(pm.h,eh);
  }
  inline friend void put(const OM_pmap<Mesh,Descriptor,Value>& pm,
                         typename boost::graph_traits<Mesh>::vertex_descriptor k,
                         const value_type& v)
  {
    pm.mesh->property(pm.h,k) = v;
  }

  inline friend void put(const OM_pmap<Mesh,Descriptor,Value>& pm,
                         typename boost::graph_traits<Mesh>::face_descriptor k,
                         const value_type& v)
  {
    pm.mesh->property(pm.h,k) = v;
  }

  inline friend void put(const OM_pmap<Mesh,Descriptor,Value>& pm,
                         typename boost::graph_traits<Mesh>::halfedge_descriptor k,
                         const value_type& v)
  {
    pm.mesh->property(pm.h,k) = v;
  }
  inline friend void put(const OM_pmap<Mesh,Descriptor,Value>& pm,
                         typename boost::graph_traits<Mesh>::edge_descriptor k,
                         const value_type& v)
  {
    typename Mesh::EdgeHandle eh(k.idx());
    pm.mesh->property(pm.h,eh) = v;
  }

  reference operator[](key_type k)
  {
    return mesh->property(h,k);
  }

  H handle() const
  {
    return h;
  }

  H& handle()
  {
    return h;
  }

  Mesh* mesh;
  H h;
};


template <typename OpenMesh>
class OM_edge_weight_pmap
  : public boost::put_get_helper<typename OpenMesh::Scalar , OM_edge_weight_pmap<OpenMesh> >
{
public:
  typedef boost::readable_property_map_tag                         category;
  typedef typename OpenMesh::Scalar                                value_type;
  typedef value_type                                               reference;
  typedef typename boost::graph_traits<OpenMesh>::edge_descriptor  key_type;

  OM_edge_weight_pmap(const OpenMesh& sm)
    : sm_(sm)
    {}

  value_type operator[](const key_type& e) const
  {
    return sm_.calc_edge_length(e.halfedge());
  }

private:
  const OpenMesh& sm_;
};

template <typename K, typename VEF>
class OM_index_pmap : public boost::put_get_helper<unsigned int, OM_index_pmap<K,VEF> >
{
public:
  typedef boost::readable_property_map_tag category;
  typedef unsigned int                      value_type;
  typedef unsigned int                      reference;
  typedef VEF                              key_type;

  OM_index_pmap()
  {}

  value_type operator[](const key_type& vd) const
  {
    return vd.idx();
  }
};


template<typename OpenMesh, typename P>
class OM_point_pmap //: public boost::put_get_helper<bool, OM_point_pmap<OpenMesh> >
{
public:
  typedef boost::read_write_property_map_tag category;
#if defined(CGAL_USE_OM_POINTS)
  typedef typename OpenMesh::Point             value_type;
  typedef const typename OpenMesh::Point&      reference;
#else
  typedef P value_type;
  typedef P reference;
#endif
  typedef typename boost::graph_traits<OpenMesh>::vertex_descriptor key_type;

  OM_point_pmap()
    : sm_(NULL)
  {}

  OM_point_pmap(const OpenMesh& sm)
    : sm_(&sm)
    {}

  OM_point_pmap(const OM_point_pmap& pm)
    : sm_(pm.sm_)
    {}

  value_type operator[](key_type v)
  {
#if defined(CGAL_USE_OM_POINTS)
    return sm_->point(v);
#else
    CGAL_assertion(sm_!=NULL);
    typename OpenMesh::Point const& omp = sm_->point(v);
    return value_type(omp[0], omp[1], omp[2]);
#endif
  }

  inline friend reference get(const OM_point_pmap<OpenMesh,P>& pm, key_type v)
  {
    CGAL_precondition(pm.sm_!=NULL);
#if defined(CGAL_USE_OM_POINTS)
    return pm.sm_->point(v);
#else
    CGAL_assertion(pm.sm_!=NULL);
    typename OpenMesh::Point const& omp = pm.sm_->point(v);
    return value_type(omp[0], omp[1], omp[2]);
#endif
  }

  inline friend void put(const OM_point_pmap<OpenMesh,P>& pm, key_type v, const value_type& p)
  {
    CGAL_precondition(pm.sm_!=NULL);
#if defined(CGAL_USE_OM_POINTS)
    const_cast<OpenMesh&>(*pm.sm_).set_point(v,p);
#else
    const_cast<OpenMesh&>(*pm.sm_).set_point
      (v, typename OpenMesh::Point((float)p[0], (float)p[1], (float)p[2]));
#endif
  }

  private:
  const OpenMesh* sm_;
};
} // CGAL
#endif // CGAL_BOOST_GRAPH_PROPERTIES_OPENMESH_H

// overloads and specializations in the boost namespace
namespace boost {

//
// edge_weight
//


template <typename K>
struct property_map<OPEN_MESH_CLASS, boost::edge_weight_t >
{
  typedef OPEN_MESH_CLASS Mesh;
  typedef CGAL::OM_edge_weight_pmap<Mesh> type;
  typedef CGAL::OM_edge_weight_pmap<Mesh> const_type;
};



//
// vertex_index
//

template <typename K>
struct property_map<OPEN_MESH_CLASS, boost::vertex_index_t >
{
  typedef OPEN_MESH_CLASS Mesh;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<Mesh>::vertex_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<Mesh>::vertex_descriptor> const_type;
};


//
// face_index
//

template <typename K>
struct property_map<OPEN_MESH_CLASS, boost::face_index_t >
{
  typedef OPEN_MESH_CLASS Mesh;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<Mesh>::face_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<Mesh>::face_descriptor> const_type;
};

//
// edge_index
//

template <typename K>
struct property_map<OPEN_MESH_CLASS, boost::edge_index_t >
{
  typedef OPEN_MESH_CLASS Mesh;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<Mesh>::edge_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<Mesh>::edge_descriptor> const_type;
};

//
// halfedge_index
//

template <typename K>
struct property_map<OPEN_MESH_CLASS, boost::halfedge_index_t >
{
  typedef OPEN_MESH_CLASS Mesh;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<Mesh>::halfedge_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<Mesh>::halfedge_descriptor> const_type;
};


template<typename K>
struct property_map<OPEN_MESH_CLASS, boost::vertex_point_t >
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 P;
  typedef OPEN_MESH_CLASS Mesh;
  typedef CGAL::OM_point_pmap<Mesh, P> type;
  typedef type const_type;
};
} // namespace boost

namespace CGAL{

template<typename K>
struct graph_has_property<OPEN_MESH_CLASS, boost::edge_weight_t>
  : CGAL::Tag_true{};
template<typename K>
struct graph_has_property<OPEN_MESH_CLASS, boost::vertex_index_t>
  : CGAL::Tag_true{};
template<typename K>
struct graph_has_property<OPEN_MESH_CLASS, boost::face_index_t>
  : CGAL::Tag_true{};
template<typename K>
struct graph_has_property<OPEN_MESH_CLASS, boost::edge_index_t>
  : CGAL::Tag_true{};
template<typename K>
struct graph_has_property<OPEN_MESH_CLASS, boost::halfedge_index_t>
  : CGAL::Tag_true{};
template<typename K>
struct graph_has_property<OPEN_MESH_CLASS, boost::vertex_point_t>
  : CGAL::Tag_true{};
} //end CGAL

namespace OpenMesh {


template <typename K>
typename boost::property_map<OPEN_MESH_CLASS, boost::edge_weight_t>::const_type
get(boost::edge_weight_t, const OPEN_MESH_CLASS& sm)
{
  typedef OPEN_MESH_CLASS Mesh;
  return CGAL::OM_edge_weight_pmap<Mesh>(sm);
}

template <typename K>
typename OPEN_MESH_CLASS::Scalar
get(boost::edge_weight_t, const OPEN_MESH_CLASS& sm,
    const typename boost::graph_traits<OPEN_MESH_CLASS >::edge_descriptor& e)
{
  typedef OPEN_MESH_CLASS Mesh;
  return CGAL::OM_edge_weight_pmap<Mesh>(sm)[e];
}


template <typename K>
CGAL::OM_index_pmap<K, typename boost::graph_traits<OPEN_MESH_CLASS >::vertex_descriptor>
get(const boost::vertex_index_t&, const OPEN_MESH_CLASS&)
{
  typedef OPEN_MESH_CLASS Mesh;
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<Mesh>::vertex_descriptor>();
}

template <typename K>
typename boost::property_map<OPEN_MESH_CLASS, boost::face_index_t>::const_type
get(const boost::face_index_t&, const OPEN_MESH_CLASS&)
{
  typedef OPEN_MESH_CLASS Mesh;
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<Mesh>::face_descriptor>();
}

template <typename K>
CGAL::OM_index_pmap<K, typename boost::graph_traits<OPEN_MESH_CLASS >::edge_descriptor>
get(const boost::edge_index_t&, const OPEN_MESH_CLASS&)
{
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<OPEN_MESH_CLASS >::edge_descriptor>();
}

template <typename K>
CGAL::OM_index_pmap<K, typename boost::graph_traits<OPEN_MESH_CLASS >::halfedge_descriptor>
get(const boost::halfedge_index_t&, const OPEN_MESH_CLASS&)
{
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<OPEN_MESH_CLASS >::halfedge_descriptor>();
}

template<typename K>
CGAL::OM_point_pmap<OPEN_MESH_CLASS,
                    typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3>
get(boost::vertex_point_t, const OPEN_MESH_CLASS& g)
{
  typedef typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 P;
  return CGAL::OM_point_pmap<OPEN_MESH_CLASS, P>(g);
}

// get for intrinsic properties
#define CGAL_OM_INTRINSIC_PROPERTY(RET, PROP, TYPE)                     \
  template<typename K>                                              \
  RET                                                                   \
  get(PROP p, const OPEN_MESH_CLASS& sm,                      \
      typename boost::graph_traits< OPEN_MESH_CLASS >::TYPE x) \
  { return get(get(p, sm), x); }                                        \

  CGAL_OM_INTRINSIC_PROPERTY(int, boost::vertex_index_t, vertex_descriptor)
  CGAL_OM_INTRINSIC_PROPERTY(int, boost::edge_index_t, edge_descriptor)
  CGAL_OM_INTRINSIC_PROPERTY(int, boost::halfedge_index_t, halfedge_descriptor)
  CGAL_OM_INTRINSIC_PROPERTY(int, boost::face_index_t, face_descriptor)
  //  CGAL_OM_INTRINSIC_PROPERTY(std::size_t, boost::halfedge_index_t, face_descriptor)
  CGAL_OM_INTRINSIC_PROPERTY(typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3, boost::vertex_point_t, vertex_descriptor)

#undef CGAL_OM_INTRINSIC_PROPERTY

// put for intrinsic properties
// only available for vertex_point

template<typename K>
void
put(boost::vertex_point_t p, OPEN_MESH_CLASS& g,
    typename boost::graph_traits< OPEN_MESH_CLASS >::vertex_descriptor vd,
    const typename K::Point& point)
{
  put(get(p,g), vd, point);
}


} // namespace OpenMesh



namespace boost {

template <typename K, typename V>
  struct property_map<OPEN_MESH_CLASS, CGAL::dynamic_vertex_property_t<V> >
{
  typedef OPEN_MESH_CLASS SM;
  typedef typename boost::graph_traits<SM>::vertex_descriptor vertex_descriptor;
  typedef CGAL::OM_pmap<SM,vertex_descriptor, V> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef type const_type;
};

template <typename K, typename V>
struct property_map<OPEN_MESH_CLASS, CGAL::dynamic_halfedge_property_t<V> >
{
  typedef OPEN_MESH_CLASS SM;
  typedef typename boost::graph_traits<SM>::halfedge_descriptor halfedge_descriptor;
  typedef CGAL::OM_pmap<SM,halfedge_descriptor, V> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef type const_type;
};

template <typename K, typename V>
struct property_map<OPEN_MESH_CLASS, CGAL::dynamic_edge_property_t<V> >
{
  typedef OPEN_MESH_CLASS SM;
  typedef typename boost::graph_traits<SM>::edge_descriptor edge_descriptor;
  typedef CGAL::OM_pmap<SM,edge_descriptor, V> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef type const_type;
};

template <typename K, typename V>
struct property_map<OPEN_MESH_CLASS, CGAL::dynamic_face_property_t<V> >
{
  typedef OPEN_MESH_CLASS SM;
  typedef typename boost::graph_traits<SM>::face_descriptor face_descriptor;
  typedef CGAL::OM_pmap<SM,face_descriptor, V> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef type const_type;
};

} // namespace boost


namespace OpenMesh {

template <typename K, typename V>
typename boost::property_map<OPEN_MESH_CLASS, CGAL::dynamic_vertex_property_t<V> >::const_type
get(CGAL::dynamic_vertex_property_t<V>, OPEN_MESH_CLASS& om)
{
  typedef OPEN_MESH_CLASS OM;
  typedef typename boost::property_map<OM, CGAL::dynamic_vertex_property_t<V> >::SMPM SMPM;
  typedef typename boost::property_map<OM, CGAL::dynamic_vertex_property_t<V> >::const_type DPM;
  return DPM(om, new SMPM(om));
}

template <typename K, typename V>
typename boost::property_map<OPEN_MESH_CLASS, CGAL::dynamic_halfedge_property_t<V> >::const_type
get(CGAL::dynamic_halfedge_property_t<V>, OPEN_MESH_CLASS& om)
{
  typedef OPEN_MESH_CLASS OM;
  typedef typename boost::property_map<OM, CGAL::dynamic_halfedge_property_t<V> >::SMPM SMPM;
  typedef typename boost::property_map<OM, CGAL::dynamic_halfedge_property_t<V> >::const_type DPM;
  return DPM(om, new SMPM(om));
}

template <typename K, typename V>
typename boost::property_map<OPEN_MESH_CLASS, CGAL::dynamic_edge_property_t<V> >::const_type
get(CGAL::dynamic_edge_property_t<V>, OPEN_MESH_CLASS& om)
{
  typedef OPEN_MESH_CLASS OM;
  typedef typename boost::property_map<OM, CGAL::dynamic_edge_property_t<V> >::SMPM SMPM;
  typedef typename boost::property_map<OM, CGAL::dynamic_edge_property_t<V> >::const_type DPM;
  return DPM(om, new SMPM(om));
}

template <typename K, typename V>
typename boost::property_map<OPEN_MESH_CLASS, CGAL::dynamic_face_property_t<V> >::const_type
get(CGAL::dynamic_face_property_t<V>, OPEN_MESH_CLASS& om)
{
  typedef OPEN_MESH_CLASS OM;
  typedef typename boost::property_map<OM, CGAL::dynamic_face_property_t<V> >::SMPM SMPM;
  typedef typename boost::property_map<OM, CGAL::dynamic_face_property_t<V> >::const_type DPM;
  return DPM(om, new SMPM(om));
}

// implementation detail: required by Dynamic_property_map_deleter
template <typename Pmap, typename K>
void remove_property(Pmap pm, OPEN_MESH_CLASS& om)
{
  om.remove_property(pm.handle());
}

} // namespace OpenMesh

#undef OPEN_MESH_CLASS
