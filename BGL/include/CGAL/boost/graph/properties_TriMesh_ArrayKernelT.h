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
// 
//
// Author(s)     : Philipp Möller


#ifndef CGAL_PROPERTIES_TRIMESH_ARRAYKERNELT_H
#define CGAL_PROPERTIES_TRIMESH_ARRAYKERNELT_H

#include <CGAL/assertions.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <boost/mpl/if.hpp>


// overloads and specializations in the boost namespace
namespace boost {

//
// edge_weight
//
template <typename K>
struct property_map<OpenMesh::TriMesh_ArrayKernelT<K>, boost::edge_weight_t >
{
  typedef OpenMesh::TriMesh_ArrayKernelT<K> Mesh;
  typedef CGAL::OM_edge_weight_pmap<Mesh> type;
  typedef CGAL::OM_edge_weight_pmap<Mesh> const_type;
};



//
// vertex_index
//
template <typename K>
struct property_map<OpenMesh::TriMesh_ArrayKernelT<K>, boost::vertex_index_t >
{
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::vertex_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::vertex_descriptor> const_type;
};


//
// face_index
//

template <typename K>
struct property_map<OpenMesh::TriMesh_ArrayKernelT<K>, boost::face_index_t >
{
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::face_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::face_descriptor> const_type;
};

//
// edge_index
//

template <typename K>
struct property_map<OpenMesh::TriMesh_ArrayKernelT<K>, boost::edge_index_t >
{
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::edge_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::edge_descriptor> const_type;
};

//
// halfedge_index
//

template <typename K>
struct property_map<OpenMesh::TriMesh_ArrayKernelT<K>, boost::halfedge_index_t >
{
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::halfedge_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::halfedge_descriptor> const_type;
};


template<typename K>
struct property_map<OpenMesh::TriMesh_ArrayKernelT<K>, boost::vertex_point_t >
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 P;
  typedef OpenMesh::TriMesh_ArrayKernelT<K> Mesh;
  typedef CGAL::OM_point_pmap<Mesh,P> type;
  typedef type const_type;
};

template<typename K>
struct graph_has_property<OpenMesh::TriMesh_ArrayKernelT<K>, edge_weight_t>
  : CGAL::Tag_true{};
template<typename K>
struct graph_has_property<OpenMesh::TriMesh_ArrayKernelT<K>, vertex_index_t>
  : CGAL::Tag_true{};
template<typename K>
struct graph_has_property<OpenMesh::TriMesh_ArrayKernelT<K>, face_index_t>
  : CGAL::Tag_true{};
template<typename K>
struct graph_has_property<OpenMesh::TriMesh_ArrayKernelT<K>, edge_index_t>
  : CGAL::Tag_true{};
template<typename K>
struct graph_has_property<OpenMesh::TriMesh_ArrayKernelT<K>, halfedge_index_t>
  : CGAL::Tag_true{};
template<typename K>
struct graph_has_property<OpenMesh::TriMesh_ArrayKernelT<K>, vertex_point_t>
  : CGAL::Tag_true{};

} // namespace boost

namespace OpenMesh {


template <typename K>
typename boost::property_map<OpenMesh::TriMesh_ArrayKernelT<K>, boost::edge_weight_t>::const_type
get(boost::edge_weight_t, const OpenMesh::TriMesh_ArrayKernelT<K>& sm)
{
  return CGAL::OM_edge_weight_pmap<K>(sm);
}

template <typename K>
typename OpenMesh::TriMesh_ArrayKernelT<K>::Scalar
get(boost::edge_weight_t, const OpenMesh::TriMesh_ArrayKernelT<K>& sm, 
    const typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::edge_descriptor& e)
{
  return CGAL::OM_edge_weight_pmap<K>(sm)[e];
}


template <typename K>
CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::vertex_descriptor>
get(const boost::vertex_index_t&, const OpenMesh::TriMesh_ArrayKernelT<K>&)
{
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::vertex_descriptor>();
}

template <typename K>
typename boost::property_map<OpenMesh::TriMesh_ArrayKernelT<K>, boost::face_index_t>::const_type
get(const boost::face_index_t&, const OpenMesh::TriMesh_ArrayKernelT<K>&)
{
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::face_descriptor>();
}

template <typename K>
CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::edge_descriptor>
get(const boost::edge_index_t&, const OpenMesh::TriMesh_ArrayKernelT<K>&)
{
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::edge_descriptor>();
}

template <typename K>
CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::halfedge_descriptor>
get(const boost::halfedge_index_t&, const OpenMesh::TriMesh_ArrayKernelT<K>&)
{
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::TriMesh_ArrayKernelT<K> >::halfedge_descriptor>();
}

template<typename K>
CGAL::OM_point_pmap<OpenMesh::TriMesh_ArrayKernelT<K>,
                    typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3>
get(boost::vertex_point_t, const OpenMesh::TriMesh_ArrayKernelT<K>& g) 
{
  typedef typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 P;
  typedef OpenMesh::TriMesh_ArrayKernelT<K> Mesh;
  return CGAL::OM_point_pmap<Mesh, P>(g);
}

// get for intrinsic properties
#define CGAL_OM_INTRINSIC_PROPERTY(RET, PROP, TYPE)                     \
  template<typename K>                                              \
  RET                                                                   \
  get(PROP p, const OpenMesh::TriMesh_ArrayKernelT<K>& sm,                      \
      typename boost::graph_traits< OpenMesh::TriMesh_ArrayKernelT<K> >::TYPE x) \
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
put(boost::vertex_point_t p, OpenMesh::TriMesh_ArrayKernelT<K>& g,
    typename boost::graph_traits< OpenMesh::TriMesh_ArrayKernelT<K> >::vertex_descriptor vd,
    const typename K::Point& point) 
{
  put(get(p,g), vd, point);
}


} // namespace OpenMesh




namespace CGAL {
namespace internal {

template <typename K, typename V>
struct dynamic_property_map<OpenMesh::TriMesh_ArrayKernelT<K>, vertex_property_t<V> >
{
  typedef OpenMesh::TriMesh_ArrayKernelT<K> SM;
  typedef typename boost::graph_traits<SM>::vertex_descriptor vertex_descriptor;
  typedef CGAL::OM_pmap<SM,vertex_descriptor, V> type;
  typedef type const_type;
};

template <typename K, typename V>
struct dynamic_property_map<OpenMesh::TriMesh_ArrayKernelT<K>, edge_property_t<V> >
{
  typedef OpenMesh::TriMesh_ArrayKernelT<K> SM;
  typedef typename boost::graph_traits<SM>::edge_descriptor edge_descriptor;
  typedef CGAL::OM_pmap<SM,edge_descriptor, V> type;
  typedef type const_type;
};

template <typename K, typename V>
struct dynamic_property_map<OpenMesh::TriMesh_ArrayKernelT<K>, halfedge_property_t<V> >
{
  typedef OpenMesh::TriMesh_ArrayKernelT<K> SM;
  typedef typename boost::graph_traits<SM>::halfedge_descriptor halfedge_descriptor;
  typedef CGAL::OM_pmap<SM,halfedge_descriptor, V> type;
  typedef type const_type;
};

template <typename K, typename V>
struct dynamic_property_map<OpenMesh::TriMesh_ArrayKernelT<K>, face_property_t<V> >
{
  typedef OpenMesh::TriMesh_ArrayKernelT<K> SM;
  typedef typename boost::graph_traits<SM>::face_descriptor face_descriptor;
  typedef CGAL::OM_pmap<SM,face_descriptor, V> type;
  typedef type const_type;
};



template <typename K, typename V>
typename dynamic_property_map<OpenMesh::TriMesh_ArrayKernelT<K>, vertex_property_t<V> >::const_type
add_property(vertex_property_t<V>, OpenMesh::TriMesh_ArrayKernelT<K>& om)
{
  typedef OpenMesh::TriMesh_ArrayKernelT<K> OM;
  typedef typename boost::graph_traits<OM>::vertex_descriptor vertex_descriptor;
  return CGAL::OM_pmap<OM,vertex_descriptor, V>(om);
}

template <typename K, typename V>
typename dynamic_property_map<OpenMesh::TriMesh_ArrayKernelT<K>, halfedge_property_t<V> >::const_type
add_property(halfedge_property_t<V>, OpenMesh::TriMesh_ArrayKernelT<K>& om)
{
  typedef OpenMesh::TriMesh_ArrayKernelT<K> OM;
  typedef typename boost::graph_traits<OM>::halfedge_descriptor halfedge_descriptor;
  return CGAL::OM_pmap<OM,halfedge_descriptor, V>(om);
}

template <typename K, typename V>
typename dynamic_property_map<OpenMesh::TriMesh_ArrayKernelT<K>, edge_property_t<V> >::const_type
add_property(edge_property_t<V>, OpenMesh::TriMesh_ArrayKernelT<K>& om)
{
  typedef OpenMesh::TriMesh_ArrayKernelT<K> OM;
  typedef typename boost::graph_traits<OM>::edge_descriptor edge_descriptor;
  return CGAL::OM_pmap<OM,edge_descriptor, V>(om);
}

template <typename K, typename V>
typename dynamic_property_map<OpenMesh::TriMesh_ArrayKernelT<K>, face_property_t<V> >::const_type
add_property(face_property_t<V>, OpenMesh::TriMesh_ArrayKernelT<K>& om)
{
  typedef OpenMesh::TriMesh_ArrayKernelT<K> OM;
  typedef typename boost::graph_traits<OM>::face_descriptor face_descriptor;
  return CGAL::OM_pmap<OM,face_descriptor, V>(om);
}

template <typename Pmap, typename K>
void remove_property(Pmap pm, OpenMesh::TriMesh_ArrayKernelT<K>& om)
{
  om.remove_property(pm.handle());
}
} // namespace internal
} // namespace CGAL

#endif /* CGAL_PROPERTIES_TRIMESH_ARRAYKERNELT_H */
