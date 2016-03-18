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



#endif /* CGAL_PROPERTIES_TRIMESH_ARRAYKERNELT_H */
