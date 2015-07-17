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


#ifndef CGAL_PROPERTIES_POLYMESH_ARRAYKERNELT_H
#define CGAL_PROPERTIES_POLYMESH_ARRAYKERNELT_H

#include <CGAL/assertions.h>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <boost/mpl/if.hpp>

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
  
  OM_pmap(Mesh& m)
    : mesh(m)
  {}
  
  OM_pmap(Mesh& m, H h)
    : mesh(m), h(h)
  {}
  
  inline friend reference get(const OM_pmap<Mesh,Descriptor,Value>& pm, key_type k)
  {
    return pm.mesh.property(pm.h,k);
  }

  inline friend void put(const OM_pmap<Mesh,Descriptor,Value>& pm, key_type k, const value_type& v)
  {
    pm.mesh.property(pm.h,k) = v;
  }

  reference operator[](key_type k)
  {
    return mesh.property(h,k);
  }

  H handle() const
  {
    return h;
  }

  H& handle()
  {
    return h;
  }

  Mesh& mesh;
  H h;
};


template <typename K>
class OM_edge_weight_pmap 
  : public boost::put_get_helper<typename OpenMesh::PolyMesh_ArrayKernelT<K>::Scalar , OM_edge_weight_pmap<K> >
{
public:
  typedef boost::readable_property_map_tag                         category;
  typedef typename OpenMesh::PolyMesh_ArrayKernelT<K>::Scalar      value_type;
  typedef value_type                                               reference;
  typedef typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor key_type;

  OM_edge_weight_pmap(const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
    : sm_(sm)
    {}

  value_type operator[](const key_type& e) const
  {
    return sm_.calc_edge_length(e.halfedge());
  }

private:
  const OpenMesh::PolyMesh_ArrayKernelT<K>& sm_;
};

template <typename K, typename VEF>
class OM_index_pmap : public boost::put_get_helper<unsigned int, OM_index_pmap<K,VEF> >
{
public:
  typedef boost::readable_property_map_tag category;
  typedef unsigned int                      value_type;
  typedef unsigned int                      reference;
  typedef VEF                              key_type;

  value_type operator[](const key_type& vd) const
  {
    return vd.idx();
  }
};


  template<typename K, typename P>
class OM_point_pmap //: public boost::put_get_helper<bool, OM_point_pmap<K> >
{
public:
  typedef boost::read_write_property_map_tag category;
#if defined(CGAL_USE_OM_POINTS)
  typedef typename OpenMesh::PolyMesh_ArrayKernelT<K>::Point             value_type;
  typedef const typename OpenMesh::PolyMesh_ArrayKernelT<K>::Point&      reference;
#else
  typedef P value_type;
  typedef P reference;
#endif
  typedef typename boost::graph_traits< OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor key_type;
    
  OM_point_pmap()
    : sm_(NULL)
  {}

  OM_point_pmap(const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
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
    typename OpenMesh::PolyMesh_ArrayKernelT<K>::Point const& omp = sm_->point(v);
    return value_type(omp[0], omp[1], omp[2]);
#endif
  }

  inline friend reference get(const OM_point_pmap<K,P>& pm, key_type v)
  {
    CGAL_precondition(pm.sm_!=NULL);
#if defined(CGAL_USE_OM_POINTS)
    return pm.sm_->point(v);
#else
    CGAL_assertion(pm.sm_!=NULL);
    typename OpenMesh::PolyMesh_ArrayKernelT<K>::Point const& omp = pm.sm_->point(v);
    return value_type(omp[0], omp[1], omp[2]);
#endif
  }

  inline friend void put(const OM_point_pmap<K,P>& pm, key_type v, const value_type& p)
  {
    CGAL_precondition(pm.sm_!=NULL);
#if defined(CGAL_USE_OM_POINTS)
    const_cast<OpenMesh::PolyMesh_ArrayKernelT<K>&>(*pm.sm_).set_point(v,p);
#else
    const_cast<OpenMesh::PolyMesh_ArrayKernelT<K>&>(*pm.sm_).set_point
      (v, typename OpenMesh::PolyMesh_ArrayKernelT<K>::Point((float)p[0], (float)p[1], (float)p[2]));
#endif
  }

  private:
  const OpenMesh::PolyMesh_ArrayKernelT<K>* sm_;
};


} // CGAL

// overloads and specializations in the boost namespace
namespace boost {

//
// edge_weight
//


template <typename K>
struct property_map<OpenMesh::PolyMesh_ArrayKernelT<K>, boost::edge_weight_t >
{
  typedef CGAL::OM_edge_weight_pmap<K> type;
  typedef CGAL::OM_edge_weight_pmap<K> const_type;
};



//
// vertex_index
//

template <typename K>
struct property_map<OpenMesh::PolyMesh_ArrayKernelT<K>, boost::vertex_index_t >
{
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor> const_type;
};


//
// face_index
//

template <typename K>
struct property_map<OpenMesh::PolyMesh_ArrayKernelT<K>, boost::face_index_t >
{
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor> const_type;
};

//
// edge_index
//

template <typename K>
struct property_map<OpenMesh::PolyMesh_ArrayKernelT<K>, boost::edge_index_t >
{
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor> const_type;
};

//
// halfedge_index
//

template <typename K>
struct property_map<OpenMesh::PolyMesh_ArrayKernelT<K>, boost::halfedge_index_t >
{
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor> type;
  typedef CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor> const_type;
};


template<typename K>
struct property_map<OpenMesh::PolyMesh_ArrayKernelT<K>, boost::vertex_point_t >
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 P;
  typedef CGAL::OM_point_pmap<K,P> type;
  typedef type const_type;
};

} // namespace boost

namespace boost {


template <typename K>
typename boost::property_map<OpenMesh::PolyMesh_ArrayKernelT<K>, boost::edge_weight_t>::const_type
get(boost::edge_weight_t, const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return CGAL::OM_edge_weight_pmap<K>(sm);
}

template <typename K>
typename OpenMesh::PolyMesh_ArrayKernelT<K>::Scalar
get(boost::edge_weight_t, const OpenMesh::PolyMesh_ArrayKernelT<K>& sm, 
    const typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor& e)
{
  return CGAL::OM_edge_weight_pmap<K>(sm)[e];
}


template <typename K>
CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor>
get(const boost::vertex_index_t&, const OpenMesh::PolyMesh_ArrayKernelT<K>&)
{
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor>();
}

template <typename K>
typename boost::property_map<OpenMesh::PolyMesh_ArrayKernelT<K>, boost::face_index_t>::const_type
get(const boost::face_index_t&, const OpenMesh::PolyMesh_ArrayKernelT<K>&)
{
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor>();
}

template <typename K>
CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor>
get(const boost::edge_index_t&, const OpenMesh::PolyMesh_ArrayKernelT<K>&)
{
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor>();
}

template <typename K>
CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor>
get(const boost::halfedge_index_t&, const OpenMesh::PolyMesh_ArrayKernelT<K>&)
{
  return CGAL::OM_index_pmap<K, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor>();
}

template<typename K>
CGAL::OM_point_pmap<K,typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3>
get(boost::vertex_point_t, const OpenMesh::PolyMesh_ArrayKernelT<K>& g) 
{
  typedef typename CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 P;
  return CGAL::OM_point_pmap<K,P>(g);
}

// get for intrinsic properties
#define CGAL_OM_INTRINSIC_PROPERTY(RET, PROP, TYPE)                     \
  template<typename K>                                              \
  RET                                                                   \
  get(PROP p, const OpenMesh::PolyMesh_ArrayKernelT<K>& sm,                      \
      typename boost::graph_traits< OpenMesh::PolyMesh_ArrayKernelT<K> >::TYPE x) \
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
put(boost::vertex_point_t p, OpenMesh::PolyMesh_ArrayKernelT<K>& g,
    typename boost::graph_traits< OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor vd,
    const typename K::Point& point) 
{
  put(get(p,g), vd, point);
}


} // namespace boost



#endif /* CGAL_PROPERTIES_POLYMESH_ARRAYKERNELT_H */
