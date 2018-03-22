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


#ifndef CGAL_PROPERTIES_SURFACE_MESH_H
#define CGAL_PROPERTIES_SURFACE_MESH_H

#ifndef DOXYGEN_RUNNING

#include <CGAL/license/Surface_mesh.h>

#include <CGAL/assertions.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/Properties.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/number_utils.h>

#include <CGAL/boost/graph/properties.h>

#include <boost/cstdint.hpp>

namespace CGAL {

template <typename Point>
class SM_edge_weight_pmap 
  : public boost::put_get_helper<typename CGAL::Kernel_traits<Point>::type::FT, SM_edge_weight_pmap<Point> >
{
  typedef CGAL::Surface_mesh<Point> SM;
public:
  typedef boost::readable_property_map_tag                category;
  typedef typename CGAL::Kernel_traits<Point>::type::FT   value_type;
  typedef value_type                                      reference;
  typedef typename SM::Edge_index                        key_type;

  SM_edge_weight_pmap(const CGAL::Surface_mesh<Point>& sm)
    : pm_(sm. template property_map<
            typename SM::Vertex_index,
            typename SM::Point >("v:point").first),
      sm_(sm)
    {}

  value_type operator[](const key_type& e) const
  {
    return approximate_sqrt(CGAL::squared_distance(pm_[source(e, sm_)],
                                                   pm_[target(e, sm_)]));
  }

private:
   typename SM::template Property_map< typename SM::Vertex_index,
                                       typename SM::Point > pm_;
  const SM& sm_;
};


template <typename K, typename VEF>
class SM_index_pmap : public boost::put_get_helper<boost::uint32_t, SM_index_pmap<K,VEF> >
{
public:
  typedef boost::readable_property_map_tag category;
  typedef boost::uint32_t                  value_type;
  typedef boost::uint32_t                  reference;
  typedef VEF                              key_type;

  value_type operator[](const key_type& vd) const
  {
    return vd;
  }
};

} // CGAL

namespace boost {
//
// edge_weight
//

template <typename Point>
struct property_map<CGAL::Surface_mesh<Point>, boost::edge_weight_t >
{
  typedef CGAL::SM_edge_weight_pmap<Point> type;
  typedef CGAL::SM_edge_weight_pmap<Point> const_type;
};
}
namespace CGAL{
template <typename Point>
typename boost::property_map<CGAL::Surface_mesh<Point>, boost::edge_weight_t>::const_type
get(boost::edge_weight_t, const CGAL::Surface_mesh<Point>& sm)
{
  return CGAL::SM_edge_weight_pmap<Point>(sm);
}

// forward declarations, see <CGAL/Surface_mesh.h>
class SM_Vertex_index;
class SM_Edge_index;
class SM_Halfedge_index;
class SM_Face_index;

template <typename Point>
typename CGAL::Kernel_traits<Point>::type::FT
get(boost::edge_weight_t, const CGAL::Surface_mesh<Point>& sm,
    const SM_Edge_index& e)
{
  return CGAL::SM_edge_weight_pmap<Point>(sm)[e];
}
}
//
// vertex_index
//

namespace boost{
template <typename Point>
struct property_map<CGAL::Surface_mesh<Point>, boost::vertex_index_t >
{
  typedef CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::vertex_descriptor> type;
  typedef CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::vertex_descriptor> const_type;
};
}
namespace CGAL{

template <typename Point>
CGAL::SM_index_pmap<Point, SM_Vertex_index>
get(const boost::vertex_index_t&, const CGAL::Surface_mesh<Point>&)
{
  return CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::vertex_descriptor>();
}
}
//
// face_index
//
namespace boost{
template <typename Point>
struct property_map<CGAL::Surface_mesh<Point>, boost::face_index_t >
{
  typedef CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::face_descriptor> type;
  typedef CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::face_descriptor> const_type;
};
}
namespace CGAL{

template <typename Point>
CGAL::SM_index_pmap<Point, SM_Face_index>
get(const boost::face_index_t&, const CGAL::Surface_mesh<Point>&)
{
  return CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::face_descriptor>();
}
}
//
// edge_index
//
namespace boost{
template <typename Point>
struct property_map<CGAL::Surface_mesh<Point>, boost::edge_index_t >
{
  typedef CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::edge_descriptor> type;
  typedef CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::edge_descriptor> const_type;
};
}
namespace CGAL{
template <typename Point>
CGAL::SM_index_pmap<Point, SM_Edge_index>
get(const boost::edge_index_t&, const CGAL::Surface_mesh<Point>&)
{
  return CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::edge_descriptor>();
}
}
//
// halfedge_index
//
namespace boost{
template <typename Point>
struct property_map<CGAL::Surface_mesh<Point>, boost::halfedge_index_t >
{
  typedef CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::halfedge_descriptor> type;
  typedef CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::halfedge_descriptor> const_type;
};
}
namespace CGAL{

template <typename Point>
CGAL::SM_index_pmap<Point, SM_Halfedge_index>
get(const boost::halfedge_index_t&, const CGAL::Surface_mesh<Point>&)
{
  return CGAL::SM_index_pmap<Point, typename boost::graph_traits<CGAL::Surface_mesh<Point> >::halfedge_descriptor>();
}
}
//
// vertex_point
//
namespace boost{
template<typename P>
struct property_map<CGAL::Surface_mesh<P>, CGAL::vertex_point_t >
{
  typedef CGAL::Surface_mesh<P> SM;

  typedef typename
    SM::template Property_map< typename SM::Vertex_index,
                               P
                               > type;
  
  typedef type const_type;

};
}
namespace CGAL{

namespace internal {
  template <typename Point>
  struct Get_vertex_point_map_for_Surface_mesh_return_type {
    typedef typename boost::property_map
    <
      CGAL::Surface_mesh<Point>,
      CGAL::vertex_point_t
      >::const_type type;
  };
} // end namespace internal

template<typename Point>
typename
boost::lazy_disable_if
<
  boost::is_const<Point>,
  internal::Get_vertex_point_map_for_Surface_mesh_return_type<Point>
>::type
get(CGAL::vertex_point_t, const CGAL::Surface_mesh<Point>& g) {
  return g.points();
}

namespace internal {
  template <typename Point>
  struct Get_graph_traits_of_SM {
    typedef boost::graph_traits< CGAL::Surface_mesh<Point> > type;
  };
} // end namespace internal

// get for intrinsic properties
#define CGAL_SM_INTRINSIC_PROPERTY(RET, PROP, TYPE)                     \
 template<typename Point>                                              \
 RET                                                                   \
 get(PROP p, const CGAL::Surface_mesh<Point>& sm,                      \
     const TYPE& x)                                                \
 { return get(get(p, sm), x); }                                        \

CGAL_SM_INTRINSIC_PROPERTY(boost::uint32_t, boost::vertex_index_t,
SM_Vertex_index)
CGAL_SM_INTRINSIC_PROPERTY(boost::uint32_t, boost::edge_index_t,
SM_Edge_index)
CGAL_SM_INTRINSIC_PROPERTY(boost::uint32_t, boost::halfedge_index_t,
SM_Halfedge_index)
CGAL_SM_INTRINSIC_PROPERTY(boost::uint32_t, boost::face_index_t,
SM_Face_index)
CGAL_SM_INTRINSIC_PROPERTY(Point&, CGAL::vertex_point_t, SM_Vertex_index)

#undef CGAL_SM_INTRINSIC_PROPERTY

// put for intrinsic properties
// only available for vertex_point
template<typename Point>
void
put(CGAL::vertex_point_t p, const CGAL::Surface_mesh<Point>& g,
    typename boost::graph_traits< CGAL::Surface_mesh<Point> >::vertex_descriptor x,
    const Point& point) {
  typedef CGAL::Surface_mesh<Point> SM;
  CGAL_assertion(g.is_valid(x));
  typename SM::template Property_map< typename boost::graph_traits<SM>::vertex_descriptor,
                    Point> prop = get(p, g);
  prop[x] = point;
}

template<typename Point>
struct graph_has_property<CGAL::Surface_mesh<Point>, boost::vertex_index_t>
  : CGAL::Tag_true {};
template<typename Point>
struct graph_has_property<CGAL::Surface_mesh<Point>, boost::edge_index_t>
  : CGAL::Tag_true {};
template<typename Point>
struct graph_has_property<CGAL::Surface_mesh<Point>, boost::halfedge_index_t>
  : CGAL::Tag_true {};
template<typename Point>
struct graph_has_property<CGAL::Surface_mesh<Point>, boost::face_index_t>
  : CGAL::Tag_true {};
template<typename Point>
struct graph_has_property<CGAL::Surface_mesh<Point>, CGAL::vertex_point_t>
  : CGAL::Tag_true {};
template<typename Point>
struct graph_has_property<CGAL::Surface_mesh<Point>, boost::edge_weight_t>
  : CGAL::Tag_true {};
} // CGAL

namespace boost
{



template <typename Point, typename T>
struct property_map<CGAL::Surface_mesh<Point>, CGAL::dynamic_vertex_property_t<T> >
{
  typedef CGAL::Surface_mesh<Point> SM;
  typedef typename SM:: template Property_map<typename SM::Vertex_index,T> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef type const_type;
};

template <typename Point, typename T>
struct property_map<CGAL::Surface_mesh<Point>, CGAL::dynamic_face_property_t<T> >
{
  typedef CGAL::Surface_mesh<Point> SM;
  typedef typename SM:: template Property_map<typename SM::Face_index,T> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef type const_type;
};

template <typename Point, typename T>
struct property_map<CGAL::Surface_mesh<Point>, CGAL::dynamic_halfedge_property_t<T> >
{
  typedef CGAL::Surface_mesh<Point> SM;
  typedef typename SM:: template Property_map<typename SM::Halfedge_index,T> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef type const_type;
};

template <typename Point, typename T>
struct property_map<CGAL::Surface_mesh<Point>, CGAL::dynamic_edge_property_t<T> >
{
  typedef CGAL::Surface_mesh<Point> SM;
  typedef typename SM:: template Property_map<typename SM::Edge_index,T> SMPM;
  typedef CGAL::internal::Dynamic<SM, SMPM> type;
  typedef type const_type;
};

} // nmamespace boost

namespace CGAL {

template <typename Point, typename T>
typename boost::property_map<CGAL::Surface_mesh<Point>, dynamic_vertex_property_t<T> >::const_type
get(dynamic_vertex_property_t<T>, Surface_mesh<Point>& sm)
{
  typedef typename boost::property_map<Surface_mesh<Point>, dynamic_vertex_property_t<T> >::SMPM SMPM;
  typedef typename boost::property_map<Surface_mesh<Point>, dynamic_vertex_property_t<T> >::const_type DPM;
  return DPM(sm, new SMPM(sm.template add_property_map<typename Surface_mesh<Point>::Vertex_index, T>(std::string()).first));
}

template <typename Point, typename T>
typename boost::property_map<Surface_mesh<Point>, dynamic_face_property_t<T> >::const_type
get(dynamic_face_property_t<T>, Surface_mesh<Point>& sm)
{
  typedef typename boost::property_map<Surface_mesh<Point>, dynamic_face_property_t<T> >::SMPM SMPM;
  typedef typename boost::property_map<Surface_mesh<Point>, dynamic_face_property_t<T> >::const_type DPM;
  return DPM(sm, new SMPM(sm.template add_property_map<typename Surface_mesh<Point>::Face_index, T>(std::string()).first));
}

template <typename Point, typename T>
typename boost::property_map<Surface_mesh<Point>, dynamic_edge_property_t<T> >::const_type
get(dynamic_edge_property_t<T>, Surface_mesh<Point>& sm)
{
  typedef typename boost::property_map<Surface_mesh<Point>, dynamic_edge_property_t<T> >::SMPM SMPM;
  typedef typename boost::property_map<Surface_mesh<Point>, dynamic_edge_property_t<T> >::const_type DPM;
  return DPM(sm, new SMPM(sm.template add_property_map<typename Surface_mesh<Point>::Edge_index, T>(std::string()).first));
}

template <typename Point, typename T>
typename boost::property_map<Surface_mesh<Point>, dynamic_halfedge_property_t<T> >::const_type
get(dynamic_halfedge_property_t<T>, Surface_mesh<Point>& sm)
{
  typedef typename boost::property_map<Surface_mesh<Point>, dynamic_halfedge_property_t<T> >::SMPM SMPM;
  typedef typename boost::property_map<Surface_mesh<Point>, dynamic_halfedge_property_t<T> >::const_type DPM;
  return DPM(sm, new SMPM(sm.template add_property_map<typename Surface_mesh<Point>::Halfedge_index, T>(std::string()).first));
}

// implementation detail: required by Dynamic_property_map_deleter
template <typename Pmap,typename P>
void
remove_property(Pmap pm, CGAL::Surface_mesh<P>& sm)
{
  return sm.remove_property_map(pm);
}


template <typename P, typename Property_tag>
struct Get_pmap_of_surface_mesh {
  typedef typename boost::property_map<Surface_mesh<P>, Property_tag >::type type;
};


} // namespace CGAL

#include <CGAL/boost/graph/properties_Surface_mesh_time_stamp.h>
#include <CGAL/boost/graph/properties_Surface_mesh_features.h>

#endif // DOXYGEN_RUNNING

#endif /* CGAL_PROPERTIES_SURFACE_MESH_H */
