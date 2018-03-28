// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Andreas Fabri, Fernando Cacciola


#ifndef CGAL_BOOST_GRAPH_BGL_PROPERTIES_H
#define CGAL_BOOST_GRAPH_BGL_PROPERTIES_H

#include <CGAL/property_map.h>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <CGAL/Dynamic_property_map.h>

#include <CGAL/basic.h>
#include <string>

namespace CGAL{
/// \ingroup PkgBGLProperties
/// \brief graph_has_property is used to indicate if
/// a model of `HalfedgeGraph` or `FaceGraph`
/// has an internal property associated with the
/// given `PropertyTag`.
///
/// It inherits from `CGAL::Tag_true` if there is a
/// default internal property map for the
/// corresponding property tag and from
/// `CGAL::Tag_false` otherwise.
///
/// \tparam Graph a model of `HalfedgeGraph` or `FaceGraph`
/// \tparam PropertyTag the type of a property tag
/// referring to the property of interest.
///
template<typename Graph, typename PropertyTag>
struct graph_has_property
#ifndef DOXYGEN_RUNNING
    : CGAL::Tag_false
#endif
{};
}
/// Boost Namespace
namespace boost {

/// \ingroup PkgBGLProperties
/// @{

/// A property tag which refers to the geometric embedding property
/// of a vertex of a \ref HalfedgeGraph.
enum vertex_point_t          { vertex_point          };
enum vertex_external_index_t { vertex_external_index } ;

/// A property tag which refers to the property
/// of a halfedge of being a border halfedge.
enum edge_external_index_t   { edge_external_index   } ;

/// A property tag which identifies the *index* property of
/// a halfedge of a \ref HalfedgeGraph.
enum halfedge_index_t        { halfedge_index        };
enum halfedge_external_index_t   { halfedge_external_index   } ;

/// A property tag which identifies the *index* property of
/// a face of a \ref FaceGraph.
enum face_index_t            { face_index            };
enum face_external_index_t   { face_external_index   } ;

  
struct cgal_no_property
{
  typedef bool type;
  typedef const bool const_type;
};

/// @}

// Introduce those two tags so we can use BOOST_INSTALL_PROPERTY
// macro. This is dangerous because we now rely on implementation
// details.
struct halfedge_property_tag { };
struct face_property_tag { };

BOOST_INSTALL_PROPERTY(vertex, point);
BOOST_INSTALL_PROPERTY(vertex, external_index);
BOOST_INSTALL_PROPERTY(halfedge, external_index);
BOOST_INSTALL_PROPERTY(edge, external_index);
BOOST_INSTALL_PROPERTY(face, index);
BOOST_INSTALL_PROPERTY(face, external_index);
} // boost

namespace CGAL {
using boost::vertex_point_t;
using boost::vertex_point;
using boost::vertex_external_index_t;
using boost::vertex_external_index;
using boost::halfedge_index_t;
using boost::halfedge_index;
using boost::halfedge_external_index_t;
using boost::halfedge_external_index;
using boost::edge_external_index_t;
using boost::edge_external_index;
using boost::face_index_t;
using boost::face_index;
using boost::face_external_index_t;
using boost::face_external_index;
} // CGAL

namespace CGAL{
namespace helpers {

// matches read-write property maps
template <class PolygonMesh, class FaceIndexMap, class Tag>
void init_face_indices(PolygonMesh& pm,
                       FaceIndexMap& fid,
                       boost::read_write_property_map_tag,
                       Tag)
{
  typename boost::property_traits<FaceIndexMap>::value_type i = 0;
  BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::face_descriptor fd,
                faces(pm))
  {
    put(fid, fd, i);
    ++i;
  }
}
template <class PolygonMesh, class VertexIndexMap, class Tag>
void init_vertex_indices(PolygonMesh& pm,
                         VertexIndexMap& vid,
                         boost::read_write_property_map_tag,
                         Tag)
{
  typename boost::property_traits<VertexIndexMap>::value_type i = 0;
  BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd,
                vertices(pm))
  {
    put(vid, vd, i);
    ++i;
  }
}
template <class PolygonMesh, class HalfedgeIndexMap, class Tag>
void init_halfedge_indices(PolygonMesh& pm,
                           HalfedgeIndexMap& hid,
                           boost::read_write_property_map_tag,
                           Tag)
{
  typename boost::property_traits<HalfedgeIndexMap>::value_type i = 0;
  BOOST_FOREACH(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor hd,
                halfedges(pm))
  {
    put(hid, hd, i);
    ++i;
  }
}

// matches mutable Lvalue property maps
template <class PolygonMesh, class FaceIndexMap>
void init_face_indices(PolygonMesh& pm,
                       FaceIndexMap& fid,
                       boost::lvalue_property_map_tag,
                       boost::false_type)
{
  init_face_indices(pm, fid,
    boost::read_write_property_map_tag(), boost::false_type());
}
template <class PolygonMesh, class VertexIndexMap>
void init_vertex_indices(PolygonMesh& pm,
                         VertexIndexMap& vid,
                         boost::lvalue_property_map_tag,
                         boost::false_type)
{
  init_vertex_indices(pm, vid,
    boost::read_write_property_map_tag(), boost::false_type());
}
template <class PolygonMesh, class HalfedgeIndexMap>
void init_halfedge_indices(PolygonMesh& pm,
                         HalfedgeIndexMap& hid,
                         boost::lvalue_property_map_tag,
                         boost::false_type)
{
  init_halfedge_indices(pm, hid,
    boost::read_write_property_map_tag(), boost::false_type());
}

// matches all other types of property map
template <class PolygonMesh, class FaceIndexMap, class MapTag, class Tag>
void init_face_indices(PolygonMesh&, FaceIndexMap, MapTag, Tag)
{}
template <class PolygonMesh, class VertexIndexMap, class MapTag, class Tag>
void init_vertex_indices(PolygonMesh&, VertexIndexMap, MapTag, Tag)
{}
template <class PolygonMesh, class HalfedgeIndexMap, class MapTag, class Tag>
void init_halfedge_indices(PolygonMesh&, HalfedgeIndexMap, MapTag, Tag)
{}

template <class PolygonMesh, class FaceIndexMap>
void init_face_indices(PolygonMesh& pm, FaceIndexMap fid)
{
  init_face_indices(pm, fid,
                    typename boost::property_traits<FaceIndexMap>::category(),
                    typename boost::is_const<
                      typename boost::remove_reference<
                        typename boost::property_traits<FaceIndexMap>::reference
                            >::type >::type() );
}

template <class PolygonMesh, class VertexIndexMap>
void init_vertex_indices(PolygonMesh& pm, VertexIndexMap vid)
{
  init_vertex_indices(pm, vid,
                      typename boost::property_traits<VertexIndexMap>::category(),
                      typename boost::is_const<
                        typename boost::remove_reference<
                          typename boost::property_traits<VertexIndexMap>::reference
                            >::type >::type() );
}

template <class PolygonMesh, class HalfedgeIndexMap>
void init_halfedge_indices(PolygonMesh& pm, HalfedgeIndexMap hid)
{
  init_halfedge_indices(pm, hid,
                        typename boost::property_traits<HalfedgeIndexMap>::category(),
                        typename boost::is_const<
                          typename boost::remove_reference<
                            typename boost::property_traits<HalfedgeIndexMap>::reference
                              >::type >::type() );
}

} //namespace helpers

namespace internal {
  
  template<typename Polyhedron, typename Handle>
struct Index_accessor
    : boost::put_get_helper< std::size_t&, Index_accessor<Polyhedron,Handle> >
{
  typedef boost::lvalue_property_map_tag category;
  typedef std::size_t&                   reference;
  typedef std::size_t                    value_type;
  typedef Handle                         key_type;

  reference operator[](Handle h) const { return h->id(); }
};

template<typename Handle>
struct Edge_index_accessor
  : boost::put_get_helper< std::size_t, Edge_index_accessor<Handle> >
{
  typedef boost::readable_property_map_tag category;
  typedef std::size_t                      reference;
  typedef std::size_t                      value_type;
  typedef Handle                           key_type;

  reference operator[](Handle h) const { return h.id(); }
};

template<typename Handle, typename ValueType, typename Reference>
struct Point_accessor
  : boost::put_get_helper< Reference, Point_accessor<Handle, ValueType, Reference> >
{
  typedef boost::lvalue_property_map_tag category;
  typedef Reference                      reference;
  typedef ValueType                      value_type;
  typedef Handle                         key_type;

  reference operator[](Handle h) const { return h->point(); }
};

} // namespace internal

// Needed by PMP::detec_features and Mesh_3
enum vertex_feature_degree_t    { vertex_feature_degree };
enum edge_is_feature_t          { edge_is_feature };

enum vertex_time_stamp_t        { vertex_time_stamp};
enum halfedge_time_stamp_t      { halfedge_time_stamp};
enum face_time_stamp_t          { face_time_stamp};

template <typename ID>
struct vertex_incident_patches_t {
  typedef ID type;
};

template <typename ID>
struct face_patch_id_t {
  typedef ID type;
};

} // namespace CGAL


#endif // CGAL_BOOST_GRAPH_BGL_PROPERTIES_H
