// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Fernando Cacciola


#ifndef CGAL_BOOST_GRAPH_BGL_PROPERTIES_H
#define CGAL_BOOST_GRAPH_BGL_PROPERTIES_H

#include <CGAL/property_map.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/basic.h>

#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>

#include <string>
#include <vector>
#include <type_traits>

namespace CGAL {

template<typename Graph, typename PropertyTag>
struct graph_has_property : CGAL::Tag_false { };

} // namespace CGAL

namespace boost {

enum vertex_point_t          { vertex_point          };

// vertex_index_t is defined in boost
enum vertex_external_index_t { vertex_external_index };

enum halfedge_index_t          { halfedge_index };
enum halfedge_external_index_t { halfedge_external_index };

// edge_index_t is defined in boost
enum edge_external_index_t   { edge_external_index   };

enum face_index_t            { face_index            };
enum face_external_index_t   { face_external_index   };

struct cgal_no_property
{
  typedef bool type;
  typedef const bool const_type;
};

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
using boost::vertex_index_t;
using boost::vertex_index;
using boost::vertex_external_index_t;
using boost::vertex_external_index;
using boost::halfedge_index_t;
using boost::halfedge_index;
using boost::halfedge_external_index_t;
using boost::halfedge_external_index;
using boost::edge_index_t;
using boost::edge_index;
using boost::edge_external_index_t;
using boost::edge_external_index;
using boost::face_index_t;
using boost::face_index;
using boost::face_external_index_t;
using boost::face_external_index;
} // CGAL

namespace CGAL {
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

template<typename Handle, typename ValueType, typename Reference,
         bool is_const = std::is_const<
                           typename std::remove_reference<Reference>::type >::value>
struct Point_accessor
  : boost::put_get_helper< Reference, Point_accessor<Handle, ValueType, Reference> >
{
  typedef boost::lvalue_property_map_tag category;
  typedef Reference                      reference;
  typedef ValueType                      value_type;
  typedef Handle                         key_type;

  reference operator[](Handle h) const { return h->point(); }
};

// partial specialization for const map to make them constructible from non-const map
template<typename Handle, typename ValueType, typename ConstReference>
struct Point_accessor<Handle, ValueType, ConstReference, true>
  : boost::put_get_helper< ConstReference, Point_accessor<Handle, ValueType, ConstReference, true> >
{
  typedef boost::lvalue_property_map_tag category;
  typedef ConstReference                      reference;
  typedef ValueType                      value_type;
  typedef Handle                         key_type;

  typedef typename boost::mpl::if_< boost::is_reference<ConstReference>,
                                    ValueType&,
                                    ValueType >::type Reference;

  Point_accessor() {}
  Point_accessor(Point_accessor<Handle, ValueType, Reference, false>) {}

  reference operator[](Handle h) const { return h->point(); }
};

// this one is basically 'readable_property_map_tag'
template <typename PropertyMap,
          typename PropertyMapCategory = typename boost::property_traits<PropertyMap>::category>
struct Is_writable_property_map : CGAL::Tag_false { };

template <typename PropertyMap>
struct Is_writable_property_map<PropertyMap, boost::writable_property_map_tag> : CGAL::Tag_true { };

template <typename PropertyMap>
struct Is_writable_property_map<PropertyMap, boost::read_write_property_map_tag> : CGAL::Tag_true { };

// 'lvalue_pmap_tag' is annoying, because the property map is allowed to be non-mutable,
// but boost::lvalue_property_map_tag is defined as:
//   struct lvalue_property_map_tag : public read_write_property_map_tag
// so we can't just check that 'writable_property_map_tag' is a base of the the lvalue tag.
//
// This checks if the reference is non-const, which is not completely correct: map[key] returning
// a non-const reference doesn't mean that 'put(map, key, val)' exists, which is what a writable
// property map must define.
template <typename PropertyMap>
struct Is_writable_property_map<PropertyMap, boost::lvalue_property_map_tag>
  : boost::mpl::if_c<std::is_const<typename std::remove_reference<
                       typename boost::property_traits<PropertyMap>::reference>::type>::value,
                     CGAL::Tag_false, CGAL::Tag_true>::type
{ };

} // namespace internal

// Needed by PMP::detec_features and Mesh_3
enum vertex_feature_degree_t    { vertex_feature_degree };
enum edge_is_feature_t          { edge_is_feature };

enum vertex_time_stamp_t        { vertex_time_stamp};
enum halfedge_time_stamp_t      { halfedge_time_stamp};
enum face_time_stamp_t          { face_time_stamp};

template<typename ID>
struct vertex_incident_patches_t {
  typedef ID type;
};

template<typename ID>
struct face_patch_id_t {
  typedef ID type;
};

} // namespace CGAL


#endif // CGAL_BOOST_GRAPH_BGL_PROPERTIES_H
