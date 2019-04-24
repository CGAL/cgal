// Copyright (c) 2019 GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Mael Rouxel-Labbé

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/internal/Has_member_id.h>
#include <CGAL/boost/graph/properties.h>

#ifndef CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS
  #error CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS is not defined
#endif

#ifndef CGAL_2D_TRIANGULATION
  #error CGAL_2D_TRIANGULATION is not defined
#endif

// note only the properties below are protected by the macro,
// the rest of the file is the shared implementation of properties for all 2D triangulations
#ifndef CGAL_BOOST_GRAPH_PROPERTIES_2D_TRIANGULATION_H
#define CGAL_BOOST_GRAPH_PROPERTIES_2D_TRIANGULATION_H

namespace CGAL {
namespace detail {

template <class Tr>
struct T2_halfedge_descriptor;

template <class Tr>
struct T2_edge_descriptor;

template <typename Tr>
class T2_vertex_point_map
{
public:
  typedef boost::lvalue_property_map_tag                           category;
  typedef typename Tr::Point                                       value_type;
  typedef value_type&                                              reference;
  typedef typename Tr::Vertex_handle                               key_type;

  T2_vertex_point_map() { }

  friend reference get(T2_vertex_point_map<Tr>, key_type vh)
  {
    return vh->point();
  }
  friend void put(T2_vertex_point_map<Tr>, key_type vh, const value_type& v)
  {
    vh->point() = v;
  }

  reference operator[](key_type vh) const { return vh->point(); }
};

template <typename Tr>
class T2_edge_weight_map
  : public boost::put_get_helper<typename Tr::Geom_traits::FT,
                                 T2_edge_weight_map<Tr> >
{
public:
  typedef boost::readable_property_map_tag                        category;
  typedef typename Tr::Geom_traits::FT                            value_type;
  typedef value_type                                              reference;
  typedef typename Tr::Edge                                       key_type;

  T2_edge_weight_map(const Tr& tr_) : tr(tr_) { }

  value_type operator[](key_type e) const { return approximate_sqrt(tr.segment(e).squared_length()); }

private:
  const Tr& tr;
};

template <typename Tr>
class T2_vertex_id_map
  : public boost::put_get_helper<int, T2_vertex_id_map<Tr> >
{
public:
  typedef boost::readable_property_map_tag                         category;
  typedef int                                                      value_type;
  typedef int                                                      reference;
  typedef typename Tr::Vertex_handle                               key_type;

  T2_vertex_id_map() { }

  value_type operator[](key_type vh) const { return vh->id(); }
};

template <typename Tr>
class T2_halfedge_id_map
  : public boost::put_get_helper<int, T2_halfedge_id_map<Tr> >
{
public:
  typedef boost::readable_property_map_tag                         category;
  typedef int                                                      value_type;
  typedef int                                                      reference;
  typedef CGAL::detail::T2_halfedge_descriptor<Tr>                 key_type;

  T2_halfedge_id_map() { }

  value_type operator[](key_type h) const { return (3 * h.first->id()) + h.second; }
};

template <typename Tr>
class T2_edge_id_map
  : public boost::put_get_helper<int, T2_edge_id_map<Tr> >
{
public:
  typedef boost::readable_property_map_tag                         category;
  typedef int                                                      value_type;
  typedef int                                                      reference;
  typedef CGAL::detail::T2_edge_descriptor<Tr>                     key_type;
  typedef typename Tr::Face_handle                                 face_descriptor;

  T2_edge_id_map() { }

  value_type operator[](key_type e) const
  {
    const face_descriptor f1 = e.first;
    const face_descriptor f2 = f1->neighbor(e.second);

    if(f1->id() < f2->id())
      return 3 * f1->id() + e.second;
    else
      return 3 * f2->id() + f2->index(f1);
  }
};

template <typename Tr>
class T2_face_id_map
  : public boost::put_get_helper<int, T2_face_id_map<Tr> >
{
public:
  typedef boost::readable_property_map_tag                         category;
  typedef int                                                      value_type;
  typedef int                                                      reference;
  typedef typename Tr::Face_handle                                 key_type;

  T2_face_id_map() { }

  value_type operator[](key_type f) const { return f->id(); }
};

template <typename Tr, typename Tag>
struct T2_property_map { };

template <typename Tr>
struct T2_property_map<Tr, boost::vertex_point_t>
{
  typedef detail::T2_vertex_point_map<Tr> type;
  typedef detail::T2_vertex_point_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::edge_weight_t>
{
  typedef detail::T2_edge_weight_map<Tr> type;
  typedef detail::T2_edge_weight_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::vertex_index_t>
{
  typedef detail::T2_vertex_id_map<Tr> type;
  typedef detail::T2_vertex_id_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::halfedge_index_t>
{
  typedef detail::T2_halfedge_id_map<Tr> type;
  typedef detail::T2_halfedge_id_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::edge_index_t>
{
  typedef detail::T2_edge_id_map<Tr> type;
  typedef detail::T2_edge_id_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::face_index_t>
{
  typedef detail::T2_face_id_map<Tr> type;
  typedef detail::T2_face_id_map<Tr> const_type;
};

} // end namespace detail

} // CGAL

#endif // CGAL_BOOST_GRAPH_PROPERTIES_2D_TRIANGULATION_H

// overloads and specializations in the boost namespace
namespace boost {

// g++ 'enumeral_type' in template unification not implemented workaround
template <CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS, class Tag>
struct property_map<CGAL_2D_TRIANGULATION, Tag>
{
  typedef typename CGAL::detail::T2_property_map<CGAL_2D_TRIANGULATION, Tag>    map_gen;
  typedef typename map_gen::type                                                type;
  typedef typename map_gen::const_type                                          const_type;
};

// see struct property_map in Polyehdron for an explanation
template <CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS, class Tag>
struct property_map<const CGAL_2D_TRIANGULATION, Tag>
{
  typedef typename CGAL::detail::T2_property_map<CGAL_2D_TRIANGULATION, Tag>    map_gen;
  typedef typename map_gen::type                                                type;
  typedef typename map_gen::const_type                                          const_type;
};

} // end namespace boost

namespace CGAL {

template <CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct graph_has_property<CGAL_2D_TRIANGULATION, boost::vertex_point_t>
  : CGAL::Tag_true{};
template<CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct graph_has_property<CGAL_2D_TRIANGULATION, boost::edge_weight_t>
  : CGAL::Tag_true{};

template<CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct graph_has_property<CGAL_2D_TRIANGULATION, boost::vertex_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL_2D_TRIANGULATION::Vertex
      >::value
    >
{};
template<CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct graph_has_property<CGAL_2D_TRIANGULATION, boost::halfedge_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL_2D_TRIANGULATION::Face
      >::value
    >
{};
template<CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct graph_has_property<CGAL_2D_TRIANGULATION, boost::edge_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL_2D_TRIANGULATION::Face
      >::value
    >
{};
template<CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
struct graph_has_property<CGAL_2D_TRIANGULATION, boost::face_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL_2D_TRIANGULATION::Face
      >::value
    >
{};

// property maps
template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline detail::T2_vertex_point_map< CGAL_2D_TRIANGULATION >
get(boost::vertex_point_t, const CGAL_2D_TRIANGULATION&)
{
  detail::T2_vertex_point_map< CGAL_2D_TRIANGULATION > m;
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline detail::T2_edge_weight_map< CGAL_2D_TRIANGULATION >
get(boost::edge_weight_t, const CGAL_2D_TRIANGULATION& g)
{
  detail::T2_edge_weight_map< CGAL_2D_TRIANGULATION > m(g);
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline detail::T2_vertex_id_map< CGAL_2D_TRIANGULATION >
get(boost::vertex_index_t, const CGAL_2D_TRIANGULATION&)
{
  detail::T2_vertex_id_map< CGAL_2D_TRIANGULATION > m;
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline detail::T2_halfedge_id_map< CGAL_2D_TRIANGULATION >
get(boost::halfedge_index_t, const CGAL_2D_TRIANGULATION&)
{
  detail::T2_halfedge_id_map< CGAL_2D_TRIANGULATION > m;
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline detail::T2_edge_id_map< CGAL_2D_TRIANGULATION >
get(boost::edge_index_t, const CGAL_2D_TRIANGULATION&)
{
  detail::T2_edge_id_map< CGAL_2D_TRIANGULATION > m;
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline detail::T2_face_id_map< CGAL_2D_TRIANGULATION >
get(boost::face_index_t, const CGAL_2D_TRIANGULATION&)
{
  detail::T2_face_id_map< CGAL_2D_TRIANGULATION > m;
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS, class PropertyTag, class Key>
inline
typename boost::property_traits<
  typename boost::property_map< CGAL_2D_TRIANGULATION, PropertyTag>::const_type>::value_type
get(PropertyTag p, const CGAL_2D_TRIANGULATION& g, const Key& key)
{
  return get(get(p, g), key);
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS, class PropertyTag, class Key, class Value>
inline void
put(PropertyTag p, CGAL_2D_TRIANGULATION& g, const Key& key, const Value& value)
{
  typedef typename boost::property_map<CGAL_2D_TRIANGULATION, PropertyTag>::type Map;
  Map pmap = get(p, g);
  put(pmap, key, value);
}

} // namespace CGAL

#undef CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS
#undef CGAL_2D_TRIANGULATION
