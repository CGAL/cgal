// Copyright (c) 2019 GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
namespace internal {

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
  typedef CGAL::internal::T2_edge_descriptor<Tr>                  key_type;

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

  T2_vertex_id_map(const Tr& tr) : tr(tr) { }

  value_type operator[](key_type v) const {
    CGAL_precondition(!tr.is_infinite(v));
    return v->id();
  }

  const Tr& tr;
};

template <typename Tr>
class T2_halfedge_id_map
  : public boost::put_get_helper<int, T2_halfedge_id_map<Tr> >
{
public:
  typedef boost::readable_property_map_tag                         category;
  typedef int                                                      value_type;
  typedef int                                                      reference;
  typedef CGAL::internal::T2_halfedge_descriptor<Tr>               key_type;

  typedef typename Tr::Face_handle                                 Face_handle;

  T2_halfedge_id_map(const Tr& tr) : tr(tr) { }

  // Halfedge id is twice the edge id, and +0/+1 depending whether
  // h.first is such that h.first < opposite(h).first
  value_type operator[](key_type h) const
  {
    const Face_handle f1 = h.first;
    const Face_handle f2 = f1->neighbor(h.second);
    CGAL_assertion(!tr.is_infinite(f1) || !tr.is_infinite(f2));

    if(tr.is_infinite(f1))
      return 2*(f2->edge_id(f2->index(f1)));
    else if(tr.is_infinite(f2))
      return 2*(f1->edge_id(h.second)) + 1;
    else if(f1->id() < f2->id())
      return 2*(f1->edge_id(h.second));
    else
      return 2*(f1->edge_id(h.second)) + 1;
  }

private:
  const Tr& tr;
};

template <typename Tr>
class T2_edge_id_map
  : public boost::put_get_helper<int, T2_edge_id_map<Tr> >
{
public:
  typedef boost::readable_property_map_tag                         category;
  typedef int                                                      value_type;
  typedef int                                                      reference;
  typedef CGAL::internal::T2_edge_descriptor<Tr>                   key_type;
  typedef typename Tr::Face_handle                                 Face_handle;

  T2_edge_id_map(const Tr& tr) : tr(tr) { }

  value_type operator[](key_type e) const
  {
    const Face_handle f1 = e.first;
    const Face_handle f2 = f1->neighbor(e.second);
    CGAL_assertion(!tr.is_infinite(f1) || !tr.is_infinite(f2));

    if(tr.is_infinite(f1))
      return f2->edge_id(f2->index(f1));
    else
      return f1->edge_id(e.second);
  }

private:
  const Tr& tr;
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

  T2_face_id_map(const Tr& tr) : tr(tr) { }

  value_type operator[](key_type f) const {
    CGAL_precondition(!tr.is_infinite(f));
    return f->id();
  }

private:
  const Tr& tr;
};

template <typename Tr, typename Tag>
struct T2_property_map { };

template <typename Tr>
struct T2_property_map<Tr, boost::vertex_point_t>
{
  typedef internal::T2_vertex_point_map<Tr> type;
  typedef internal::T2_vertex_point_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::edge_weight_t>
{
  typedef internal::T2_edge_weight_map<Tr> type;
  typedef internal::T2_edge_weight_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::vertex_index_t>
{
  typedef internal::T2_vertex_id_map<Tr> type;
  typedef internal::T2_vertex_id_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::halfedge_index_t>
{
  typedef internal::T2_halfedge_id_map<Tr> type;
  typedef internal::T2_halfedge_id_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::edge_index_t>
{
  typedef internal::T2_edge_id_map<Tr> type;
  typedef internal::T2_edge_id_map<Tr> const_type;
};

template <typename Tr>
struct T2_property_map<Tr, boost::face_index_t>
{
  typedef internal::T2_face_id_map<Tr> type;
  typedef internal::T2_face_id_map<Tr> const_type;
};

} // end namespace internal
} // CGAL

#endif // CGAL_BOOST_GRAPH_PROPERTIES_2D_TRIANGULATION_H

// overloads and specializations in the boost namespace
namespace boost {

#define CGAL_PM_SPECIALIZATION(TAG) \
template <CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS> \
struct property_map<CGAL_2D_TRIANGULATION, TAG> \
{ \
  typedef typename CGAL::internal::T2_property_map<CGAL_2D_TRIANGULATION, TAG>  map_gen; \
  typedef typename map_gen::type                                                type; \
  typedef typename map_gen::const_type                                          const_type; \
}; \
\
template <CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS> \
struct property_map<const CGAL_2D_TRIANGULATION, TAG> \
{ \
  typedef typename CGAL::internal::T2_property_map<CGAL_2D_TRIANGULATION, TAG>  map_gen; \
  typedef typename map_gen::type                                                type; \
  typedef typename map_gen::const_type                                          const_type; \
};

CGAL_PM_SPECIALIZATION(vertex_point_t)
CGAL_PM_SPECIALIZATION(edge_weight_t)
CGAL_PM_SPECIALIZATION(vertex_index_t)
CGAL_PM_SPECIALIZATION(halfedge_index_t)
CGAL_PM_SPECIALIZATION(edge_index_t)
CGAL_PM_SPECIALIZATION(face_index_t)

#undef CGAL_PM_SPECIALIZATION

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
inline internal::T2_vertex_point_map< CGAL_2D_TRIANGULATION >
get(boost::vertex_point_t, const CGAL_2D_TRIANGULATION&)
{
  internal::T2_vertex_point_map< CGAL_2D_TRIANGULATION > m;
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline internal::T2_edge_weight_map< CGAL_2D_TRIANGULATION >
get(boost::edge_weight_t, const CGAL_2D_TRIANGULATION& g)
{
  internal::T2_edge_weight_map< CGAL_2D_TRIANGULATION > m(g);
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline internal::T2_vertex_id_map< CGAL_2D_TRIANGULATION >
get(boost::vertex_index_t, const CGAL_2D_TRIANGULATION& g)
{
  internal::T2_vertex_id_map< CGAL_2D_TRIANGULATION > m(g);
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline internal::T2_halfedge_id_map< CGAL_2D_TRIANGULATION >
get(boost::halfedge_index_t, const CGAL_2D_TRIANGULATION& g)
{
  internal::T2_halfedge_id_map< CGAL_2D_TRIANGULATION > m(g);
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline internal::T2_edge_id_map< CGAL_2D_TRIANGULATION >
get(boost::edge_index_t, const CGAL_2D_TRIANGULATION& g)
{
  internal::T2_edge_id_map< CGAL_2D_TRIANGULATION > m(g);
  return m;
}

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
inline internal::T2_face_id_map< CGAL_2D_TRIANGULATION >
get(boost::face_index_t, const CGAL_2D_TRIANGULATION& g)
{
  internal::T2_face_id_map< CGAL_2D_TRIANGULATION > m(g);
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

template < CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS >
void set_triangulation_ids(CGAL_2D_TRIANGULATION& g)
{
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::edge_descriptor     edge_descriptor;
  typedef typename boost::graph_traits< CGAL_2D_TRIANGULATION >::face_descriptor     face_descriptor;

  int vid = 0;
  for(vertex_descriptor vd : vertices(g))
    vd->id() = vid++;

  int eid = 0;
  for(edge_descriptor ed : edges(g))
  {
    halfedge_descriptor hd = halfedge(ed, g);
    face_descriptor fd = face(hd, g);
    if(fd != boost::graph_traits< CGAL_2D_TRIANGULATION >::null_face())
      fd->edge_id(hd.second) = eid;

    halfedge_descriptor opp_hd = opposite(hd, g);
    face_descriptor opp_fd = face(opp_hd, g);
    if(opp_fd != boost::graph_traits< CGAL_2D_TRIANGULATION >::null_face())
      opp_fd->edge_id(opp_hd.second) = eid;

    ++eid;
  }

  int fid = 0;
  for(face_descriptor fd : faces(g))
    fd->id() = fid++;
}

} // namespace CGAL

#undef CGAL_2D_TRIANGULATION_TEMPLATE_PARAMETERS
#undef CGAL_2D_TRIANGULATION
