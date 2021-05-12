// Copyright (c) 2019  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PROPERTIES_TRIANGULATION_DATA_STRUCTURE_2_H
#define CGAL_PROPERTIES_TRIANGULATION_DATA_STRUCTURE_2_H

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/boost/graph/internal/graph_traits_2D_triangulation_helper.h>
#include <CGAL/boost/graph/internal/Has_member_id.h>

#include <CGAL/boost/graph/Named_function_parameters.h>

#include <boost/graph/properties.hpp>

namespace CGAL {
namespace internal {

// property maps
template <class VB, class FB>
class TDS2_vertex_point_map
{
public:
  typedef boost::lvalue_property_map_tag                                      category;
  typedef typename VB::Point                                                  value_type;
  typedef value_type&                                                         reference;
  typedef typename CGAL::Triangulation_data_structure_2<VB,FB>::Vertex_handle key_type;

  friend reference get(TDS2_vertex_point_map<VB,FB>, key_type vh) { return vh->point(); }
  friend void put(TDS2_vertex_point_map<VB,FB>, key_type vh, reference v) { vh->point() = v; }
  reference operator[](key_type vh) const { return vh->point(); }
};

template <class VB, class FB>
class TDS2_edge_weight_map
  : public boost::put_get_helper<typename VB::FT, TDS2_edge_weight_map<VB, FB> >
{
public:
  typedef boost::readable_property_map_tag                           category;
  typedef typename VB::FT                                            value_type;
  typedef value_type                                                 reference;
  typedef typename CGAL::Triangulation_data_structure_2<VB,FB>::Edge key_type;

  TDS2_edge_weight_map(const CGAL::Triangulation_data_structure_2<VB,FB>& tds_) : tds(tds_) { }

  value_type operator[](key_type e) const { return approximate_sqrt(tds.segment(e).squared_length()); }

private:
  const CGAL::Triangulation_data_structure_2<VB,FB>& tds;
};

template <class VB, class FB>
class TDS2_vertex_id_map
  : public boost::put_get_helper<int, TDS2_vertex_id_map<VB, FB> >
{
public:
  typedef boost::readable_property_map_tag                                    category;
  typedef int                                                                 value_type;
  typedef int                                                                 reference;
  typedef typename CGAL::Triangulation_data_structure_2<VB,FB>::Vertex_handle key_type;

  TDS2_vertex_id_map() {}

  long operator[](key_type vh) const { return vh->id(); }
};

template <class VB, class FB>
class TDS2_halfedge_id_map
  : public boost::put_get_helper<int, TDS2_halfedge_id_map<VB, FB> >
{
  typedef typename CGAL::Triangulation_data_structure_2<VB,FB>     TDS;

public:
  typedef boost::readable_property_map_tag                         category;
  typedef int                                                      value_type;
  typedef int                                                      reference;
  typedef CGAL::internal::T2_halfedge_descriptor<TDS>              key_type;

  typedef typename TDS::Face_handle                                face_descriptor;

  TDS2_halfedge_id_map() { }

  // Halfedge id is twice the edge id, and +0/+1 depending whether
  // h.first is such that h.first < opposite(h).first --> different ids
  value_type operator[](key_type h) const
  {
    const face_descriptor f1 = h.first;
    const face_descriptor f2 = f1->neighbor(h.second);

    if(f1->id() < f2->id())
      return 2*(3 * f1->id() + h.second);
    else
      return 2*(3 * f2->id() + f2->index(f1)) + 1;
  }
};

template <class VB, class FB>
class TDS2_edge_id_map
  : public boost::put_get_helper<int, TDS2_edge_id_map<VB, FB> >
{
  typedef typename CGAL::Triangulation_data_structure_2<VB,FB>       TDS;

public:
  typedef boost::readable_property_map_tag                           category;
  typedef int                                                        value_type;
  typedef int                                                        reference;
  typedef CGAL::internal::T2_edge_descriptor<TDS>                    key_type;

  typedef typename TDS::Face_handle                                  Face_handle;

  TDS2_edge_id_map() {}

  value_type operator[](key_type h) const
  {
    const Face_handle f1 = h.first;
    const Face_handle f2 = f1->neighbor(h.second);

    if(f1->id() < f2->id())
      return 3 * f1->id() + h.second;
    else
      return 3 * f2->id() + f2->index(f1);
  }
};

template <class VB, class FB>
class TDS2_face_id_map
  : public boost::put_get_helper<int, TDS2_face_id_map<VB, FB> >
{
  typedef typename CGAL::Triangulation_data_structure_2<VB,FB>     TDS;

public:
  typedef boost::readable_property_map_tag                         category;
  typedef int                                                      value_type;
  typedef int                                                      reference;
  typedef typename TDS::Face_handle                                key_type;

  TDS2_face_id_map() { }

  value_type operator[](key_type f) const { return f->id(); }
};

template <class VB, class FB, class Tag>
struct TDS2_property_map { };

template <class VB, class FB>
struct TDS2_property_map<VB, FB, boost::vertex_point_t>
{
  typedef internal::TDS2_vertex_point_map<VB,FB> type;
  typedef internal::TDS2_vertex_point_map<VB,FB> const_type;
};

template <class VB, class FB>
struct TDS2_property_map<VB, FB, boost::edge_weight_t>
{
  typedef internal::TDS2_edge_weight_map<VB,FB> type;
  typedef internal::TDS2_edge_weight_map<VB,FB> const_type;
};

template <class VB, class FB>
struct TDS2_property_map<VB, FB, boost::vertex_index_t>
{
  typedef internal::TDS2_vertex_id_map<VB,FB> type;
  typedef internal::TDS2_vertex_id_map<VB,FB> const_type;
};

template <class VB, class FB>
struct TDS2_property_map<VB, FB, boost::halfedge_index_t>
{
  typedef internal::TDS2_vertex_id_map<VB,FB> type;
  typedef internal::TDS2_vertex_id_map<VB,FB> const_type;
};

template <class VB, class FB>
struct TDS2_property_map<VB, FB, boost::edge_index_t>
{
  typedef internal::TDS2_edge_id_map<VB,FB> type;
  typedef internal::TDS2_edge_id_map<VB,FB> const_type;
};

template <class VB, class FB>
struct TDS2_property_map<VB, FB, boost::face_index_t>
{
  typedef internal::TDS2_vertex_id_map<VB,FB> type;
  typedef internal::TDS2_vertex_id_map<VB,FB> const_type;
};

} // end namespace internal

template <class VB, class FB >
struct graph_has_property<CGAL::Triangulation_data_structure_2<VB, FB>, boost::vertex_point_t>
  : CGAL::Tag_true{};
template<class VB, class FB >
struct graph_has_property<CGAL::Triangulation_data_structure_2<VB, FB>, boost::edge_weight_t>
  : CGAL::Tag_true{};

template<class VB, class FB >
struct graph_has_property<CGAL::Triangulation_data_structure_2<VB, FB>, boost::vertex_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL::Triangulation_data_structure_2<VB, FB>::Vertex
      >::value
    >
{};
template<class VB, class FB >
struct graph_has_property<CGAL::Triangulation_data_structure_2<VB, FB>, boost::halfedge_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL::Triangulation_data_structure_2<VB, FB>::Face
      >::value
    >
{};
template<class VB, class FB >
struct graph_has_property<CGAL::Triangulation_data_structure_2<VB, FB>, boost::edge_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL::Triangulation_data_structure_2<VB, FB>::Face
      >::value
    >
{};
template<class VB, class FB >
struct graph_has_property<CGAL::Triangulation_data_structure_2<VB, FB>, boost::face_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL::Triangulation_data_structure_2<VB, FB>::Face
      >::value
    >
{};

template <class VB, class FB>
inline internal::TDS2_vertex_point_map<VB,FB>
get(boost::vertex_point_t, const Triangulation_data_structure_2<VB,FB>&)
{
  internal::TDS2_vertex_point_map<VB,FB> m;
  return m;
}

template <class VB, class FB>
inline internal::TDS2_edge_weight_map<VB,FB>
get(boost::edge_weight_t, const Triangulation_data_structure_2<VB,FB>& g)
{
  internal::TDS2_edge_weight_map<VB,FB> m(g);
  return m;
}

template <class VB, class FB>
inline internal::TDS2_vertex_id_map<VB,FB>
get(boost::vertex_index_t, const Triangulation_data_structure_2<VB,FB>&)
{
  internal::TDS2_vertex_id_map<VB,FB> m;
  return m;
}

template <class VB, class FB>
inline internal::TDS2_halfedge_id_map<VB,FB>
get(boost::halfedge_index_t, const Triangulation_data_structure_2<VB,FB>&)
{
  internal::TDS2_halfedge_id_map<VB,FB> m;
  return m;
}

template <class VB, class FB>
inline internal::TDS2_edge_id_map<VB,FB>
get(boost::edge_index_t, const Triangulation_data_structure_2<VB,FB>&)
{
  internal::TDS2_edge_id_map<VB,FB> m;
  return m;
}

template <class VB, class FB>
inline internal::TDS2_face_id_map<VB,FB>
get(boost::face_index_t, const Triangulation_data_structure_2<VB,FB>&)
{
  internal::TDS2_face_id_map<VB,FB> m;
  return m;
}

} // namespace CGAL

namespace boost {

#define CGAL_PM_SPECIALIZATION(TAG) \
template <class VB, class FB> \
struct property_map<CGAL::Triangulation_data_structure_2<VB,FB>, TAG> \
{ \
  typedef typename CGAL::internal::TDS2_property_map<VB, FB, TAG> map_gen; \
  typedef typename map_gen::type type; \
  typedef typename map_gen::const_type const_type; \
}; \
\
template <class VB, class FB> \
struct property_map<const CGAL::Triangulation_data_structure_2<VB,FB>, TAG> \
{ \
  typedef typename CGAL::internal::TDS2_property_map<VB, FB, TAG> map_gen; \
  typedef typename map_gen::type type; \
  typedef typename map_gen::const_type const_type; \
};

CGAL_PM_SPECIALIZATION(vertex_point_t)
CGAL_PM_SPECIALIZATION(edge_weight_t)
CGAL_PM_SPECIALIZATION(vertex_index_t)
CGAL_PM_SPECIALIZATION(halfedge_index_t)
CGAL_PM_SPECIALIZATION(edge_index_t)
CGAL_PM_SPECIALIZATION(face_index_t)

#undef CGAL_PM_SPECIALIZATION

} // namespace boost

namespace CGAL {

template <class VB, class FB, class PropertyTag, class Key>
inline
typename boost::property_traits<
typename boost::property_map<Triangulation_data_structure_2<VB,FB>,PropertyTag>::const_type>::value_type
get(PropertyTag p, const Triangulation_data_structure_2<VB,FB>& g, const Key& key)
{
  return get(get(p, g), key);
}

template <class VB, class FB, class PropertyTag, class Key,class Value>
inline void
put(PropertyTag p, Triangulation_data_structure_2<VB,FB>& g,
    const Key& key, const Value& value)
{
  typedef typename boost::property_map<Triangulation_data_structure_2<VB,FB>, PropertyTag>::type Map;
  Map pmap = get(p, g);
  put(pmap, key, value);
}

} // namespace CGAL

namespace boost {

// What are those needed for ???
template <typename VB, typename FB>
struct edge_property_type<CGAL::Triangulation_data_structure_2<VB,FB> > {
  typedef void type;
};

template <typename VB, typename FB>
struct vertex_property_type<CGAL::Triangulation_data_structure_2<VB,FB> > {
  typedef void type;
};

} // namespace boost

#endif /* CGAL_PROPERTIES_TRIANGULATION_DATA_STRUCTURE_2_H */
