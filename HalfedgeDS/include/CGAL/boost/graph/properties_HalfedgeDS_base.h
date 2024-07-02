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

// note only the properties below are protected by the macro,
// the rest of the file is the shared implementation of properties for
// Polyhedron and HalfedgeDS_default
#ifndef CGAL_BOOST_GRAPH_PROPERTIES_HALFEDGEDS_BASE_H
#define CGAL_BOOST_GRAPH_PROPERTIES_HALFEDGEDS_BASE_H

#include <CGAL/boost/graph/properties.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/number_utils.h>
#include <memory>
#include <CGAL/boost/graph/internal/Has_member_id.h>
#include <CGAL/Distance_3/Point_3_Point_3.h>
#include <CGAL/Dynamic_property_map.h>

namespace CGAL {

namespace internal {

template<class Handle>
class HDS_index_map_external
  : public boost::put_get_helper<std::size_t&, HDS_index_map_external<Handle> >
{
public:
  typedef boost::lvalue_property_map_tag   category;
  typedef std::size_t                      value_type;
  typedef std::size_t&                     reference;
  typedef Handle                           key_type;

private:
  typedef CGAL::Unique_hash_map<key_type,std::size_t> Map;

public:
  template <typename InputIterator>
  HDS_index_map_external(InputIterator begin, InputIterator end, std::size_t max)
    : map_(new Map(begin, end, 0, std::size_t(-1), max)) {}

  reference operator[](const key_type& k) const { return (*map_)[k]; }
private:
   std::shared_ptr<Map> map_;
};

// Special case for edges.
template<class Polyhedron>
class HDS_edge_index_map_external
  : public boost::put_get_helper<std::size_t&, HDS_edge_index_map_external<Polyhedron> >
{
public:
  typedef boost::lvalue_property_map_tag                            category;
  typedef std::size_t                                               value_type;
  typedef std::size_t&                                              reference;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor key_type;

private:
  typedef CGAL::Unique_hash_map<key_type,std::size_t> Map;

public:
  HDS_edge_index_map_external(Polyhedron& p)
    : map_(new Map(std::size_t(-1), num_halfedges(p)))
  {
    unsigned int data = 0;
    typename boost::graph_traits<Polyhedron>::edge_iterator it, end;
    for(boost::tie(it, end) = edges(p); it != end; ++it, ++data)
      (*map_)[*it] = data;
  }

  reference operator[](const key_type& k) const { return (*map_)[k]; }
private:
  std::shared_ptr<Map> map_;
};

template<typename Handle, typename FT>
struct HDS_wrap_squared
{
  typedef FT value_type;
  typedef FT reference;
  typedef Handle key_type;
  typedef boost::readable_property_map_tag category;

  template<typename E>
  FT operator[](const E& e) const {
    return approximate_sqrt(CGAL::squared_distance(e.halfedge()->vertex()->point(),
                                                   e.halfedge()->opposite()->vertex()->point()));
  }

  friend inline
  value_type get(const HDS_wrap_squared& m, const key_type k)
  {
    return m[k];
  }
};

}

// the tag we dispatch on from property_map<G, Property>
template <class HDS, class Tag>
struct HDS_property_map {};

} // end of CGAL::internal namespace

#endif // CGAL_BOOST_GRAPH_PROPERTIES_HALFEDGEDS_BASE_H

#if !defined(CGAL_HDS_TMPLT) || ! defined(CGAL_HDS_CLASS)
#error CGAL_HDS_TMPLT or CGAL_HDS_CLASS is not defined
#endif

namespace CGAL {

// generalized 2-ary get functions
template<class CGAL_HDS_TMPLT, class PropertyTag>
typename boost::property_map<CGAL_HDS_CLASS, PropertyTag >::const_type
get(PropertyTag,CGAL_HDS_CLASS const&)
{ return typename boost::property_map<CGAL_HDS_CLASS, PropertyTag >::const_type(); }

template<class CGAL_HDS_TMPLT, class PropertyTag>
typename boost::property_map<CGAL_HDS_CLASS, PropertyTag >::type
get(PropertyTag,CGAL_HDS_CLASS&)
{ return typename boost::property_map<CGAL_HDS_CLASS, PropertyTag >::type(); }


// generalized 3-ary get functions
template<class CGAL_HDS_TMPLT, class PropertyTag, class Key,
         class F = std::enable_if_t<!std::is_same_v<PropertyTag, dynamic_vertex_property_t<Key>> &&
                                    !std::is_same_v<PropertyTag, dynamic_halfedge_property_t<Key>> &&
                                    !std::is_same_v<PropertyTag, dynamic_edge_property_t<Key>> &&
                                    !std::is_same_v<PropertyTag, dynamic_face_property_t<Key>>
                                   >
>
typename boost::property_traits< typename boost::property_map<CGAL_HDS_CLASS, PropertyTag >::type >::reference
get(PropertyTag p,CGAL_HDS_CLASS& g, const Key& key)
{
  return get(get(p, g), key);
}

template<class CGAL_HDS_TMPLT, class PropertyTag, class Key>
typename boost::property_traits< typename boost::property_map<CGAL_HDS_CLASS, PropertyTag >::const_type >::reference
get(PropertyTag p,CGAL_HDS_CLASS const& g, const Key& key)
{ return get(get(p, g), key); }



#define DECLARE_HDS_DYNAMIC_PM(TAG, DESCRIPTOR)                                  \
template <typename CGAL_HDS_TMPLT, class T>                                      \
typename boost::property_map<CGAL_HDS_CLASS, TAG >::const_type                   \
get(const TAG&, const CGAL_HDS_CLASS&, const T& dv = T())                        \
{                                                                                \
  typedef typename boost::graph_traits< CGAL_HDS_CLASS >::DESCRIPTOR descriptor; \
  return internal::Dynamic_property_map<descriptor,T>(dv);                       \
}

DECLARE_HDS_DYNAMIC_PM(dynamic_vertex_property_t<T>, vertex_descriptor)
DECLARE_HDS_DYNAMIC_PM(dynamic_halfedge_property_t<T>, halfedge_descriptor)
DECLARE_HDS_DYNAMIC_PM(dynamic_edge_property_t<T>, edge_descriptor)
DECLARE_HDS_DYNAMIC_PM(dynamic_face_property_t<T>, face_descriptor)

#undef DECLARE_HDS_DYNAMIC_PM

// generalized put
template<class CGAL_HDS_TMPLT, class PropertyTag, class Key,class Value>
void put(PropertyTag p,CGAL_HDS_CLASS& g, const Key& key, const Value& value)
{
  typedef typename boost::property_map<CGAL_HDS_CLASS, PropertyTag>::type Map;
  Map pmap = get(p, g);
  put(pmap, key, value);
}

// specialization needs to be repeated for halfedge, vertex, face
#define DECLARE_HDS_INDEX_PM(ENTITY, TAG, ACCESSOR)               \
  template<class CGAL_HDS_TMPLT>                                  \
  struct HDS_property_map<CGAL_HDS_CLASS,                         \
                          boost::ENTITY##TAG> {                   \
    struct bind_ {                                                \
      typedef internal::ACCESSOR##_accessor<                      \
        CGAL_HDS_CLASS,                                           \
        typename boost::graph_traits< CGAL_HDS_CLASS              \
                                    >::ENTITY##_descriptor > type;\
      typedef type const_type;                                    \
    };                                                            \
  };

DECLARE_HDS_INDEX_PM(halfedge, _index_t, Index)
DECLARE_HDS_INDEX_PM(vertex, _index_t, Index)
DECLARE_HDS_INDEX_PM(face, _index_t, Index)

} // end of CGAL namespace

#undef DECLARE_HDS_INDEX_PM

namespace CGAL {
// not done with macros, because HDS_edge::id does not return a
// reference
template<class CGAL_HDS_TMPLT>
struct HDS_property_map<CGAL_HDS_CLASS, boost::edge_index_t>
{
  struct bind_
  {
    typedef internal::Edge_index_accessor<
      typename boost::graph_traits<
        CGAL_HDS_CLASS
        >::edge_descriptor > type;
    typedef type const_type;
  };
};

template<class CGAL_HDS_TMPLT>
struct HDS_property_map<CGAL_HDS_CLASS, boost::edge_weight_t>
{
  struct bind_
  {
    typedef typename CGAL_HDS_CLASS::Traits::FT FT;
    typedef typename boost::graph_traits<CGAL_HDS_CLASS >::edge_descriptor edge_descriptor;
    typedef internal::HDS_wrap_squared<edge_descriptor,FT> type;
    typedef type const_type;
  };
};

template<class CGAL_HDS_TMPLT>
struct HDS_property_map<CGAL_HDS_CLASS,vertex_point_t>
{
  struct bind_
  {
    typedef internal::Point_accessor<
      typename boost::graph_traits<
        CGAL_HDS_CLASS
        >::vertex_descriptor,
      typename Gt::Point_3, typename Gt::Point_3&> type;

    typedef internal::Point_accessor<
      typename boost::graph_traits<
        CGAL_HDS_CLASS
        >::vertex_descriptor,
      typename Gt::Point_3, const typename Gt::Point_3&> const_type;
  };
};

//
// external indices
//

template<class CGAL_HDS_TMPLT>
struct HDS_property_map<CGAL_HDS_CLASS, edge_external_index_t>
{
  struct bind_
  {
    typedef internal::HDS_edge_index_map_external<
      CGAL_HDS_CLASS
      > type;
    typedef type const_type;
  };
};

template<class CGAL_HDS_TMPLT>
struct HDS_property_map<CGAL_HDS_CLASS, halfedge_external_index_t>
{
  struct bind_
  {
    typedef internal::HDS_index_map_external<
      typename boost::graph_traits<
        CGAL_HDS_CLASS
        >::halfedge_descriptor > type;
    typedef type const_type;
  };
};


template<class CGAL_HDS_TMPLT>
struct HDS_property_map<CGAL_HDS_CLASS, vertex_external_index_t>
{
  struct bind_
  {
    typedef internal::HDS_index_map_external<
      typename boost::graph_traits<
        CGAL_HDS_CLASS
        >::vertex_descriptor > type;
    typedef type const_type;
  };
};

template<class CGAL_HDS_TMPLT>
struct HDS_property_map<CGAL_HDS_CLASS, face_external_index_t>
{
  struct bind_
  {
    typedef internal::HDS_index_map_external<
      typename boost::graph_traits<
        CGAL_HDS_CLASS
        >::face_descriptor > type;
    typedef type const_type;
  };
};

template<class CGAL_HDS_TMPLT>
typename boost::property_map<CGAL_HDS_CLASS, boost::edge_external_index_t >::const_type
get(boost::edge_external_index_t,CGAL_HDS_CLASS const& p)
{
  return typename boost::property_map<CGAL_HDS_CLASS, boost::edge_external_index_t >::const_type(
    const_cast<CGAL_HDS_CLASS& >(p));
}

template<class CGAL_HDS_TMPLT>
typename boost::property_map<CGAL_HDS_CLASS, boost::halfedge_external_index_t >::const_type
get(boost::halfedge_external_index_t,CGAL_HDS_CLASS const& p)
{
 CGAL_HDS_CLASS& ncp = const_cast<CGAL_HDS_CLASS&>(p);

  return typename boost::property_map<CGAL_HDS_CLASS, boost::halfedge_external_index_t >::const_type(
    ncp.halfedges_begin(), ncp.halfedges_end(), ncp.size_of_halfedges());
}

template<class CGAL_HDS_TMPLT>
typename boost::property_map<CGAL_HDS_CLASS, boost::vertex_external_index_t >::const_type
get(boost::vertex_external_index_t,CGAL_HDS_CLASS const& p)
{
 CGAL_HDS_CLASS& ncp = const_cast<CGAL_HDS_CLASS&>(p);

  return typename boost::property_map<CGAL_HDS_CLASS, boost::vertex_external_index_t >::const_type(
    ncp.vertices_begin(), ncp.vertices_end(), ncp.size_of_vertices());
}

template<class CGAL_HDS_TMPLT>
typename boost::property_map<CGAL_HDS_CLASS, boost::face_external_index_t >::const_type
get(boost::face_external_index_t,CGAL_HDS_CLASS const& p)
{
 CGAL_HDS_CLASS& ncp = const_cast<CGAL_HDS_CLASS&>(p);

  return typename boost::property_map<CGAL_HDS_CLASS, boost::face_external_index_t >::const_type(
    ncp.facets_begin(), ncp.facets_end(), ncp.size_of_facets());
}

// the same blurb for non-const

template<class CGAL_HDS_TMPLT>
typename boost::property_map<CGAL_HDS_CLASS, boost::edge_external_index_t >::type
get(boost::edge_external_index_t,CGAL_HDS_CLASS& p)
{
  return typename boost::property_map<CGAL_HDS_CLASS, boost::edge_external_index_t >::type(
    p);
}

template<class CGAL_HDS_TMPLT>
typename boost::property_map<CGAL_HDS_CLASS, boost::halfedge_external_index_t >::type
get(boost::halfedge_external_index_t,CGAL_HDS_CLASS & ncp)
{
  return typename boost::property_map<CGAL_HDS_CLASS, boost::halfedge_external_index_t >::type(
    ncp.halfedges_begin(), ncp.halfedges_end(), ncp.size_of_halfedges());
}

template<class CGAL_HDS_TMPLT>
typename boost::property_map<CGAL_HDS_CLASS, boost::vertex_external_index_t >::type
get(boost::vertex_external_index_t,CGAL_HDS_CLASS & ncp)
{
  return typename boost::property_map<CGAL_HDS_CLASS, boost::vertex_external_index_t >::type(
    ncp.vertices_begin(), ncp.vertices_end(), ncp.size_of_vertices());
}


template<class CGAL_HDS_TMPLT>
typename boost::property_map<CGAL_HDS_CLASS, boost::face_external_index_t >::type
get(boost::face_external_index_t,CGAL_HDS_CLASS & ncp)
{
  return typename boost::property_map<CGAL_HDS_CLASS, boost::face_external_index_t >::type(
    ncp.facets_begin(), ncp.facets_end(), ncp.size_of_facets());
}

} // end of CGAL namespace


namespace boost {

// property_map dispatcher into Polyhedron
template<class CGAL_HDS_TMPLT, class Tag>
struct property_map<CGAL_HDS_CLASS, Tag>
{
  typedef typename CGAL::HDS_property_map<CGAL_HDS_CLASS, Tag>::
      bind_ map_gen;
  typedef typename map_gen::type       type;
  typedef typename map_gen::const_type const_type;
};

// property_map dispatcher into const Polyhedron
template<class CGAL_HDS_TMPLT, class Tag>
struct property_map<const CGAL_HDS_CLASS, Tag>
{
  typedef typename CGAL::HDS_property_map<CGAL_HDS_CLASS, Tag>::
    bind_ map_gen;
  typedef typename map_gen::type       type;
  typedef typename map_gen::const_type const_type;
};

template<class CGAL_HDS_TMPLT, class T>
struct property_map<CGAL_HDS_CLASS, CGAL::dynamic_vertex_property_t<T> >
{
  typedef CGAL_HDS_CLASS G;
  typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
  typedef CGAL::internal::Dynamic_property_map<vertex_descriptor,T> type;
  typedef type const_type;
};

template<class CGAL_HDS_TMPLT, class T>
struct property_map<CGAL_HDS_CLASS, CGAL::dynamic_halfedge_property_t<T> >
{
  typedef CGAL_HDS_CLASS G;
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<halfedge_descriptor,T> type;
  typedef type const_type;
};

template<class CGAL_HDS_TMPLT, class T>
struct property_map<CGAL_HDS_CLASS, CGAL::dynamic_edge_property_t<T> >
{
  typedef CGAL_HDS_CLASS G;
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<edge_descriptor,T> type;
  typedef type const_type;
};

template<class CGAL_HDS_TMPLT, class T>
struct property_map<CGAL_HDS_CLASS, CGAL::dynamic_face_property_t<T> >
{
  typedef CGAL_HDS_CLASS G;
  typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;
  typedef CGAL::internal::Dynamic_property_map<face_descriptor,T> type;
  typedef type const_type;
};

// What are those needed for ???
template<class CGAL_HDS_TMPLT>
struct edge_property_type<CGAL_HDS_CLASS >
{
  typedef edge_weight_t type;
};

template<class CGAL_HDS_TMPLT>
struct vertex_property_type<CGAL_HDS_CLASS >
{
  typedef CGAL::vertex_point_t type;
};

template<class CGAL_HDS_TMPLT>
struct vertex_property_type<const CGAL_HDS_CLASS >
{
  typedef CGAL::vertex_point_t type;
};

} // end of boost namespace

namespace CGAL{
template<class CGAL_HDS_TMPLT>
struct graph_has_property<CGAL_HDS_CLASS, boost::vertex_point_t>
  : CGAL::Tag_true {};

template<class CGAL_HDS_TMPLT>
struct graph_has_property<CGAL_HDS_CLASS, boost::edge_weight_t>
  : CGAL::Tag_true {};

template<class CGAL_HDS_TMPLT>
struct graph_has_property<CGAL_HDS_CLASS, boost::edge_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename boost::graph_traits<CGAL_HDS_CLASS >::edge_descriptor
      >::value
    >
{};

template<class CGAL_HDS_TMPLT>
struct graph_has_property<CGAL_HDS_CLASS, boost::face_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL_HDS_CLASS::Facet
      >::value
    >
{};

template<class CGAL_HDS_TMPLT>
struct graph_has_property<CGAL_HDS_CLASS, boost::halfedge_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL_HDS_CLASS::Halfedge
      >::value
    >
{};

template<class CGAL_HDS_TMPLT>
struct graph_has_property<CGAL_HDS_CLASS, boost::vertex_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL_HDS_CLASS::Vertex
      >::value
    >
{};
}// end of CGAL namespace

#undef CGAL_HDS_TMPLT
#undef CGAL_HDS_CLASS
