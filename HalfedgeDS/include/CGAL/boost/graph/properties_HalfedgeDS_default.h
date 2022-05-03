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

#ifndef CGAL_BOOST_GRAPH_PROPERTIES_HALFEDGEDS_DEFAULT_H
#define CGAL_BOOST_GRAPH_PROPERTIES_HALFEDGEDS_DEFAULT_H

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/graph_traits_HalfedgeDS_default.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/number_utils.h>
#include <memory>
#include <CGAL/boost/graph/internal/Has_member_id.h>
#include <CGAL/Distance_3/Point_3_Point_3.h>



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

} // internal

// the tag we dispatch on from property_map<G, Property>
template <class Tag>
struct HDS_property_map {};

} // namespace CGAL

namespace CGAL {
    /*
// generalized 2-ary get functions
template<class Gt, class I,  class A, class PropertyTag>
typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag >::const_type
get(PropertyTag,CGAL::HalfedgeDS_default<Gt,I,A> const&)
{ return typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag >::const_type(); }

template<class Gt, class I,  class A, class PropertyTag>
typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag >::type
get(PropertyTag,CGAL::HalfedgeDS_default<Gt,I,A>&)
{ return typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag >::type(); }
*/
// generalized 3-ary get functions
template<class Gt, class I,  class A, class PropertyTag, class Key>
typename boost::property_traits< typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag >::type >::reference
get(PropertyTag p,CGAL::HalfedgeDS_default<Gt,I,A>& g, const Key& key)
{ return get(get(p, g), key); }

template<class Gt, class I,  class A, class PropertyTag, class Key>
typename boost::property_traits< typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag >::const_type >::reference
get(PropertyTag p,CGAL::HalfedgeDS_default<Gt,I,A> const& g, const Key& key)
{ return get(get(p, g), key); }



#define CGAL_POLYHEDRON_DYNAMIC_PM(TAG, DESCRIPTOR) \
template <typename T, typename Gt, typename I,  typename A> \
typename boost::property_map<HalfedgeDS_default<Gt,I,A>, TAG >::const_type \
get(const TAG&, const HalfedgeDS_default<Gt,I,A>&) \
{ \
  typedef typename boost::graph_traits< HalfedgeDS_default<Gt,I,A> >::DESCRIPTOR descriptor; \
  return internal::Dynamic_property_map<descriptor,T>(); \
}

CGAL_POLYHEDRON_DYNAMIC_PM(dynamic_vertex_property_t<T>, vertex_descriptor)
CGAL_POLYHEDRON_DYNAMIC_PM(dynamic_halfedge_property_t<T>, halfedge_descriptor)
CGAL_POLYHEDRON_DYNAMIC_PM(dynamic_edge_property_t<T>, edge_descriptor)
CGAL_POLYHEDRON_DYNAMIC_PM(dynamic_face_property_t<T>, face_descriptor)


#undef CGAL_POLYHEDRON_DYNAMIC_PM


// generalized put
template<class Gt, class I,  class A, class PropertyTag, class Key,class Value>
void put(PropertyTag p,CGAL::HalfedgeDS_default<Gt,I,A>& g, const Key& key, const Value& value)
{
  typedef typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag>::type Map;
  Map pmap = get(p, g);
  put(pmap, key, value);
}

} // CGAL

// specialization needs to be repeated for halfedge, vertex, face
#define CGAL_POLYHEDRON_INDEX_PM(ENTITY, TAG, ACCESSOR)                 \
  namespace CGAL {                                                      \
  template<> struct HDS_property_map<boost::ENTITY##TAG> {  \
  template<class Gt, class I,  class A>                 \
  struct bind_ {                                                        \
    typedef internal::ACCESSOR##_accessor<                              \
    CGAL::HalfedgeDS_default<Gt, I, A>, \
    typename boost::graph_traits< CGAL::HalfedgeDS_default<Gt, I, A>     \
                                  >::ENTITY##_descriptor > type;        \
    typedef type const_type;                                            \
  };                                                                    \
  };                                                                    \
  } //CGAL

CGAL_POLYHEDRON_INDEX_PM(halfedge, _index_t, Index)
CGAL_POLYHEDRON_INDEX_PM(vertex, _index_t, Index)
CGAL_POLYHEDRON_INDEX_PM(face, _index_t, Index)

#undef CGAL_POLYHEDRON_INDEX_PM

namespace CGAL {
// not done with macros, because HDS_edge::id does not return a
// reference
template <>
struct HDS_property_map<boost::edge_index_t>
{
  template<class Gt, class I,  class A>
  struct bind_
  {
    typedef internal::Edge_index_accessor<
      typename boost::graph_traits<
        CGAL::HalfedgeDS_default<Gt, I, A>
        >::edge_descriptor > type;
    typedef type const_type;
  };
};

template <>
struct HDS_property_map<boost::edge_weight_t>
{
  template<class Gt, class I,  class A>
  struct bind_
  {
    typedef typename CGAL::HalfedgeDS_default<Gt,I,A>::Traits::FT FT;
    typedef typename boost::graph_traits<CGAL::HalfedgeDS_default<Gt,I,A> >::edge_descriptor edge_descriptor;
    typedef internal::HDS_wrap_squared<edge_descriptor,FT> type;
    typedef type const_type;
  };
};
/*
// already defined in line 448 ??    why not the same for Polyhedron?
template <>
struct HDS_property_map<vertex_point_t>
{
  template<class Gt, class I,  class A>
  struct bind_
  {
    typedef internal::Point_accessor<
      typename boost::graph_traits<
        CGAL::HalfedgeDS_default<Gt,I,A>
        >::vertex_descriptor,
      typename Gt::Point_3, typename Gt::Point_3&> type;

    typedef internal::Point_accessor<
      typename boost::graph_traits<
        CGAL::HalfedgeDS_default<Gt,I,A>
        >::vertex_descriptor,
      typename Gt::Point_3, const typename Gt::Point_3&> const_type;
  };
};
*/
//
// external indices
//

template <>
struct HDS_property_map<edge_external_index_t>
{
  template<class Gt, class I,  class A>
  struct bind_
  {
    typedef internal::HDS_edge_index_map_external<
      CGAL::HalfedgeDS_default<Gt,I,A>
      > type;
    typedef type const_type;
  };
};

template <>
struct HDS_property_map<halfedge_external_index_t>
{
  template<class Gt, class I,  class A>
  struct bind_
  {
    typedef internal::HDS_index_map_external<
      typename boost::graph_traits<
        CGAL::HalfedgeDS_default<Gt,I,A>
        >::halfedge_descriptor > type;
    typedef type const_type;
  };
};


template <>
struct HDS_property_map<vertex_external_index_t>
{
  template<class Gt, class I,  class A>
  struct bind_
  {
    typedef internal::HDS_index_map_external<
      typename boost::graph_traits<
        CGAL::HalfedgeDS_default<Gt,I,A>
        >::vertex_descriptor > type;
    typedef type const_type;
  };
};

template <>
struct HDS_property_map<face_external_index_t>
{
  template<class Gt, class I,  class A>
  struct bind_
  {
    typedef internal::HDS_index_map_external<
      typename boost::graph_traits<
        CGAL::HalfedgeDS_default<Gt,I,A>
        >::face_descriptor > type;
    typedef type const_type;
  };
};

} // CGAL

namespace CGAL {

template<class Gt, class I,  class A>
typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::edge_external_index_t >::const_type
get(boost::edge_external_index_t,CGAL::HalfedgeDS_default<Gt,I,A> const& p)
{
  return typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::edge_external_index_t >::const_type(
    const_cast<CGAL::HalfedgeDS_default<Gt,I,A>& >(p));
}

template<class Gt, class I,  class A>
typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::halfedge_external_index_t >::const_type
get(boost::halfedge_external_index_t,CGAL::HalfedgeDS_default<Gt,I,A> const& p)
{
 CGAL::HalfedgeDS_default<Gt,I,A>& ncp = const_cast<CGAL::HalfedgeDS_default<Gt,I,A>&>(p);

  return typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::halfedge_external_index_t >::const_type(
    ncp.halfedges_begin(), ncp.halfedges_end(), ncp.size_of_halfedges());
}

template<class Gt, class I,  class A>
typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::vertex_external_index_t >::const_type
get(boost::vertex_external_index_t,CGAL::HalfedgeDS_default<Gt,I,A> const& p)
{
 CGAL::HalfedgeDS_default<Gt,I,A>& ncp = const_cast<CGAL::HalfedgeDS_default<Gt,I,A>&>(p);

  return typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::vertex_external_index_t >::const_type(
    ncp.vertices_begin(), ncp.vertices_end(), ncp.size_of_vertices());
}

template<class Gt, class I,  class A>
typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::face_external_index_t >::const_type
get(boost::face_external_index_t,CGAL::HalfedgeDS_default<Gt,I,A> const& p)
{
 CGAL::HalfedgeDS_default<Gt,I,A>& ncp = const_cast<CGAL::HalfedgeDS_default<Gt,I,A>&>(p);

  return typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::face_external_index_t >::const_type(
    ncp.facets_begin(), ncp.facets_end(), ncp.size_of_facets());
}

// the same blurb for non-const

template<class Gt, class I,  class A>
typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::edge_external_index_t >::type
get(boost::edge_external_index_t,CGAL::HalfedgeDS_default<Gt,I,A>& p)
{
  return typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::edge_external_index_t >::type(
    p);
}

template<class Gt, class I,  class A>
typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::halfedge_external_index_t >::type
get(boost::halfedge_external_index_t,CGAL::HalfedgeDS_default<Gt,I,A> & ncp)
{
  return typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::halfedge_external_index_t >::type(
    ncp.halfedges_begin(), ncp.halfedges_end(), ncp.size_of_halfedges());
}

template<class Gt, class I,  class A>
typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::vertex_external_index_t >::type
get(boost::vertex_external_index_t,CGAL::HalfedgeDS_default<Gt,I,A> & ncp)
{
  return typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::vertex_external_index_t >::type(
    ncp.vertices_begin(), ncp.vertices_end(), ncp.size_of_vertices());
}


template<class Gt, class I,  class A>
typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::face_external_index_t >::type
get(boost::face_external_index_t,CGAL::HalfedgeDS_default<Gt,I,A> & ncp)
{
  return typename boost::property_map<CGAL::HalfedgeDS_default<Gt,I,A>, boost::face_external_index_t >::type(
    ncp.facets_begin(), ncp.facets_end(), ncp.size_of_facets());
}



} // namespace CGAL


namespace boost {

// property_map dispatcher into Polyhedron
template<class Gt, class I,  class A, class Tag>
struct property_map<CGAL::HalfedgeDS_default<Gt,I,A>, Tag>
{
  typedef typename CGAL::HDS_property_map<Tag>::
      template bind_<Gt,I,A> map_gen;
  typedef typename map_gen::type       type;
  typedef typename map_gen::const_type const_type;
};

// property_map dispatcher into const Polyhedron
template<class Gt, class I,  class A, class Tag>
struct property_map<const CGAL::HalfedgeDS_default<Gt,I,A>, Tag>
{
  typedef typename CGAL::HDS_property_map<Tag>::
      template bind_<Gt,I,A> map_gen;
  typedef typename map_gen::type       type;
  typedef typename map_gen::const_type const_type;
};

template<class Gt, class I,  class A, class T>
struct property_map<CGAL::HalfedgeDS_default<Gt,I,A>, CGAL::dynamic_vertex_property_t<T> >
{
  typedef CGAL::HalfedgeDS_default<Gt,I,A> G;
  typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
  typedef CGAL::internal::Dynamic_property_map<vertex_descriptor,T> type;
  typedef type const_type;
};

template<class Gt, class I,  class A, class T>
struct property_map<CGAL::HalfedgeDS_default<Gt,I,A>, CGAL::dynamic_halfedge_property_t<T> >
{
  typedef CGAL::HalfedgeDS_default<Gt,I,A> G;
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<halfedge_descriptor,T> type;
  typedef type const_type;
};

template<class Gt, class I,  class A, class T>
struct property_map<CGAL::HalfedgeDS_default<Gt,I,A>, CGAL::dynamic_edge_property_t<T> >
{
  typedef CGAL::HalfedgeDS_default<Gt,I,A> G;
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef CGAL::internal::Dynamic_property_map<edge_descriptor,T> type;
  typedef type const_type;
};

template<class Gt, class I,  class A, class T>
struct property_map<CGAL::HalfedgeDS_default<Gt,I,A>, CGAL::dynamic_face_property_t<T> >
{
  typedef CGAL::HalfedgeDS_default<Gt,I,A> G;
  typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;
  typedef CGAL::internal::Dynamic_property_map<face_descriptor,T> type;
  typedef type const_type;
};

// What are those needed for ???
template<class Gt, class I,  class A>
struct edge_property_type<CGAL::HalfedgeDS_default<Gt,I,A> >
{
  typedef edge_weight_t type;
};

template<class Gt, class I,  class A>
struct vertex_property_type<CGAL::HalfedgeDS_default<Gt,I,A> >
{
  typedef CGAL::vertex_point_t type;
};

template<class Gt, class I,  class A>
struct vertex_property_type<const CGAL::HalfedgeDS_default<Gt,I,A> >
{
  typedef CGAL::vertex_point_t type;
};


} // namespace boost

namespace CGAL{
template<class Gt, class I,  class A>
struct graph_has_property<CGAL::HalfedgeDS_default<Gt,I,A>, boost::vertex_point_t>
  : CGAL::Tag_true {};

template<class Gt, class I,  class A>
struct graph_has_property<CGAL::HalfedgeDS_default<Gt,I,A>, boost::edge_weight_t>
  : CGAL::Tag_true {};

template<class Gt, class I,  class A>
struct graph_has_property<CGAL::HalfedgeDS_default<Gt,I,A>, boost::edge_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename boost::graph_traits<CGAL::HalfedgeDS_default<Gt,I,A> >::edge_descriptor
      >::value
    >
{};

template<class Gt, class I,  class A>
struct graph_has_property<CGAL::HalfedgeDS_default<Gt,I,A>, boost::face_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL::HalfedgeDS_default<Gt,I,A>::Facet
      >::value
    >
{};

template<class Gt, class I,  class A>
struct graph_has_property<CGAL::HalfedgeDS_default<Gt,I,A>, boost::halfedge_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL::HalfedgeDS_default<Gt,I,A>::Halfedge
      >::value
    >
{};

template<class Gt, class I,  class A>
struct graph_has_property<CGAL::HalfedgeDS_default<Gt,I,A>, boost::vertex_index_t>
  : CGAL::Boolean_tag<
      CGAL::internal::Has_member_id<
        typename CGAL::HalfedgeDS_default<Gt,I,A>::Vertex
      >::value
    >
{};
}// end CGAL
#undef CGAL_HDS_PARAM_




#endif // CGAL_BOOST_GRAPH_PROPERTIES_HALFEDGEDS_DEFAULT_H
