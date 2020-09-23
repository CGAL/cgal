// Copyright (c) 2018  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_GRAPH_TRAITS_HALFEDGEDS_DEFAULT_H
#define CGAL_GRAPH_TRAITS_HALFEDGEDS_DEFAULT_H

#include <CGAL/boost/graph/graph_traits_HalfedgeDS.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/boost/graph/properties.h>

namespace CGAL {

template <class Traits_, class HalfedgeDSItems,
          class Alloc>
class HalfedgeDS_default;
} // namespace CGAL

namespace boost {



template<class T, class I, class A>
struct graph_traits< CGAL::HalfedgeDS_default<T,I,A> >
   : CGAL::HDS_graph_traits< CGAL::HalfedgeDS_default<T,I,A> >
{
  typedef typename T::Point_3 vertex_property_type;
};

template<class T, class I, class A>
struct graph_traits< CGAL::HalfedgeDS_default<T,I,A> const >
   : CGAL::HDS_graph_traits< CGAL::HalfedgeDS_default<T,I,A> > // See NOTE above!
{};

} // namespace boost

namespace CGAL {

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertices_size_type
num_vertices(const HalfedgeDS_default<T,I,A>& p)
{
  return p.size_of_vertices();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::edges_size_type
num_edges(const HalfedgeDS_default<T,I,A>& p)
{
  return p.size_of_halfedges() / 2;
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::degree_size_type
degree(typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_descriptor v
       , const HalfedgeDS_default<T,I,A>&)
{
  return v->vertex_degree();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::degree_size_type
out_degree(typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_descriptor v
           , const HalfedgeDS_default<T,I,A>&)
{
  return v->vertex_degree();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::degree_size_type
in_degree(typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_descriptor v
          , const HalfedgeDS_default<T,I,A>&)
{
  return v->vertex_degree();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::degree_size_type
degree(typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::face_descriptor f
       , const HalfedgeDS_default<T,I,A>&)
{
  return f->facet_degree();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_descriptor
source(typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::edge_descriptor e
       , const HalfedgeDS_default<T,I,A> & )
{
  return e.halfedge()->opposite()->vertex();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_descriptor
target(typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::edge_descriptor e
       , const HalfedgeDS_default<T,I,A> & )
{
  return e.halfedge()->vertex();
}

template<class T, class I, class A>
std::pair<
  typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::edge_descriptor
  , bool>
edge(typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_descriptor u
     , typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_descriptor v
     , const HalfedgeDS_default<T,I,A> &)
{
  typedef HalfedgeDS_default<T,I,A> P;
  typedef boost::graph_traits< P > Traits;
  typedef typename Traits::halfedge_descriptor halfedge_descriptor;
    typedef typename Traits::edge_descriptor edge;

  // circulate around the inedges of u
  halfedge_descriptor c(u->halfedge()), d(u->halfedge());
  if(c != halfedge_descriptor()) {
    do {
      if(c->opposite()->vertex() == v) {
        return std::make_pair(edge(c->opposite()), true);
      }
      c = c->next()->opposite();
    } while (c != d);
  }

  return std::make_pair(edge(), false);
}

template<class T, class I, class A>
inline Iterator_range<typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_iterator>
vertices( const HalfedgeDS_default<T,I,A>& p)
{
  typedef typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_iterator Iter;
  HalfedgeDS_default<T,I,A>& ncp = const_cast<HalfedgeDS_default<T,I,A>&>(p);
  return make_range( Iter(ncp.vertices_begin()), Iter(ncp.vertices_end()) );
}

template<class T, class I, class A>
inline Iterator_range<typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::edge_iterator>
edges( const HalfedgeDS_default<T,I,A>& p)
{
  typedef typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::edge_iterator_i Iter_i;
  typedef typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::edge_iterator Iter;
  HalfedgeDS_default<T,I,A>& ncp = const_cast<HalfedgeDS_default<T,I,A>&>(p);
  return make_range( Iter(Iter_i(ncp.halfedges_begin())), Iter(Iter_i(ncp.halfedges_end()) ));
}

template<class T, class I, class A>
inline Iterator_range<typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::in_edge_iterator>
in_edges( typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_descriptor u
          , const HalfedgeDS_default<T,I,A>& p)
{
  typedef typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::in_edge_iterator Iter;
  return make_range(Iter(halfedge(u,p),p), Iter(halfedge(u,p),p,1));
}

template<class T, class I, class A>
inline Iterator_range<typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::out_edge_iterator>
out_edges( typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertex_descriptor u
           , const HalfedgeDS_default<T,I,A>& p)
{
  typedef typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::out_edge_iterator Iter;
  return make_range(Iter(halfedge(u,p),p), Iter(halfedge(u,p),p,1));
}

//
// MutableHalfedgeGraph
//

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::vertex_descriptor
add_vertex(HalfedgeDS_default<T,I,A>& g)
{
  return g.vertices_push_back(typename HalfedgeDS_default<T,I,A>::Vertex());
}


template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::vertex_descriptor
add_vertex(const typename boost::graph_traits<HalfedgeDS_default<T,I,A> >::vertex_property_type& p
           , HalfedgeDS_default<T,I,A>& g)
{
  return g.vertices_push_back(typename HalfedgeDS_default<T,I,A>::Vertex(p));
}

template<class T, class I, class A>
void
remove_vertex(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::vertex_descriptor v
             , HalfedgeDS_default<T,I,A>& g)
{
  g.vertices_erase(v);
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::edge_descriptor
add_edge(HalfedgeDS_default<T,I,A>& g)
{
  return typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::edge_descriptor(
    g.edges_push_back(typename HalfedgeDS_default<T,I,A>::Halfedge(),
                            typename HalfedgeDS_default<T,I,A>::Halfedge()));
}

template<class T, class I, class A>
void
remove_edge(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::edge_descriptor e
            , HalfedgeDS_default<T,I,A>& g)
{
  g.edges_erase(e.halfedge());
}


template<class T, class I, class A>
void
set_target(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h1
         , typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::vertex_descriptor v
         , HalfedgeDS_default<T,I,A>&)
{
  // set_face has become private in the halfedge provided by
  // polyhedron for unknown reasons, although it used to be public
  // once.

  // We sneak in anyway. Inheritance can't keep us out.
  typedef typename HalfedgeDS_default<T,I,A>::Halfedge::Base Sneak;
  static_cast<Sneak&>(*h1).set_vertex(v);
}

template<class T, class I, class A>
void
set_next(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h1
         , typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h2
         , HalfedgeDS_default<T,I,A>&)
{
  typedef typename HalfedgeDS_default<T,I,A>::Halfedge::Base Sneak;
  static_cast<Sneak&>(*h1).set_next(h2);
  static_cast<Sneak&>(*h2).set_prev(h1);
}

//
// MutableFaceGraph
//

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::face_descriptor
add_face(HalfedgeDS_default<T,I,A>& g)
{
  return g.faces_push_back(typename HalfedgeDS_default<T,I,A>::Face());
}

template<class T, class I, class A, class InputIterator>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::face_descriptor
add_face(InputIterator begin, InputIterator end, HalfedgeDS_default<T,I,A>& g)
{
  return g.add_facet(begin, end);
}

template<class T, class I, class A>
void
remove_face(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::face_descriptor f
            , HalfedgeDS_default<T,I,A>& g)
{
  g.faces_erase(f);
}

template<class T, class I, class A>
void
set_face(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h
  , typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::face_descriptor f
  , const HalfedgeDS_default<T,I,A>&)
{
  // set_face has become private in the halfedge provided by
  // polyhedron for unknown reasons, although it used to be public
  // once.

  // We sneak in anyway. Inheritance can't keep us out.
  typedef typename HalfedgeDS_default<T,I,A>::Halfedge::Base Sneak;
  static_cast<Sneak&>(*h).set_face(f);
}

template<class T, class I, class A>
void
set_halfedge(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::face_descriptor f
  , typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h
  , HalfedgeDS_default<T,I,A>& g)
{
  typedef HalfedgeDS_default<T,I,A> Hds;
  CGAL::HalfedgeDS_decorator<Hds> D(g);
  D.set_face_halfedge(f, h);
}

template<class T, class I, class A>
void
set_halfedge(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::vertex_descriptor v
  , typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h
  , const HalfedgeDS_default<T,I,A>&)
{
  typedef typename HalfedgeDS_default<T,I,A>::Vertex::Base Sneak;
  static_cast<Sneak&>(*v).set_halfedge(h);
}


//
// HalfedgeGraph
//
template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::edge_descriptor
edge(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h
     , const HalfedgeDS_default<T,I,A>&)
{
  return typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::edge_descriptor(h);
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor
halfedge(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::edge_descriptor e
         , const HalfedgeDS_default<T,I,A>&)
{
  return e.halfedge();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor
halfedge(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::vertex_descriptor v
         , const HalfedgeDS_default<T,I,A>&)
{
  return v->halfedge();
}

template<class T, class I, class A>
std::pair< typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor
           , bool>
halfedge(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::vertex_descriptor u
         , typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::vertex_descriptor v
         , const HalfedgeDS_default<T,I,A>& g)
{
  std::pair< typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::edge_descriptor
             , bool> e = edge(u, v, g);
  return std::make_pair(e.first.halfedge(), e.second);
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor
opposite(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h
         , const HalfedgeDS_default<T,I,A>&)
{
  return h->opposite();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::vertex_descriptor
source(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h
         , const HalfedgeDS_default<T,I,A>& g)
{
  return target(opposite(h, g), g);
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::vertex_descriptor
target(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h
         , const HalfedgeDS_default<T,I,A>&)
{
  return h->vertex();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor
next(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor outedge
     , const HalfedgeDS_default<T,I,A>&)
{
  return outedge->next();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor
prev(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor outedge
     , const HalfedgeDS_default<T,I,A>&)
{
  return outedge->prev();
}


//
// HalfedgeListGraph
//

template<class T, class I, class A>
Iterator_range<typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_iterator>
halfedges(const HalfedgeDS_default<T,I,A>& p)
{
  typedef typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_iterator Iter;
  HalfedgeDS_default<T,I,A>& ncp = const_cast<HalfedgeDS_default<T,I,A>&>(p);
  return make_range(Iter(ncp.halfedges_begin()), Iter(ncp.halfedges_end()));
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedges_size_type
num_halfedges(const HalfedgeDS_default<T,I,A>& p)
{
  return p.size_of_halfedges();
}

// FaceGraph
template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::face_descriptor
face(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor h
     , const HalfedgeDS_default<T,I,A>&)
{
  return h->face();
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::halfedge_descriptor
halfedge(typename boost::graph_traits< HalfedgeDS_default<T,I,A> >::face_descriptor f
         , const HalfedgeDS_default<T,I,A>&)
{
  return f->halfedge();
}

template<class T, class I, class A>
inline Iterator_range<typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::face_iterator >
faces(const HalfedgeDS_default<T,I,A>& p)
{
  typedef typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::face_iterator face_iterator;
  HalfedgeDS_default<T,I,A>& ncp = const_cast<HalfedgeDS_default<T,I,A>&>(p);
  return make_range( face_iterator(ncp.faces_begin()), face_iterator(ncp.faces_end()));
}

template<class T, class I, class A>
typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::faces_size_type
num_faces(const HalfedgeDS_default<T,I,A>& p)
{
  return p.size_of_faces();
}

template <class T>
struct HDS_property_map;

template <>
struct HDS_property_map<vertex_point_t>
{
  template<class T, class I, class A>
  struct bind_
  {
    typedef internal::Point_accessor<
      typename boost::graph_traits<
        HalfedgeDS_default<T, I, A>
        >::vertex_descriptor,
      typename T::Point_3, typename T::Point_3&> type;

    typedef internal::Point_accessor<
      typename boost::graph_traits<
        HalfedgeDS_default<T, I, A>
        >::vertex_descriptor,
      typename T::Point_3, const typename T::Point_3&> const_type;
  };
};

template<class T, class I, class A>
void reserve(HalfedgeDS_default<T,I,A>& p,
             typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::vertices_size_type nv,
             typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::edges_size_type ne,
             typename boost::graph_traits< HalfedgeDS_default<T,I,A> const>::faces_size_type nf)
{
  p.reserve(nv, 2*ne, nf);
}

}// namespace CGAL
namespace boost {

#define CGAL_PM_SPECIALIZATION(TAG) \
template<class T, class I, class A> \
struct property_map<CGAL::HalfedgeDS_default<T,I,A>, TAG> \
{\
  typedef typename CGAL::HDS_property_map<TAG>:: \
      template bind_<T,I,A> map_gen; \
  typedef typename map_gen::type       type; \
  typedef typename map_gen::const_type const_type; \
};

CGAL_PM_SPECIALIZATION(vertex_point_t)

#undef CGAL_PM_SPECIALIZATION

} // namespace boost

namespace CGAL {

// generalized 2-ary get functions
template<class Gt, class I, class A, class PropertyTag>
typename boost::property_map< CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag >::const_type
get(PropertyTag, CGAL::HalfedgeDS_default<Gt,I,A> const&)
{ return typename boost::property_map< CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag >::const_type(); }

template<class Gt, class I, class A, class PropertyTag>
typename boost::property_map< CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag >::type
get(PropertyTag, CGAL::HalfedgeDS_default<Gt,I,A>&)
{ return typename boost::property_map< CGAL::HalfedgeDS_default<Gt,I,A>, PropertyTag >::type(); }


} // namespace CGAL
#endif
