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
//
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYHEDRON_3_H
#define CGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYHEDRON_3_H

#include <CGAL/boost/graph/graph_traits_HalfedgeDS.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/HalfedgeDS_decorator.h>

#define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS

//
// NOTE: The BGL algorithms are NOT const-correct: i.e., they take a "G const&"
// but instantiate "graph_traits<G>" instead of "graph_traits<G const>"
// This is known Boost bug which will eventually be fixed, but in the meantime we need
// to coerce both const and non-const specializations.
// That is, HDS_graph_traits<G const> is really the same as HDS_graph_traits<G>
// so graph_traits<G> is also the same as graph_traits<G const>.
// Therefore, while, for instance, "graph_traits<G const>::vertex_descriptor"
// is conceptually const-correct, it actually corresponds to the non-const handle,
// hence the const_cast<> used below in the functions implementation.
//

namespace boost
{

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >
   : CGAL::HDS_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >
{
  typedef typename Gt::Point_3 vertex_property_type;
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const >
   : CGAL::HDS_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> > // See NOTE above!
{};

} // namespace boost

namespace CGAL {

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertices_size_type
num_vertices(const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  return p.size_of_vertices();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edges_size_type
num_edges(const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  return p.size_of_halfedges() / 2;
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::degree_size_type
degree(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor v
       , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  return v->vertex_degree();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::degree_size_type
out_degree(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor v
           , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  return v->vertex_degree();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::degree_size_type
in_degree(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor v
          , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  return v->vertex_degree();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor
source(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor e
       , const CGAL::Polyhedron_3<Gt,I,HDS,A> & )
{
  return e.halfedge()->opposite()->vertex();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor
target(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor e
       , const CGAL::Polyhedron_3<Gt,I,HDS,A> & )
{
  return e.halfedge()->vertex();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
std::pair<
  typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_descriptor
  , bool>
edge(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor u
     , typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor v
     , const CGAL::Polyhedron_3<Gt,I,HDS,A> &)
{
  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> P;
  typedef typename P::Halfedge_around_vertex_circulator Circ;
  typedef boost::graph_traits< P > Traits;
  typedef typename Traits::edge_descriptor edge;

  // circulate around the inedges of u
  Circ c(u->halfedge()), d(u->halfedge());
  if(c != 0) {
    do {
      if(c->opposite()->vertex() == v) {
        return std::make_pair(edge(c->opposite()), true);
      }
    } while (++c != d);
  }

  return std::make_pair(edge(), false);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline Iterator_range<typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_iterator>
vertices( const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  typedef typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_iterator Iter;
  CGAL::Polyhedron_3<Gt,I,HDS,A>& ncp = const_cast<CGAL::Polyhedron_3<Gt,I,HDS,A>&>(p);
  return make_range( Iter(ncp.vertices_begin()), Iter(ncp.vertices_end()) );
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline Iterator_range<typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_iterator>
edges( const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  typedef typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_iterator_i Iter_i;
  typedef typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::edge_iterator Iter;
  CGAL::Polyhedron_3<Gt,I,HDS,A>& ncp = const_cast<CGAL::Polyhedron_3<Gt,I,HDS,A>&>(p);
  return make_range( Iter(Iter_i(ncp.halfedges_begin())), Iter(Iter_i(ncp.halfedges_end()) ));
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline Iterator_range<typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::in_edge_iterator>
in_edges( typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor u
          , const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  typedef typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::in_edge_iterator Iter;
  return make_range(Iter(halfedge(u,p),p), Iter(halfedge(u,p),p,1));
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline Iterator_range<typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::out_edge_iterator>
out_edges( typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::vertex_descriptor u
           , const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  typedef typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::out_edge_iterator Iter;
  return make_range(Iter(halfedge(u,p),p), Iter(halfedge(u,p),p,1));
}

//
// MutableHalfedgeGraph
// 

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor
add_vertex(CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{
  return g.hds().vertices_push_back(typename CGAL::Polyhedron_3<Gt,I,HDS,A>::Vertex());
}


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor
add_vertex(const typename boost::graph_traits<CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_property_type& p
           , CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{
  return g.hds().vertices_push_back(typename CGAL::Polyhedron_3<Gt,I,HDS,A>::Vertex(p));
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
void
remove_vertex(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor v
             , CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{
  g.hds().vertices_erase(v);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor
add_edge(CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{ 
  return typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor(
    g.hds().edges_push_back(typename CGAL::Polyhedron_3<Gt,I,HDS,A>::Halfedge(), 
                            typename CGAL::Polyhedron_3<Gt,I,HDS,A>::Halfedge()));
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
void
remove_edge(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor e
            , CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{
  g.hds().edges_erase(e.halfedge());
}


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
void
set_target(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h1
         , typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor v
         , CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  // set_face has become private in the halfedge provided by
  // polyhedron for unknown reasons, although it used to be public
  // once.

  // We sneak in anyway. Inheritance can't keep us out.
  typedef typename CGAL::Polyhedron_3<Gt,I,HDS,A>::Halfedge::Base Sneak;
  static_cast<Sneak&>(*h1).set_vertex(v);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
void
set_next(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h1
         , typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h2
         , CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  typedef typename CGAL::Polyhedron_3<Gt,I,HDS,A>::Halfedge::Base Sneak;
  static_cast<Sneak&>(*h1).set_next(h2);
  static_cast<Sneak&>(*h2).set_prev(h1);
}

//
// MutableFaceGraph 
//

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::face_descriptor
add_face(CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{
  return g.hds().faces_push_back(typename CGAL::Polyhedron_3<Gt,I,HDS,A>::HalfedgeDS::Face());
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A, class InputIterator>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt, I, HDS, A> >::face_descriptor
add_face(InputIterator begin, InputIterator end, CGAL::Polyhedron_3<Gt, I, HDS, A>& g)
{
  return g.hds().add_facet(begin, end);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
void
remove_face(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::face_descriptor f
            , CGAL::Polyhedron_3<Gt,I,HDS,A>& g) 
{
  g.hds().faces_erase(f);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
void
set_face(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h
  , typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::face_descriptor f
  , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  // set_face has become private in the halfedge provided by
  // polyhedron for unknown reasons, although it used to be public
  // once.

  // We sneak in anyway. Inheritance can't keep us out.
  typedef typename CGAL::Polyhedron_3<Gt,I,HDS,A>::Halfedge::Base Sneak;
  static_cast<Sneak&>(*h).set_face(f);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
void
set_halfedge(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::face_descriptor f
  , typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h
  , CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{
  typedef typename CGAL::Polyhedron_3<Gt, I, HDS, A>::Halfedge_data_structure Hds;
  CGAL::HalfedgeDS_decorator<Hds> D(g.hds());
  D.set_face_halfedge(f, h);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
void
set_halfedge(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor v
  , typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h
  , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  typedef typename CGAL::Polyhedron_3<Gt,I,HDS,A>::Vertex::Base Sneak;
  static_cast<Sneak&>(*v).set_halfedge(h);
}


//
// HalfedgeGraph
//
template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor
edge(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h
     , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{ 
  return typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor(h);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor
halfedge(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor e
         , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{ 
  return e.halfedge();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor
halfedge(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor v
         , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{ 
  return v->halfedge();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
std::pair< typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor
           , bool>
halfedge(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor u
         , typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor v
         , const CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{ 
  std::pair< typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::edge_descriptor
             , bool> e = edge(u, v, g);
  return std::make_pair(e.first.halfedge(), e.second);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor
opposite(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h
         , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  return h->opposite();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor
source(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h
         , const CGAL::Polyhedron_3<Gt,I,HDS,A>& g)
{
  return target(opposite(h, g), g);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::vertex_descriptor
target(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h
         , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  return h->vertex();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor
next(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor outedge
     , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  return outedge->next();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor
prev(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor outedge
     , const CGAL::Polyhedron_3<Gt,I,HDS,A>&)
{
  return outedge->prev();
}


//
// HalfedgeListGraph
//

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
Iterator_range<typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_iterator>
halfedges(const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  typedef typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_iterator Iter;
  CGAL::Polyhedron_3<Gt,I,HDS,A>& ncp = const_cast<CGAL::Polyhedron_3<Gt,I,HDS,A>&>(p);
  return make_range(Iter(ncp.halfedges_begin()), Iter(ncp.halfedges_end()));
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedges_size_type
num_halfedges(const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  return p.size_of_halfedges();
}

// FaceGraph
template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::face_descriptor
face(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor h
     , const CGAL::Polyhedron_3<Gt,I,HDS,A>&) 
{
  return h->face();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::halfedge_descriptor
halfedge(typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >::face_descriptor f
         , const CGAL::Polyhedron_3<Gt,I,HDS,A>&) 
{
  return f->halfedge();
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline Iterator_range<typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::face_iterator >
faces(const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  typedef typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::face_iterator face_iterator;
  CGAL::Polyhedron_3<Gt,I,HDS,A>& ncp = const_cast<CGAL::Polyhedron_3<Gt,I,HDS,A>&>(p);
  return make_range( face_iterator(ncp.facets_begin()), face_iterator(ncp.facets_end()));
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
typename boost::graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> const>::faces_size_type
num_faces(const CGAL::Polyhedron_3<Gt,I,HDS,A>& p)
{
  return p.size_of_facets();
}



template<class Gt, class I, CGAL_HDS_PARAM_, class A>
bool is_valid(const CGAL::Polyhedron_3<Gt,I,HDS,A>& p, bool verbose = false)
{
  return p.is_valid(verbose);
}
} // namespace CGAL


#ifndef CGAL_NO_DEPRECATED_CODE

namespace CGAL {
template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct halfedge_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >
   : CGAL::HDS_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> >
{
  typedef CGAL::HDS_graph_traits< CGAL::Polyhedron_3<Gt,I,HDS,A> > Base;
  typedef typename Gt::Point_3 Point;
  typedef typename Base::edge_iterator undirected_edge_iterator;
};
} // namespace CGAL
#include <CGAL/boost/graph/backward_compatibility_functions.h>

namespace boost {
  // The following functions were defined in the namespace boost
  using CGAL::vertices;
  using CGAL::edges;
  using CGAL::num_vertices;
  using CGAL::num_edges;
  using CGAL::out_edges;
  using CGAL::in_edges;
  using CGAL::target;
  using CGAL::source;
} // namespace boost

#endif //CGAL_NO_DEPRECATED_CODE

#undef CGAL_HDS_PARAM_

#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#endif // CGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYHEDRON_3_H
