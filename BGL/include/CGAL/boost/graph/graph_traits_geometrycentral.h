// Copyright (c) 2024  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_TRAITS_GEOMETRYCENTRAL_H
#define CGAL_BOOST_GRAPH_TRAITS_GEOMETRYCENTRAL_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>


#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include <CGAL/boost/graph/iterator.h>

namespace  std {
template <>
struct iterator_traits<geometrycentral::RangeIteratorBase<geometrycentral::surface::VertexRangeF>>
{
  typedef std::forward_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;
  typedef geometrycentral::surface::Vertex value_type;
  typedef value_type* pointer;
  typedef value_type & reference;
};

template <>
struct iterator_traits<geometrycentral::RangeIteratorBase<geometrycentral::surface::HalfedgeRangeF>>
{
  typedef std::forward_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;
  typedef geometrycentral::surface::Halfedge value_type;
  typedef value_type* pointer;
  typedef value_type & reference;
};

template <>
struct iterator_traits<geometrycentral::RangeIteratorBase<geometrycentral::surface::EdgeRangeF>>
{
  typedef std::forward_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;
  typedef geometrycentral::surface::Edge value_type;
  typedef value_type* pointer;
  typedef value_type & reference;
};

template <>
struct iterator_traits<geometrycentral::RangeIteratorBase<geometrycentral::surface::FaceRangeF>>
{
  typedef std::forward_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;
  typedef geometrycentral::surface::Face value_type;
  typedef value_type* pointer;
  typedef value_type & reference;
};
}

namespace boost {

template <>
struct graph_traits< geometrycentral::surface::ManifoldSurfaceMesh >
{
private:
  typedef geometrycentral::surface::ManifoldSurfaceMesh SM;

  struct SM_graph_traversal_category : public virtual boost::bidirectional_graph_tag,
                                       public virtual boost::vertex_list_graph_tag,
                                       public virtual boost::edge_list_graph_tag,
                                       public virtual boost::adjacency_graph_tag
  {};

public:
  // Graph
  typedef geometrycentral::surface::Vertex                                 vertex_descriptor;
  typedef geometrycentral::surface::Edge  edge_descriptor;
  typedef boost::undirected_tag                                            directed_category;
  typedef boost::disallow_parallel_edge_tag                                edge_parallel_category;
  typedef SM_graph_traversal_category                                      traversal_category;

  // HalfedgeGraph
  typedef geometrycentral::surface::Halfedge              halfedge_descriptor;

   // FaceGraph
  typedef geometrycentral::surface::Face   face_descriptor;

  // VertexListGraph
  typedef typename geometrycentral::RangeIteratorBase<geometrycentral::surface::VertexRangeF>   vertex_iterator;
  typedef typename std::size_t              vertices_size_type;
  // EdgeListGraph
  typedef typename geometrycentral::RangeIteratorBase<geometrycentral::surface::EdgeRangeF>     edge_iterator;

  typedef typename std::size_t              edges_size_type;
  // HalfEdgeListGraph
  typedef typename geometrycentral::RangeIteratorBase<geometrycentral::surface::HalfedgeRangeF> halfedge_iterator;
  typedef typename std::size_t              halfedges_size_type;
  // FaceListGraph
  typedef typename geometrycentral::RangeIteratorBase<geometrycentral::surface::FaceRangeF>     face_iterator;
  typedef typename std::size_t              faces_size_type;

  // IncidenceGraph
  typedef typename std::size_t              degree_size_type;


  typedef CGAL::In_edge_iterator<SM> in_edge_iterator;
  typedef CGAL::Out_edge_iterator<SM> out_edge_iterator;

  typedef CGAL::Vertex_around_target_iterator<SM> adjacency_iterator;


  // nulls
  static vertex_descriptor   null_vertex() { return vertex_descriptor(); }
  static face_descriptor     null_face()   { return face_descriptor(); }
  static halfedge_descriptor null_halfedge()   { return halfedge_descriptor(); }
};

template<>
struct graph_traits< const geometrycentral::surface::ManifoldSurfaceMesh >
  : public graph_traits< geometrycentral::surface::ManifoldSurfaceMesh >
{ };

} // namespace boost

namespace geometrycentral {

namespace surface {


// Forward declarations
typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor
halfedge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
         const geometrycentral::surface::ManifoldSurfaceMesh& );


// Declarations

typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertices_size_type
inline num_vertices(const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  return sm.nVertices();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edges_size_type
inline num_edges(const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  return sm.nEdges();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::degree_size_type
inline degree(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
              const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return v.degree();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::degree_size_type
inline degree(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::face_descriptor f,
              const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return f.degree();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::degree_size_type
inline out_degree(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
                  const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return v.degree();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::degree_size_type
inline in_degree(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
                 const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return v.degree();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor
inline source(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_descriptor e,
              const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return e.halfedge().tailVertex();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor
inline source(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
              const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return h.tailVertex();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor
inline target(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_descriptor e,
              const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  return e.halfedge().tipVertex();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor
inline target(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
              const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return h.tipVertex();
}


CGAL::Iterator_range<typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_iterator>
inline vertices(const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  geometrycentral::surface::VertexSet vs = const_cast<geometrycentral::surface::ManifoldSurfaceMesh&>(sm).vertices();
  return CGAL::make_range(vs.begin(), vs.end());
}


CGAL::Iterator_range<typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_iterator>
inline edges(const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  geometrycentral::surface::EdgeSet es = const_cast<geometrycentral::surface::ManifoldSurfaceMesh&>(sm).edges();
  return CGAL::make_range(es.begin(), es.end());
}


CGAL::Iterator_range<typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::in_edge_iterator>
inline in_edges(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
                const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  typedef typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::in_edge_iterator Iter;

  return make_range(Iter(halfedge(v,sm),sm), Iter(halfedge(v,sm),sm,1));
}



CGAL::Iterator_range<typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::out_edge_iterator>
inline out_edges(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
                 const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  typedef typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::out_edge_iterator Iter;
  return make_range(Iter(halfedge(v,sm),sm), Iter(halfedge(v,sm),sm,1));
}



CGAL::Iterator_range<typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::adjacency_iterator>
inline adjacent_vertices(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
                         const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  return CGAL::vertices_around_target(v,sm);
}


std::pair<typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_descriptor,
          bool>
inline edge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor u,
            typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
            const geometrycentral::surface::ManifoldSurfaceMesh& sm) {
  typedef typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_descriptor edge_descriptor;
  for(auto e : u.adjacentEdges()) {
    if(e.otherVertex(u) == v){
      return std::make_pair(e, true);
    }
  }
  return std::make_pair(edge_descriptor(), false);
}


//
// HalfedgeGraph
//

typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor
inline next(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
            const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return h.next();
}

/*

typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor
prev(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
     const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  return sm.prev(h);
}
*/


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor
inline opposite(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
                const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return h.twin();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_descriptor
inline edge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
            const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return h.edge();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor
inline halfedge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_descriptor e,
                const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return e.halfedge();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor
inline halfedge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
                const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return v.halfedge().twin();
}


//
// HalfedgeListGraph
//

// CGAL::Iterator_range<typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_iterator>
geometrycentral::surface::HalfedgeSet
inline halfedges(const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  return const_cast<geometrycentral::surface::ManifoldSurfaceMesh&>(sm).halfedges();
}



typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedges_size_type
inline num_halfedges(const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  return sm.nHalfedges();
}


#if 0
//
// MutableHalfedgeGraph
//
template<typename P>
void
set_next(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h1,
         typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h2,
         geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  sm.set_next(h1, h2);
}



template<typename P>
void
set_target(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
           typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
           geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  sm.set_target(h, v);
}


template<typename P>
void
set_halfedge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
             typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
             geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  sm.set_halfedge(v, h);
}


template<typename P>
typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_descriptor
add_edge(geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  return sm.edge(sm.add_edge());
}

#endif



void
inline collect_garbage(geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  sm.compress();
}


//
// FaceGraph
//
typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor
inline halfedge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::face_descriptor f,
                const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return f.halfedge();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::face_descriptor
inline face(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
            const geometrycentral::surface::ManifoldSurfaceMesh& )
{
  return h.face();
}

#if 0

//
// MutableFaceGraph
//
template<typename P>
void
set_face(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
         typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::face_descriptor f,
         geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  sm.set_face(h, f);
}


template<typename P>
void
set_halfedge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::face_descriptor f,
             typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
             geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  sm.set_halfedge(f, h);
}

#endif


//
// FaceListGraph
//

typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::faces_size_type
inline num_faces(const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  return sm.nFaces();
}


CGAL::Iterator_range<typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::face_iterator>

inline faces(const geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  geometrycentral::surface::FaceSet fs = const_cast<geometrycentral::surface::ManifoldSurfaceMesh&>(sm).faces();
  return CGAL::make_range(fs.begin(), fs.end());
}

#if 0

typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor
add_vertex(geometrycentral::surface::ManifoldSurfaceMesh& sm) {
  return sm.add_vertex();
}


typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor
add_vertex(const typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_property_type& p, geometrycentral::surface::ManifoldSurfaceMesh& sm) {
  return sm.add_vertex(p);
}

// MutableGraph
template<typename P>
void
reserve(geometrycentral::surface::ManifoldSurfaceMesh& sm,
        typename boost::graph_traits< geometrycentral::surface::ManifoldSurfaceMesh >::vertices_size_type nv,
        typename boost::graph_traits< geometrycentral::surface::ManifoldSurfaceMesh >::edges_size_type ne,
        typename boost::graph_traits< geometrycentral::surface::ManifoldSurfaceMesh >::faces_size_type nf)
{
  sm.reserve(nv, ne, nf);
}



void
remove_vertex(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
              geometrycentral::surface::ManifoldSurfaceMesh& sm) {

  sm.remove_vertex(v);
}


void
remove_edge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor u,
            typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
            geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_descriptor e = edge(u, v, sm);
  remove_edge(e,sm);
}


void
remove_edge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_descriptor e,
            geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  sm.remove_edge(e);
}


void
remove_edge(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_iterator eiter,
            geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  remove_edge(*eiter, sm);
}

template<typename P>
void
remove_face(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::face_descriptor f,
            geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  sm.remove_face(f);
}

template<typename P>
void
remove_all_elements(geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  sm.clear_without_removing_property_maps();
}

template<typename P>
typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::face_descriptor
add_face(geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  return sm.add_face();
}

template<typename P, typename InputIterator>
typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::face_descriptor
add_face(InputIterator begin, InputIterator end, geometrycentral::surface::ManifoldSurfaceMesh& sm)
{
  std::vector<typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor>
    v(begin, end);
  return sm.add_face(v);
}

template<typename P>
void normalize_border(const geometrycentral::surface::ManifoldSurfaceMesh&)
{}



bool is_valid_vertex_descriptor(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::vertex_descriptor v,
                                const geometrycentral::surface::ManifoldSurfaceMesh& g,
                                const bool verbose = false)
{
  if(!g.is_valid(v, verbose))
    return false;

  return BGL::is_valid_vertex_descriptor(v, g, verbose);
}


bool is_valid_halfedge_descriptor(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::halfedge_descriptor h,
                                  const geometrycentral::surface::ManifoldSurfaceMesh& g,
                                  const bool verbose = false)
{
  if(!g.is_valid(h, verbose))
    return false;

  return BGL::is_valid_halfedge_descriptor(h, g, verbose);
}


bool is_valid_edge_descriptor(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::edge_descriptor e,
                              const geometrycentral::surface::ManifoldSurfaceMesh& g,
                              const bool verbose = false)
{
  if(!g.is_valid(e, verbose))
    return false;

  return BGL::is_valid_edge_descriptor(e, g, verbose);
}


bool is_valid_face_descriptor(typename boost::graph_traits<geometrycentral::surface::ManifoldSurfaceMesh >::face_descriptor f,
                              const geometrycentral::surface::ManifoldSurfaceMesh& g,
                              const bool verbose = false)
{
  if(!g.is_valid(f, verbose))
    return false;

  return BGL::is_valid_face_descriptor(f, g, verbose);
}
#endif

} // namespace surface
} // namespace geometrycentral

#endif
