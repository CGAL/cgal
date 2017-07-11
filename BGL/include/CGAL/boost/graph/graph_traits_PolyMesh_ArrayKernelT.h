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
// Author(s)     : Andreas Fabri, Philipp Moeller

#ifndef CGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYMESH_ARRAYKERNELT_H
#define CGAL_BOOST_GRAPH_GRAPH_TRAITS_POLYMESH_ARRAYKERNELT_H

// include this to avoid a VC15 warning
#include <CGAL/boost/graph/named_function_params.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/iterator/transform_iterator.hpp>

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/internal/OM_iterator_from_circulator.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/assertions.h>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/hash_openmesh.h>

// http://openmesh.org/Documentation/OpenMesh-Doc-Latest/classOpenMesh_1_1Concepts_1_1KernelT.html

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4267)
#endif
namespace CGAL { namespace internal {


template<typename Halfedge_handle>
class OMesh_edge {
public:
  OMesh_edge() : halfedge_() {}
  explicit OMesh_edge(const Halfedge_handle& h) : halfedge_(h) {}
  Halfedge_handle halfedge() const { return halfedge_; }
  bool is_valid() const { return halfedge_.is_valid(); }

  bool
  operator==(const OMesh_edge& other) {
    if(halfedge_ == other.halfedge_) {
      return true;
    } else if(halfedge_ != Halfedge_handle()) {
      return opposite() == other.halfedge_;
    } else {
      return false;
    }
  }

  bool operator<(const OMesh_edge& other) const
  { 
    return this->idx() < other.idx();
  }

  bool
  operator!=(const OMesh_edge& other) { return !(*this == other); }

  Halfedge_handle
  opposite() const { return Halfedge_handle((halfedge_.idx() & 1) ? halfedge_.idx()-1 : halfedge_.idx()+1); }

  OMesh_edge
  opposite_edge() const { return OMesh_edge(Halfedge_handle((halfedge_.idx() & 1) ? halfedge_.idx()-1 : halfedge_.idx()+1)); }

  unsigned int idx() const { return halfedge_.idx() / 2; }
private:
  Halfedge_handle halfedge_;
};

template <typename Halfedge_handle, typename OMeshEdge>
struct Convert_omesh_edge
{
  typedef OMesh_edge<Halfedge_handle> result_type;
  result_type operator()(const OMeshEdge& h) const { 
    return result_type(Halfedge_handle(h.idx() * 2)); 
  }
};

template <typename Halfedge_handle>
struct Construct_omesh_edge
{
  typedef OMesh_edge<Halfedge_handle> result_type;
  template <typename T>
  result_type operator()(const T& h) const { return result_type(h); }
};

template <typename Halfedge_handle>
struct Construct_omesh_edge_opposite
{
  typedef OMesh_edge<Halfedge_handle> result_type;
  template <typename T>
  result_type operator()(const T& h) const { return result_type(h).opposite_edge(); }
};


} // internal
} // CGAL


namespace boost {

template <class K>
struct graph_traits< OpenMesh::PolyMesh_ArrayKernelT<K> >
{
private:
  typedef OpenMesh::PolyMesh_ArrayKernelT<K> SM;

  struct SM_graph_traversal_category : public virtual boost::bidirectional_graph_tag,
                                       public virtual boost::vertex_list_graph_tag,
                                       public virtual boost::edge_list_graph_tag
  {};

public:
  // Graph
  typedef typename SM::VertexHandle                                        vertex_descriptor;
  typedef typename SM::Point                                               vertex_property_type;
  typedef typename CGAL::internal::OMesh_edge<typename SM::HalfedgeHandle> edge_descriptor;
  typedef boost::undirected_tag                                            directed_category;
  typedef boost::disallow_parallel_edge_tag                                edge_parallel_category; 
  typedef SM_graph_traversal_category                                      traversal_category;

  // HalfedgeGraph
  typedef typename SM::HalfedgeHandle              halfedge_descriptor;

   // FaceGraph
  typedef typename SM::FaceHandle   face_descriptor;
  
  // VertexListGraph
  typedef typename SM::VertexIter   vertex_iterator;
  typedef unsigned int              vertices_size_type;
  // EdgeListGraph
  typedef boost::transform_iterator<
    CGAL::internal::Convert_omesh_edge<typename SM::HalfedgeHandle, typename SM::EdgeHandle>,
    typename SM::EdgeIter,
    edge_descriptor>                edge_iterator;

  typedef unsigned int              edges_size_type;
  // HalfEdgeListGraph
  typedef typename SM::HalfedgeIter halfedge_iterator;
  typedef unsigned int              halfedges_size_type;
  // FaceListGraph
  typedef typename SM::FaceIter     face_iterator;
  typedef unsigned int              faces_size_type;

  // IncidenceGraph
  typedef unsigned int              degree_size_type;

  
  typedef CGAL::In_edge_iterator<SM> in_edge_iterator;

  typedef CGAL::Out_edge_iterator<SM> out_edge_iterator;

  // nulls
  static vertex_descriptor   null_vertex() { return vertex_descriptor(); }
  static face_descriptor     null_face()   { return face_descriptor(); }
  static halfedge_descriptor     null_halfedge()   { return halfedge_descriptor(); }
};

template<typename K>
struct graph_traits< const OpenMesh::PolyMesh_ArrayKernelT<K> >  
  : public graph_traits< OpenMesh::PolyMesh_ArrayKernelT<K> >
{ };

} // namespace boost

namespace OpenMesh {

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertices_size_type
num_vertices(const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.n_vertices();
}
  

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edges_size_type
num_edges(const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.n_edges();
}
  

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::degree_size_type
degree(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v,
       const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.valence(v);
}


template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::degree_size_type
degree(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor f,
       const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.valence(f);
}

         
template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::degree_size_type
out_degree(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v,
           const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.valence(v);
}
             
  
template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::degree_size_type
in_degree(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v,
          const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.valence(v);
}
            
  
template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor
source(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor e,
       const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.from_vertex_handle(e.halfedge());
}

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor
source(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
       const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.from_vertex_handle(h);
}
         

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor
target(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor e,
       const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.to_vertex_handle(e.halfedge());
}

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor
target(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
       const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.to_vertex_handle(h);
}
    
template <typename K>
CGAL::Iterator_range<typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_iterator>
vertices(const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return CGAL::make_range(sm.vertices_sbegin(), sm.vertices_end()); 
}

 
template <typename K>
CGAL::Iterator_range<typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_iterator>
edges(const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  typedef typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_iterator iterator;
  iterator beg(sm.edges_sbegin());
  iterator end(sm.edges_end());
  return CGAL::make_range(beg,end); 
}

  
template <typename K>
CGAL::Iterator_range<typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::in_edge_iterator>
in_edges(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v,
         const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  typedef typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::in_edge_iterator Iter;

  return CGAL::make_range(Iter(halfedge(v,sm),sm), Iter(halfedge(v,sm),sm,1));
}


template <typename K>
CGAL::Iterator_range<typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::out_edge_iterator>
out_edges(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v,
          const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  typedef typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::out_edge_iterator Iter;
  return CGAL::make_range(Iter(halfedge(v,sm),sm), Iter(halfedge(v,sm),sm,1));
}


template<typename K>
std::pair<typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor,
          bool>
edge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor u, 
     typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v, 
     const OpenMesh::PolyMesh_ArrayKernelT<K>& sm) {
  typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor 
    he(sm.find_halfedge(u, v));
  return std::make_pair(he, he.is_valid());
}


//
// HalfedgeGraph
//
template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor
next(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
     const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.next_halfedge_handle(h);
}

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor
prev(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
     const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.prev_halfedge_handle(h);
}

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor
opposite(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
         const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.opposite_halfedge_handle(h);
}

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor
edge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
     const OpenMesh::PolyMesh_ArrayKernelT<K>& /*sm*/)
{
  return typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor(h);
}

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor
halfedge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor e,
         const OpenMesh::PolyMesh_ArrayKernelT<K>&)
{
  return e.halfedge();
}

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor
halfedge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v,
         const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  if(sm.halfedge_handle(v) == boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::null_halfedge()){
    return boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::null_halfedge();
  }
  // prev because OpenMesh stores out-going halfedges
  // return sm.prev_halfedge_handle(sm.halfedge_handle(v));
  return sm.opposite_halfedge_handle(sm.halfedge_handle(v));
}


template <typename K>
std::pair<
  typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor,
  bool
>
halfedge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor u,
         typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v,
         const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h = sm.find_halfedge(u, v);
  return std::make_pair(h, h.is_valid());
}



//
// HalfedgeListGraph
//
template <typename K>
CGAL::Iterator_range<typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_iterator>
halfedges(const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return CGAL::make_range(sm.halfedges_sbegin(), sm.halfedges_end());
}

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedges_size_type
num_halfedges(const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.n_halfedges();
}



//
// MutableHalfedgeGraph
//
template<typename K>
void
set_next(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h1, 
         typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h2,
         OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  sm.set_next_halfedge_handle(h1, h2);
}



template<typename K>
void
set_target(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
           typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v,
           OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  sm.set_vertex_handle(h, v);
}


template<typename K>
void
set_halfedge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v,
             typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
             OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  sm.set_halfedge_handle(v, sm.opposite_halfedge_handle(h));
}


template<typename K>
void
adjust_border_halfedge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v,
                       OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  sm.adjust_outgoing_halfedge(v);
}

template<typename K>
void
garbage_collection(OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  sm.garbage_collection();
}

template<typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor
add_edge(OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return edge(sm.new_edge(boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::null_vertex(),
                          boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::null_vertex() ), sm);
}


//
// FaceGraph
//
template<typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor
halfedge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor f,
     const OpenMesh::PolyMesh_ArrayKernelT<K>& sm) 
{
  return sm.halfedge_handle(f);
}
  
template<typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor
face(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
     const OpenMesh::PolyMesh_ArrayKernelT<K>& sm) 
{
  return sm.face_handle(h);
}



//
// MutableFaceGraph
//
template<typename K>
void
set_face(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
         typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor f,
         OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  sm.set_face_handle(h, f);
}

  
template<typename K>
void
set_halfedge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor f,
             typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor h,
             OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  sm.set_halfedge_handle(f, h);
}

template<typename K>
void
reserve(OpenMesh::PolyMesh_ArrayKernelT<K>& tm,
        typename boost::graph_traits< OpenMesh::PolyMesh_ArrayKernelT<K> >::vertices_size_type nv,
        typename boost::graph_traits< OpenMesh::PolyMesh_ArrayKernelT<K> >::edges_size_type ne,
        typename boost::graph_traits< OpenMesh::PolyMesh_ArrayKernelT<K> >::faces_size_type nf)
{
  tm.reserve(nv, ne, nf);
}

//
// FaceListGraph
//
template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::faces_size_type
num_faces(const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.n_faces();
}
  
template <typename K>
CGAL::Iterator_range<typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_iterator>
faces(const OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return CGAL::make_range(sm.faces_sbegin(), sm.faces_end()); 
}
 

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor
add_vertex(OpenMesh::PolyMesh_ArrayKernelT<K>& sm) {
  return sm.new_vertex();
}

  /*

// MutableGraph
// add a vertex with a default constructed property
template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor
add_vertex(OpenMesh::PolyMesh_ArrayKernelT<K>& sm) {
  return sm.add_vertex(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_property_type());
}

template <typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor
add_vertex(const typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_property_type& p, OpenMesh::PolyMesh_ArrayKernelT<K>& sm) {
  return sm.add_vertex(p);
}

template <typename K>
void
clear_vertex(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor, 
             OpenMesh::PolyMesh_ArrayKernelT<K>&) {
  CGAL_assert(false);
}

  */

template <typename K>
void
remove_vertex(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v, 
              OpenMesh::PolyMesh_ArrayKernelT<K>& sm) {

  sm.request_face_status();
  sm.request_vertex_status();
  sm.request_halfedge_status();
  sm.set_halfedge_handle(v, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor());
  sm.status(v).set_deleted(true);
}

  
template <typename K>
void
remove_edge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor u, 
            typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v, 
            OpenMesh::PolyMesh_ArrayKernelT<K>& sm) 
{
  typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor e = edge(u, v, sm);
  remove_edge(e,sm);
}

template <typename K>
void
remove_edge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor e, 
            OpenMesh::PolyMesh_ArrayKernelT<K>& sm) 
{
  sm.request_face_status();
  sm.request_vertex_status();
  sm.request_halfedge_status();
  sm.request_edge_status();

  typedef typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h1 = halfedge(e,sm);
  halfedge_descriptor h2 = opposite(halfedge(e,sm),sm);
  sm.status(sm.edge_handle(h1)).set_deleted(true);
  sm.status(h1).set_deleted(true);
  sm.status(h2).set_deleted(true);

}


template <typename K>
void
remove_edge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_iterator eiter, 
            OpenMesh::PolyMesh_ArrayKernelT<K>& sm) 
{
  remove_edge(*eiter, sm);
}

template<typename K>
void
remove_face(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor f, 
            OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{

  sm.request_face_status();
  sm.request_vertex_status();
  sm.request_halfedge_status();
  
  set_halfedge(f, typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::halfedge_descriptor(), sm);
  sm.status(f).set_deleted(true);
}

#if 0 // conflits with function in Euler_operations.h
template<typename K>
std::pair<typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::edge_descriptor,
          bool>
add_edge(typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v1, 
         typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor v2, 
         OpenMesh::PolyMesh_ArrayKernelT<K>& sm) {
  
  return sm.new_edge(v1, v2);
}
#endif

template<typename K>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor
add_face(OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return sm.new_face();
}

template<typename K, typename InputIterator>
typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::face_descriptor
add_face(InputIterator begin, InputIterator end, OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  std::vector<typename boost::graph_traits<OpenMesh::PolyMesh_ArrayKernelT<K> >::vertex_descriptor> 
    v(begin, end);
  return sm.add_face(v);
}

template<typename K>
bool is_valid(OpenMesh::PolyMesh_ArrayKernelT<K>& sm, bool /* verbose */ = false)
{
  return CGAL::is_valid_polygon_mesh(sm);
}

} // namespace OpenMesh

namespace CGAL {

// Overload CGAL::clear function. PolyMesh_ArrayKernel behaves
// differently from other meshes. Calling clear does not affect the
// number of vertices, edges, or faces in the mesh. To get actual
// numbers it is necessary to first collect garbage. We add an
// overlaod to get consistent behavior.
template<typename K>
void clear(OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  sm.clear();
  sm.garbage_collection(true, true, true);
  CGAL_postcondition(num_edges(sm) == 0);
  CGAL_postcondition(num_vertices(sm) == 0);
  CGAL_postcondition(num_faces(sm) == 0);
}


template<typename K>
bool read_off(std::istream& is, OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return OpenMesh::IO::read_mesh(sm, is, ".OFF");
}


template<typename K>
bool write_off(std::ostream& os, OpenMesh::PolyMesh_ArrayKernelT<K>& sm)
{
  return OpenMesh::IO::write_mesh(sm, os, ".OFF");
}

}
#ifndef CGAL_NO_DEPRECATED_CODE
#include <CGAL/boost/graph/backward_compatibility_functions.h>

namespace boost {
  // The following functions were defined in the namespace boost
  using OpenMesh::vertices;
  using OpenMesh::edges;
  using OpenMesh::num_vertices;
  using OpenMesh::num_edges;
  using OpenMesh::out_edges;
  using OpenMesh::in_edges;
  using OpenMesh::target;
  using OpenMesh::source;
} // namespace boost
#endif //CGAL_NO_DEPRECATED_CODE

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_BOOST_GRAPH_TRAITS_POLYMESH_ARRAYKERNELT_H
