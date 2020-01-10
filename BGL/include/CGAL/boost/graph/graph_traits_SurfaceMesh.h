// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
// 
//
// Author(s)     : Andreas Fabri, Philipp Moeller

#ifndef CGAL_BOOST_GRAPH_TRAITS_SURFACEMESH_H
#define CGAL_BOOST_GRAPH_TRAITS_SURFACEMESH_H

#ifndef DOXYGEN_RUNNING

#include <CGAL/license/Surface_mesh.h>

// include this to avoid a VC15 warning
#include <CGAL/boost/graph/Named_function_parameters.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <CGAL/boost/graph/iterator.h>

#include <pmp/SurfaceMesh.h>
#include <CGAL/assertions.h>



namespace CGAL { namespace internal {



class PMP_edge {
public:
  PMP_edge() : halfedge_() {}
  explicit PMP_edge(const pmp::Halfedge& h) : halfedge_(h) {}
  pmp::Halfedge halfedge() const { return halfedge_; }
  bool is_valid() const { return halfedge_.is_valid(); }

  bool
  operator==(const PMP_edge& other) const {
    if(halfedge_ == other.halfedge_) {
      return true;
    } else if(halfedge_ != pmp::Halfedge()) {
      return opposite() == other.halfedge_;
    } else {
      return false;
    }
  }

  bool operator<(const PMP_edge& other) const
  {
    return halfedge_.idx() < other.halfedge_.idx();
  }

  bool
  operator!=(const PMP_edge& other) { return !(*this == other); }

  pmp::Halfedge
  opposite() const { return pmp::Halfedge((halfedge_.idx() & 1) ? halfedge_.idx()-1 : halfedge_.idx()+1); }

  PMP_edge
  opposite_edge() const { return PMP_edge(pmp::Halfedge((halfedge_.idx() & 1) ? halfedge_.idx()-1 : halfedge_.idx()+1)); }

  unsigned int idx() const { return halfedge_.idx() / 2; }

  
friend
static std::ostream& operator<<(std::ostream& os, const PMP_edge& e)
{
  os << "Edge around " << e.halfedge_ ; 
  return os;
}

  
private:
  pmp::Halfedge halfedge_;
};
    
template <typename PMPEdge>
struct Convert_pmp_edge
{
  typedef PMP_edge result_type;
  result_type operator()(const PMPEdge& h) const {
    return result_type(pmp::Halfedge(h.idx() * 2));
  }
};


struct Construct_pmp_edge
{
  typedef PMP_edge result_type;
  template <typename T>
  result_type operator()(const T& h) const { return result_type(h); }
};

struct Construct_pmp_edge_opposite
{
  typedef PMP_edge result_type;
  template <typename T>
  result_type operator()(const T& h) const { return result_type(h).opposite_edge(); }
};

  } } // CGAL::internal

namespace boost {

template <>
struct graph_traits< pmp::SurfaceMesh >
{
private:
  typedef pmp::SurfaceMesh SM;

  struct SM_graph_traversal_category : public virtual boost::bidirectional_graph_tag,
                                       public virtual boost::vertex_list_graph_tag,
                                       public virtual boost::edge_list_graph_tag
  {};

public:
  // Graph
  typedef pmp::Vertex                                        vertex_descriptor;
  typedef pmp::Point                                         vertex_property_type;
  typedef typename CGAL::internal::PMP_edge                  edge_descriptor;
  typedef boost::undirected_tag                              directed_category;
  typedef boost::disallow_parallel_edge_tag                  edge_parallel_category; 
  typedef SM_graph_traversal_category                        traversal_category;

  // HalfedgeGraph
  typedef pmp::Halfedge                                       halfedge_descriptor;

   // FaceGraph
  typedef pmp::Face                                           face_descriptor;
  
  // VertexListGraph
  typedef typename SM::VertexIterator                         vertex_iterator;
  typedef typename std::size_t                                vertices_size_type;

  // EdgeListGraph
  typedef boost::transform_iterator<
    CGAL::internal::Convert_pmp_edge<typename pmp::Edge>,
    typename SM::EdgeIterator,
    edge_descriptor>                                          edge_iterator;

  typedef typename std::size_t                                edges_size_type;
  // HalfEdgeListGraph
  typedef typename SM::HalfedgeIterator                       halfedge_iterator;
  typedef typename std::size_t                                halfedges_size_type;
  // FaceListGraph
  typedef typename SM::FaceIterator                           face_iterator;
  typedef typename std::size_t                                faces_size_type;

  // IncidenceGraph
  typedef typename std::size_t                                 degree_size_type;

  
  typedef CGAL::In_edge_iterator<SM> in_edge_iterator;

  typedef CGAL::Out_edge_iterator<SM> out_edge_iterator;

  // nulls
  static vertex_descriptor   null_vertex() { return vertex_descriptor(); }
  static face_descriptor     null_face()   { return face_descriptor(); }
  static halfedge_descriptor null_halfedge()   { return halfedge_descriptor(); }
};

template <>
struct graph_traits< const pmp::SurfaceMesh >  
  : public graph_traits< pmp::SurfaceMesh >
{ };

} // namespace boost

namespace pmp {


typename boost::graph_traits<pmp::SurfaceMesh >::vertices_size_type
num_vertices(const pmp::SurfaceMesh& sm)
{
  return sm.vertices_size();
}
  


typename boost::graph_traits<pmp::SurfaceMesh >::edges_size_type
num_edges(const pmp::SurfaceMesh& sm)
{
  return sm.edges_size();
}
  


typename boost::graph_traits<pmp::SurfaceMesh >::degree_size_type
degree(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v,
       const pmp::SurfaceMesh& sm)
{
  return sm.valence(v);
}



typename boost::graph_traits<pmp::SurfaceMesh >::degree_size_type
degree(typename boost::graph_traits<pmp::SurfaceMesh >::face_descriptor f,
       const pmp::SurfaceMesh& sm)
{
  return sm.valence(f);
}

         

typename boost::graph_traits<pmp::SurfaceMesh >::degree_size_type
out_degree(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v,
           const pmp::SurfaceMesh& sm)
{
  return sm.valence(v);
}
             
  

typename boost::graph_traits<pmp::SurfaceMesh >::degree_size_type
in_degree(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v,
          const pmp::SurfaceMesh& sm)
{
  return sm.valence(v);
}
            
  

typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor
source(typename boost::graph_traits<pmp::SurfaceMesh >::edge_descriptor e,
       const pmp::SurfaceMesh& sm)
{
  return sm.from_vertex(e.halfedge());
}


typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor
source(typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
       const pmp::SurfaceMesh& sm)
{
  return sm.from_vertex(h);
}
         


typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor
target(typename boost::graph_traits<pmp::SurfaceMesh >::edge_descriptor e,
       const pmp::SurfaceMesh& sm)
{
  return sm.to_vertex(e.halfedge());
}


typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor
target(typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
       const pmp::SurfaceMesh& sm)
{
  return sm.to_vertex(h);
}
    

CGAL::Iterator_range<typename boost::graph_traits<pmp::SurfaceMesh >::vertex_iterator>
vertices(const pmp::SurfaceMesh& sm)
{
  return CGAL::make_range(sm.vertices_begin(), sm.vertices_end()); 
}

 

CGAL::Iterator_range<typename boost::graph_traits<pmp::SurfaceMesh >::edge_iterator>
edges(const pmp::SurfaceMesh& sm)
{
  typedef typename boost::graph_traits<pmp::SurfaceMesh>::edge_iterator iterator;
  iterator beg(sm.edges_begin());
  iterator end(sm.edges_end());
  return CGAL::make_range(beg,end); 
}


  // fwd decl
boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor
halfedge(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v,
         const pmp::SurfaceMesh& sm);
 
inline
CGAL::Iterator_range<typename boost::graph_traits<pmp::SurfaceMesh >::in_edge_iterator>
in_edges(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v,
         const pmp::SurfaceMesh& sm)
{
  typedef typename boost::graph_traits<pmp::SurfaceMesh >::in_edge_iterator Iter;

  return CGAL::make_range(Iter(halfedge(v,sm),sm), Iter(halfedge(v,sm),sm,1));
}



CGAL::Iterator_range<typename boost::graph_traits<pmp::SurfaceMesh >::out_edge_iterator>
out_edges(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v,
          const pmp::SurfaceMesh& sm)
{
  typedef typename boost::graph_traits<pmp::SurfaceMesh >::out_edge_iterator Iter;
  return CGAL::make_range(Iter(halfedge(v,sm),sm), Iter(halfedge(v,sm),sm,1));
}


template<typename P>
std::pair<typename boost::graph_traits<pmp::SurfaceMesh >::edge_descriptor,
          bool>
edge(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor u, 
     typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v, 
     const pmp::SurfaceMesh& sm) {
  typename boost::graph_traits<pmp::SurfaceMesh >::edge_descriptor 
    he(sm.halfedge(u, v));
  return std::make_pair(he, he.is_valid());
}


//
// HalfedgeGraph
//

typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor
next(typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
     const pmp::SurfaceMesh& sm)
{
  return sm.next_halfedge(h);
}


typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor
prev(typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
     const pmp::SurfaceMesh& sm)
{
  return sm.prev_halfedge(h);
}


typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor
opposite(typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
         const pmp::SurfaceMesh& sm)
{
  return sm.opposite_halfedge(h);
}


typename boost::graph_traits<pmp::SurfaceMesh >::edge_descriptor
edge(typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
     const pmp::SurfaceMesh& sm)
{
  return typename boost::graph_traits<pmp::SurfaceMesh >::edge_descriptor(h);
}


typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor
halfedge(typename boost::graph_traits<pmp::SurfaceMesh >::edge_descriptor e,
         const pmp::SurfaceMesh& sm)
{
  return e.halfedge();
}

inline
typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor
halfedge(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v,
         const pmp::SurfaceMesh& sm)
{
  if(sm.halfedge(v) == boost::graph_traits<pmp::SurfaceMesh>::null_halfedge()){
    return boost::graph_traits<pmp::SurfaceMesh>::null_halfedge();
  }
  return sm.opposite_halfedge(sm.halfedge(v));
}

#if 0
  

std::pair<
  typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor,
  bool
>
halfedge(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor u,
         typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v,
         const pmp::SurfaceMesh& sm)
{
  typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h = sm.halfedge(u, v);
  return std::make_pair(h, h.is_valid());
}

#endif

//
// HalfedgeListGraph
//

CGAL::Iterator_range<typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_iterator>
halfedges(const pmp::SurfaceMesh& sm)
{
  return CGAL::make_range(sm.halfedges_begin(), sm.halfedges_end());
}



typename boost::graph_traits<pmp::SurfaceMesh >::halfedges_size_type
num_halfedges(const pmp::SurfaceMesh& sm)
{
  return sm.halfedges_size();
}



//
// MutableHalfedgeGraph
//
template<typename P>
void
set_next_halfedge(typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h1, 
                  typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h2,
                  pmp::SurfaceMesh& sm)
{
  sm.set_next(h1, h2);
}



template<typename P>
void
set_target(typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
           typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v,
           pmp::SurfaceMesh& sm)
{
  sm.set_vertex(h, v);
}


template<typename P>
void
set_halfedge(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v,
             typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
             pmp::SurfaceMesh& sm)
{
  sm.set_halfedge(v, sm.opposite_halfedge(h));
}


template<typename P>
void
collect_garbage(pmp::SurfaceMesh& sm)
{
  sm.garbage_collection();
}

template<typename P>
typename boost::graph_traits<pmp::SurfaceMesh >::edge_descriptor
add_edge(pmp::SurfaceMesh& sm)
{
  return sm.edge(sm.add_edge());
}


//
// FaceGraph
//
template<typename P>
typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor
halfedge(typename boost::graph_traits<pmp::SurfaceMesh >::face_descriptor f,
     const pmp::SurfaceMesh& sm) 
{
  return sm.halfedge(f);
}
  
template<typename P>
typename boost::graph_traits<pmp::SurfaceMesh >::face_descriptor
face(typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
     const pmp::SurfaceMesh& sm) 
{
  return sm.face(h);
}



//
// MutableFaceGraph
//
template<typename P>
void
set_face(typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
         typename boost::graph_traits<pmp::SurfaceMesh >::face_descriptor f,
         pmp::SurfaceMesh& sm)
{
  sm.set_face(h, f);
}

  
template<typename P>
void
set_halfedge(typename boost::graph_traits<pmp::SurfaceMesh >::face_descriptor f,
             typename boost::graph_traits<pmp::SurfaceMesh >::halfedge_descriptor h,
             pmp::SurfaceMesh& sm)
{
  sm.set_halfedge(f, h);
}

 
//
// FaceListGraph
//

typename boost::graph_traits<pmp::SurfaceMesh >::faces_size_type
num_faces(const pmp::SurfaceMesh& sm)
{
  return sm.faces_size();
}


CGAL::Iterator_range<typename boost::graph_traits<pmp::SurfaceMesh >::face_iterator>
faces(const pmp::SurfaceMesh& sm)
{
  return CGAL::make_range(sm.faces_begin(), sm.faces_end()); 
}
 


typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor
add_vertex(pmp::SurfaceMesh& sm) {
    typename boost::graph_traits<pmp::SurfaceMesh >::vertex_property_type p;
  return sm.add_vertex(p);
}


typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor
add_vertex(const typename boost::graph_traits<pmp::SurfaceMesh >::vertex_property_type& p, pmp::SurfaceMesh& sm) {
  return sm.add_vertex(p);
}

// MutableGraph
template<typename P>
void
reserve(pmp::SurfaceMesh& sm,
        typename boost::graph_traits< pmp::SurfaceMesh >::vertices_size_type nv,
        typename boost::graph_traits< pmp::SurfaceMesh >::edges_size_type ne,
        typename boost::graph_traits< pmp::SurfaceMesh >::faces_size_type nf)
{
  sm.reserve(nv, ne, nf);
}



void
remove_vertex(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v, 
              pmp::SurfaceMesh& sm) {

  sm.delete_vertex(v);
}

  
#if 0
void
remove_edge(typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor u, 
            typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor v, 
            pmp::SurfaceMesh& sm) 
{
  typename boost::graph_traits<pmp::SurfaceMesh >::edge_descriptor e = edge(u, v, sm);
  sm.delete_edge(e);
}
#endif

void
remove_edge(typename boost::graph_traits<pmp::SurfaceMesh >::edge_descriptor e, 
            pmp::SurfaceMesh& sm) 
{
  sm.delete_edge(sm.edge(halfedge(e,sm)));
}



void
remove_edge(typename boost::graph_traits<pmp::SurfaceMesh >::edge_iterator eiter, 
            pmp::SurfaceMesh& sm) 
{
  remove_edge(*eiter,sm);
}


void
remove_face(typename boost::graph_traits<pmp::SurfaceMesh >::face_descriptor f, 
            pmp::SurfaceMesh& sm)
{
  sm.delete_face(f);
}


typename boost::graph_traits<pmp::SurfaceMesh >::face_descriptor
add_face(pmp::SurfaceMesh& sm)
{
    assert(false); // 
    std::vector<typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor> v;
        
  return sm.add_face(v);
}

  
template<typename P, typename InputIterator>
typename boost::graph_traits<pmp::SurfaceMesh >::face_descriptor
add_face(InputIterator begin, InputIterator end, pmp::SurfaceMesh& sm)
{
  std::vector<typename boost::graph_traits<pmp::SurfaceMesh >::vertex_descriptor> 
    v(begin, end);
  return sm.add_face(v);
}


void normalize_border(const pmp::SurfaceMesh&)
{}

} // namespace pmp

#endif // DOXYGEN_RUNNING

#endif // CGAL_BOOST_GRAPH_TRAITS_SURFACEMESH_H
