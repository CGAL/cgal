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

#include <pmp/SurfaceMesh.h>

namespace std {

template<>
struct iterator_traits<pmp::SurfaceMesh::VertexIterator> {
  typedef int difference_type;
  typedef pmp::Vertex value_type;
  typedef const pmp::Vertex& reference;
  typedef pmp::Vertex* pointer;
  typedef std::random_access_iterator_tag iterator_category;
};

template<>
struct iterator_traits<pmp::SurfaceMesh::EdgeIterator> {
  typedef int difference_type;
  typedef pmp::Edge value_type;
  typedef const pmp::Edge& reference;
  typedef pmp::Edge* pointer;
  typedef std::random_access_iterator_tag iterator_category;
};
template<>
struct iterator_traits<pmp::SurfaceMesh::HalfedgeIterator> {
  typedef int difference_type;
  typedef pmp::Halfedge value_type;
  typedef const pmp::Halfedge& reference;
  typedef pmp::Halfedge* pointer;
  typedef std::random_access_iterator_tag iterator_category;
};
template<>
struct iterator_traits<pmp::SurfaceMesh::FaceIterator> {
  typedef int difference_type;
  typedef pmp::Face value_type;
  typedef const pmp::Face& reference;
  typedef pmp::Face* pointer;
  typedef std::random_access_iterator_tag iterator_category;
};
}

namespace CGAL{
template <typename T>
void
remove_property(pmp::VertexProperty<T> pm, pmp::SurfaceMesh& sm);
}

// include this to avoid a VC15 warning
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <CGAL/boost/graph/iterator.h>

#include <CGAL/assertions.h>


namespace pmp {

  int operator-(pmp::SurfaceMesh::FaceIterator a, pmp::SurfaceMesh::FaceIterator b)
  {
    return (*a).idx() - (*b).idx();
  }

  int operator-(pmp::SurfaceMesh::HalfedgeIterator a, pmp::SurfaceMesh::HalfedgeIterator b)
  {
    return (*a).idx() - (*b).idx();
  }

  int operator-(pmp::SurfaceMesh::VertexIterator a, pmp::SurfaceMesh::VertexIterator b)
  {
    return (*a).idx() - (*b).idx();
  }
}

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
std::ostream& operator<<(std::ostream& os, const PMP_edge& e)
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
  //typedef pmp::Point                                         vertex_property_type;
  //typedef boost::edge_weight_t                               edge_property_type;
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

inline
boost::graph_traits<SurfaceMesh>::vertices_size_type
num_vertices(const SurfaceMesh& sm)
{
  return sm.vertices_size();
}


inline
boost::graph_traits<SurfaceMesh>::edges_size_type
num_edges(const SurfaceMesh& sm)
{
  return sm.edges_size();
}


inline
boost::graph_traits<SurfaceMesh>::degree_size_type
degree(boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
       const SurfaceMesh& sm)
{
  return sm.valence(v);
}


inline
boost::graph_traits<SurfaceMesh>::degree_size_type
degree(boost::graph_traits<SurfaceMesh>::face_descriptor f,
       const SurfaceMesh& sm)
{
  return sm.valence(f);
}


inline
boost::graph_traits<SurfaceMesh>::degree_size_type
out_degree(boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
           const SurfaceMesh& sm)
{
  return sm.valence(v);
}


inline
boost::graph_traits<SurfaceMesh>::degree_size_type
in_degree(boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
          const SurfaceMesh& sm)
{
  return sm.valence(v);
}


inline
boost::graph_traits<SurfaceMesh>::vertex_descriptor
source(boost::graph_traits<SurfaceMesh>::edge_descriptor e,
       const SurfaceMesh& sm)
{
  return sm.from_vertex(e.halfedge());
}


inline
boost::graph_traits<SurfaceMesh>::vertex_descriptor
source(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
       const SurfaceMesh& sm)
{
  return sm.from_vertex(h);
}


inline
boost::graph_traits<SurfaceMesh>::vertex_descriptor
target(boost::graph_traits<SurfaceMesh>::edge_descriptor e,
       const SurfaceMesh& sm)
{
  return sm.to_vertex(e.halfedge());
}


inline
boost::graph_traits<SurfaceMesh>::vertex_descriptor
target(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
       const SurfaceMesh& sm)
{
  return sm.to_vertex(h);
}


inline
CGAL::Iterator_range<boost::graph_traits<SurfaceMesh>::vertex_iterator>
vertices(const SurfaceMesh& sm)
{
  return CGAL::make_range(sm.vertices_begin(), sm.vertices_end());
}


inline
CGAL::Iterator_range<boost::graph_traits<SurfaceMesh>::edge_iterator>
edges(const SurfaceMesh& sm)
{
  typedef boost::graph_traits<SurfaceMesh>::edge_iterator iterator;
  iterator beg(sm.edges_begin());
  iterator end(sm.edges_end());
  return CGAL::make_range(beg,end);
}


  // fwd decl
boost::graph_traits<SurfaceMesh>::halfedge_descriptor
halfedge(boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
         const SurfaceMesh& sm);

inline
CGAL::Iterator_range<boost::graph_traits<SurfaceMesh>::in_edge_iterator>
in_edges(boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
         const SurfaceMesh& sm)
{
  typedef boost::graph_traits<SurfaceMesh>::in_edge_iterator Iter;

  return CGAL::make_range(Iter(halfedge(v,sm),sm), Iter(halfedge(v,sm),sm,1));
}


inline
CGAL::Iterator_range<boost::graph_traits<SurfaceMesh>::out_edge_iterator>
out_edges(boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
          const SurfaceMesh& sm)
{
  typedef boost::graph_traits<SurfaceMesh>::out_edge_iterator Iter;
  return CGAL::make_range(Iter(halfedge(v,sm),sm), Iter(halfedge(v,sm),sm,1));
}



std::pair<
  boost::graph_traits<SurfaceMesh>::halfedge_descriptor,
  bool
>
halfedge(boost::graph_traits<SurfaceMesh>::vertex_descriptor u,
         boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
         const SurfaceMesh& sm);


boost::graph_traits<SurfaceMesh>::edge_descriptor
edge(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
     const SurfaceMesh& sm);


inline
std::pair<boost::graph_traits<SurfaceMesh>::edge_descriptor,
          bool>
edge(boost::graph_traits<SurfaceMesh>::vertex_descriptor u,
     boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
     const SurfaceMesh& sm) {

  std::pair<boost::graph_traits<SurfaceMesh>::halfedge_descriptor,bool> res = halfedge(u,v,sm);

  return std::make_pair(edge(res.first,sm), res.second);
}


//
// HalfedgeGraph
//
inline
boost::graph_traits<SurfaceMesh>::halfedge_descriptor
next(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
     const SurfaceMesh& sm)
{
  return sm.next_halfedge(h);
}


inline
boost::graph_traits<SurfaceMesh>::halfedge_descriptor
prev(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
     const SurfaceMesh& sm)
{
  return sm.prev_halfedge(h);
}


inline
boost::graph_traits<SurfaceMesh>::halfedge_descriptor
opposite(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
         const SurfaceMesh& sm)
{
  return sm.opposite_halfedge(h);
}


inline
boost::graph_traits<SurfaceMesh>::edge_descriptor
edge(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
     const SurfaceMesh& sm)
{
  return boost::graph_traits<SurfaceMesh>::edge_descriptor(h);
}


inline
boost::graph_traits<SurfaceMesh>::halfedge_descriptor
halfedge(boost::graph_traits<SurfaceMesh>::edge_descriptor e,
         const SurfaceMesh& sm)
{
  return e.halfedge();
}

inline
boost::graph_traits<SurfaceMesh>::halfedge_descriptor
halfedge(boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
         const SurfaceMesh& sm)
{
  if(sm.halfedge(v) == boost::graph_traits<SurfaceMesh>::null_halfedge()){
    return boost::graph_traits<SurfaceMesh>::null_halfedge();
  }
  return sm.opposite_halfedge(sm.halfedge(v));
}


inline
std::pair<
  boost::graph_traits<SurfaceMesh>::halfedge_descriptor,
  bool
>
halfedge(boost::graph_traits<SurfaceMesh>::vertex_descriptor u,
         boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
         const SurfaceMesh& sm)
{
  SurfaceMesh::HalfedgeAroundVertexCirculator c(&sm,u), done;
  if(c){
    do {
      if(target(*c,sm)==v){
        return std::make_pair(*c,true);
      }
      ++c;
    }while(c != done);
  }
  return std::make_pair(boost::graph_traits<SurfaceMesh>::halfedge_descriptor(),false);
}


//
// HalfedgeListGraph
//
inline
CGAL::Iterator_range<boost::graph_traits<SurfaceMesh>::halfedge_iterator>
halfedges(const SurfaceMesh& sm)
{
  return CGAL::make_range(sm.halfedges_begin(), sm.halfedges_end());
}


inline
boost::graph_traits<SurfaceMesh>::halfedges_size_type
num_halfedges(const SurfaceMesh& sm)
{
  return sm.halfedges_size();
}



//
// MutableHalfedgeGraph
//

void
set_next(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h1,
         boost::graph_traits<SurfaceMesh>::halfedge_descriptor h2,
         SurfaceMesh& sm)
{
  sm.set_next_halfedge(h1, h2);
}


inline
void
set_target(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
           boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
           SurfaceMesh& sm)
{
  sm.set_vertex(h, v);
}


inline
void
set_halfedge(boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
             boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
             SurfaceMesh& sm)
{
  sm.set_halfedge(v, sm.opposite_halfedge(h));
}


inline
void
collect_garbage(SurfaceMesh& sm)
{
  sm.garbage_collection();
}

inline
boost::graph_traits<SurfaceMesh>::edge_descriptor
add_edge(SurfaceMesh& sm)
{
  return boost::graph_traits<SurfaceMesh>::edge_descriptor(sm.add_edge());
}


//
// FaceGraph
//
inline
boost::graph_traits<SurfaceMesh>::halfedge_descriptor
halfedge(boost::graph_traits<SurfaceMesh>::face_descriptor f,
     const SurfaceMesh& sm)
{
  return sm.halfedge(f);
}


inline
boost::graph_traits<SurfaceMesh>::face_descriptor
face(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
     const SurfaceMesh& sm)
{
  return sm.face(h);
}



//
// MutableFaceGraph
//
inline
void
set_face(boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
         boost::graph_traits<SurfaceMesh>::face_descriptor f,
         SurfaceMesh& sm)
{
  sm.set_face(h, f);
}


inline
void
set_halfedge(boost::graph_traits<SurfaceMesh>::face_descriptor f,
             boost::graph_traits<SurfaceMesh>::halfedge_descriptor h,
             SurfaceMesh& sm)
{
  sm.set_halfedge(f, h);
}


//
// FaceListGraph
//
inline
boost::graph_traits<SurfaceMesh>::faces_size_type
num_faces(const SurfaceMesh& sm)
{
  return sm.faces_size();
}


inline
CGAL::Iterator_range<boost::graph_traits<SurfaceMesh>::face_iterator>
faces(const SurfaceMesh& sm)
{
  return CGAL::make_range(sm.faces_begin(), sm.faces_end());
}


inline
boost::graph_traits<SurfaceMesh>::vertex_descriptor
add_vertex(SurfaceMesh& sm) {

  return sm.add_vertex();
}


inline
boost::graph_traits<SurfaceMesh>::vertex_descriptor
add_vertex(const Point& p, SurfaceMesh& sm) {
  return sm.add_vertex(p);
}

// MutableGraph
inline
void
reserve(SurfaceMesh& sm,
        boost::graph_traits< SurfaceMesh >::vertices_size_type nv,
        boost::graph_traits< SurfaceMesh >::edges_size_type ne,
        boost::graph_traits< SurfaceMesh >::faces_size_type nf)
{
  sm.reserve(nv, ne, nf);
}


inline
void
remove_vertex(boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
              SurfaceMesh& sm) {

  sm.delete_vertex(v);
}


#if 0
inline
void
remove_edge(boost::graph_traits<SurfaceMesh>::vertex_descriptor u,
            boost::graph_traits<SurfaceMesh>::vertex_descriptor v,
            SurfaceMesh& sm)
{
  boost::graph_traits<SurfaceMesh>::edge_descriptor e = edge(u, v, sm);
  sm.delete_edge(e);
}
#endif


inline
void
remove_edge(boost::graph_traits<SurfaceMesh>::edge_descriptor e,
            SurfaceMesh& sm)
{
  sm.delete_edge(sm.edge(halfedge(e,sm)));
}


inline
void
remove_edge(boost::graph_traits<SurfaceMesh>::edge_iterator eiter,
            SurfaceMesh& sm)
{
  remove_edge(*eiter,sm);
}


inline
void
remove_face(boost::graph_traits<SurfaceMesh>::face_descriptor f,
            SurfaceMesh& sm)
{
  sm.delete_face(f);
}


inline
boost::graph_traits<SurfaceMesh>::face_descriptor
add_face(SurfaceMesh& sm)
{
  return sm.add_face();
}


template<typename InputIterator>
boost::graph_traits<SurfaceMesh>::face_descriptor
add_face(InputIterator begin, InputIterator end, SurfaceMesh& sm)
{
  std::vector<boost::graph_traits<SurfaceMesh>::vertex_descriptor>
    v(begin, end);
  return sm.add_face(v);
}

inline
void normalize_border(const pmp::SurfaceMesh&)
{}

} // namespace pmp

#include <CGAL/boost/graph/properties_SurfaceMesh.h>

#endif // DOXYGEN_RUNNING

#endif // CGAL_BOOST_GRAPH_TRAITS_SURFACEMESH_H
