// Copyright (c) 2005 University of Crete (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_DIAGRAM_ADAPTOR_2_H
#define CGAL_VORONOI_DIAGRAM_ADAPTOR_2_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/iterator.h>
#include <CGAL/circulator.h>

#include <iterator>
#include <vector>

#include <CGAL/Voronoi_diagram_adaptor_2/Halfedge.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Face.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Vertex.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Circulator_adaptors.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Iterator_adaptors.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Handle_adaptor.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Validity_testers.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Dummy_iterator.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Unbounded_faces.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Degeneracy_tester_binders.h>

CGAL_BEGIN_NAMESPACE


template<class DG, class Tr>
class Voronoi_diagram_adaptor_2
{
 private:
  typedef Voronoi_diagram_adaptor_2<DG,Tr>   Self;
  typedef Triangulation_cw_ccw_2             CW_CCW_2;

 public:
  //-------
  // TYPES
  //-------

  // TYPES FOR THE DUAL GRAPH

  // the (triangulated) dual graph
  typedef DG                                   Dual_graph;
  typedef typename Dual_graph::Geom_traits     Geom_traits;  
  typedef Tr                                   Voronoi_traits;

  typedef typename Dual_graph::size_type       size_type;
  typedef size_type                            Size;

  typedef typename Dual_graph::Vertex_handle   Dual_vertex_handle;
  typedef typename Dual_graph::Face_handle     Dual_face_handle;
  typedef typename Dual_graph::Edge            Dual_edge;


  typedef typename Dual_graph::Finite_vertices_iterator 
  Dual_vertices_iterator;

  typedef typename Dual_graph::Finite_faces_iterator
  Dual_faces_iterator;

  typedef typename Dual_graph::Finite_edges_iterator
  Dual_edges_iterator;
  typedef typename Dual_graph::All_edges_iterator
  All_dual_edges_iterator;

  typedef typename Dual_graph::Edge_circulator
  Dual_edge_circulator;


  typedef typename Voronoi_traits::Edge_degeneracy_tester
  Edge_degeneracy_tester;
  typedef typename Voronoi_traits::Face_degeneracy_tester
  Face_degeneracy_tester;

 protected:
  // DEGENERACY TESTER BINDERS
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Edge_degeneracy_tester_binder<Self>
  Edge_degeneracy_tester_binder;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Face_degeneracy_tester_binder<Self>
  Face_degeneracy_tester_binder;

  // ITERATORS FOR EDGES
  typedef Filter_iterator<Dual_edges_iterator,Edge_degeneracy_tester_binder>
  Non_degenerate_edges_iterator;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Edge_iterator_adaptor
  <Self,Non_degenerate_edges_iterator>
  Edge_iterator_base;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Edge_validity_tester
  <Self,Edge_iterator_base>
  Edge_validity_tester;

 public:
  typedef Filter_iterator<Edge_iterator_base,Edge_validity_tester>
  Edge_iterator;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Halfedge_iterator_adaptor<Self>
  Halfedge_iterator;

  // THE HALFEDGE
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Halfedge<Self>     Halfedge;

 protected:
  // ITERATORS FOR FACES
  typedef Filter_iterator<Dual_vertices_iterator,Face_degeneracy_tester_binder>
  Non_degenerate_faces_iterator;

 public:
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Face_iterator_adaptor
  <Self,Non_degenerate_faces_iterator>
  Face_iterator;

  // THE FACE
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Face<Self>         Face;

 protected:
  // ITERATORS FOR VERTICES
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Vertex_validity_tester<Self>
  Vertex_validity_tester;

  typedef Filter_iterator<Dual_faces_iterator,Vertex_validity_tester>
  Non_degenerate_vertices_iterator;

 public:
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Vertex_iterator_adaptor
  <Self,Non_degenerate_vertices_iterator>
  Vertex_iterator;

  // THE VERTEX
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Vertex<Self>       Vertex;

 public:
  // HANDLES
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Handle_adaptor<Halfedge>  Halfedge_handle;
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Handle_adaptor<Vertex>    Vertex_handle;
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Handle_adaptor<Face>      Face_handle;

  // THE HOLES ITERATOR
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Dummy_iterator<Halfedge_handle>
  Holes_iterator;

  // CIRCULATORS
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Halfedge_around_vertex_circulator_adaptor
  <Halfedge_handle>
  Halfedge_around_vertex_circulator;

  typedef CGAL_VORONOI_DIAGRAM_2_NS::Ccb_halfedge_circulator_adaptor
  <Halfedge_handle>
  Ccb_halfedge_circulator;

 protected:
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Bounded_face_tester
  <Self,Non_degenerate_faces_iterator>
  Bounded_face_tester;

 public:
  typedef
  Filter_iterator<Non_degenerate_faces_iterator,Bounded_face_tester>
  Unbounded_faces_iterator;

 public:
  // PREDICATES
  //-----------
  bool is_degenerate_edge(const Dual_face_handle& f, int i) const
  {
    return edge_tester()(dual_, f, i);
  }

  bool is_degenerate_edge(const Dual_edge& e) const
  {
    return is_degenerate_edge(e.first, e.second);
  }

  bool is_degenerate_edge(const Dual_edge_circulator& ec) const
  {
    return is_degenerate_edge(*ec);
  }

  bool is_degenerate_edge(const Dual_edges_iterator& eit) const
  {
    return is_degenerate_edge(*eit);
  }

  bool is_degenerate_edge(const All_dual_edges_iterator& eit) const
  {
    return is_degenerate_edge(*eit);
  }

  bool has_empty_Voronoi_cell_interior(const Dual_vertex_handle& v) const
  {
    return face_tester()(dual(), v);
  }


public:
  struct Face_circulator {}; // 1. circulates through the Voronoi cells
			     //    that are neighbors of the given
			     //    Voronoi cell;
                             // 2. also circulates through the Voronoi
                             //    cells adjacent to a Voronoi vertex

  struct Vertex_circulator {}; // 1. circulates through the Voronoi
			       //    vertices at the boundary of a
			       //    Voronoi cell
                               // 2. circulates also through the
                               //    Voronoi vertices that are
                               //    neighbors (through edges) to a
                               //    given vertex.

  struct Edge_circulator {}; // 1. circulates through the Voronoi
			     //    vertices at the boundary of a
			     //    Voronoi cell
                             // 2. also circulates around the edges of
                             //    a Voronoi vertex

public:
  //--------------
  // CONSTRUCTORS
  //--------------
  Voronoi_diagram_adaptor_2(const Voronoi_traits& tr = Voronoi_traits())
    : dual_(), tr_(tr), bf_tester_(this) {}

  Voronoi_diagram_adaptor_2(const Dual_graph& dg,
			    const Voronoi_traits& tr = Voronoi_traits())
    : dual_(dg), tr_(tr), bf_tester_(this) {}


public:
  //------------------
  // ACCESS FUNCTIONS
  //------------------

  // DUAL
  const Dual_graph& dual() const { return dual_; }

  // VORONOI TRAITS
  const Voronoi_traits& voronoi_traits() const { return tr_; }

  // SIZE RELATED FUNCTIONS
  size_type size_of_vertices() const {
    size_type num_v = 0;
    for (Vertex_iterator it = vertices_begin();	it != vertices_end();
	 ++it, ++num_v) {}
    return num_v;
  }

  size_type size_of_faces() const {
    size_type num_f = 0;
    for (Face_iterator it = faces_begin(); it != faces_end();
	 ++it, ++num_f) {}
    return num_f;
  }

  size_type size_of_halfedges() const {
    size_type num_h = 0;
    for (Halfedge_iterator it = halfedges_begin(); it != halfedges_end();
	 ++it, ++num_h) {}
    return num_h;
  }

  size_type number_of_vertices() const  { return size_of_vertices(); }
  size_type number_of_faces() const     { return size_of_faces(); }
  size_type number_of_halfedges() const { return size_of_halfedges(); }

  // DEGENERACY TESTERS
  const Edge_degeneracy_tester& edge_tester() const {
    return tr_.edge_degeneracy_tester_object();
  }

  const Face_degeneracy_tester& face_tester() const {
    return tr_.face_degeneracy_tester_object();
  }

  // UNBOUNDED FACE
  Face_handle unbounded_face() const {
    return Face_handle(*unbounded_faces_begin());
  }

  // FACE ITERATORS
 private:
  Non_degenerate_faces_iterator non_degenerate_faces_begin() const {
    return filter_iterator( dual_.finite_vertices_end(),
			    Face_degeneracy_tester_binder(this),
			    dual_.finite_vertices_begin() );
  }

  Non_degenerate_faces_iterator non_degenerate_faces_end() const {
    return filter_iterator( dual_.finite_vertices_end(),
			    Face_degeneracy_tester_binder(this) );
  }

 public:
  Face_iterator faces_begin() const {
    return Face_iterator(this, non_degenerate_faces_begin());
  }

  Face_iterator faces_end() const {
    return Face_iterator(this, non_degenerate_faces_end());
  }

  Unbounded_faces_iterator unbounded_faces_begin() const {
    return filter_iterator( non_degenerate_faces_end(),
			    bf_tester_,
			    non_degenerate_faces_begin() );
  }

  Unbounded_faces_iterator unbounded_faces_end() const {
    return filter_iterator( non_degenerate_faces_end(),
			    bf_tester_ );
  }

  // EDGE ITERATORS
 private:
  Non_degenerate_edges_iterator non_degenerate_edges_begin() const {
    return filter_iterator( dual_.finite_edges_end(),
			    Edge_degeneracy_tester_binder(this),
			    dual_.finite_edges_begin() );
  }

  Non_degenerate_edges_iterator non_degenerate_edges_end() const {
    return filter_iterator( dual_.finite_edges_end(),
			    Edge_degeneracy_tester_binder(this) );
  }


  Edge_iterator_base edges_base_begin() const {
    return Edge_iterator_base(this, non_degenerate_edges_begin());
  }

  Edge_iterator_base edges_base_end() const {
    return Edge_iterator_base(this, non_degenerate_edges_end());
  }

 public:
  Edge_iterator edges_begin() const {
    return filter_iterator( edges_base_end(),
			    Edge_validity_tester(this),
			    edges_base_begin() );
  }

  Edge_iterator edges_end() const {
    return filter_iterator( edges_base_end(),
			    Edge_validity_tester(this) );
  }
  
  Halfedge_iterator halfedges_begin() const {
    return Halfedge_iterator(this, edges_begin());
  }

  Halfedge_iterator halfedges_end() const {
    return Halfedge_iterator(this, edges_end());
  }

  // VERTEX ITERATORS
 private:
  Non_degenerate_vertices_iterator non_degenerate_vertices_begin() const {
    return filter_iterator( dual_.finite_faces_end(),
			    Vertex_validity_tester(this),
			    dual_.finite_faces_begin() );
  }

  Non_degenerate_vertices_iterator non_degenerate_vertices_end() const {
    return filter_iterator( dual_.finite_faces_end(),
			    Vertex_validity_tester(this) );
  }

 public:
  Vertex_iterator vertices_begin() const {
    return Vertex_iterator(this, non_degenerate_vertices_begin());
  }

  Vertex_iterator vertices_end() const {
    return Vertex_iterator(this, non_degenerate_vertices_end());
  }

  // CIRCULATORS
  Ccb_halfedge_circulator ccb_halfedges(const Face_handle& f) const {
    return Ccb_halfedge_circulator(f->halfedge());
  }

  Ccb_halfedge_circulator ccb_halfedges(const Face_handle& f,
					const Halfedge_handle& he) const {
    CGAL_precondition( he.face() == f );
    return Ccb_halfedge_circulator(he);
  }


  Halfedge_around_vertex_circulator
  incident_halfedges(const Vertex_handle& v) const {
    return incident_halfedges(v, v->halfedge());
  }

  Halfedge_around_vertex_circulator
  incident_halfedges(const Vertex_handle& v, const Halfedge_handle& he) const {
    CGAL_precondition( he->vertex() == v );
    return Halfedge_around_vertex_circulator(he);
  }

  bool is_valid() const {
    bool valid = dual_.is_valid();
    for (Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it) {
      valid = valid && it->is_valid();
    }

    for (Face_iterator it = faces_begin(); it != faces_end(); ++it) {
      valid = valid && it->is_valid();
    }

    for (Halfedge_iterator it = halfedges_begin(); it != halfedges_end();
	 ++it) {
      // I HAVE TO FIGURE OUT WHY THE FOLLOWING WORKS AND DOES NOT
      // POSE ANY PROBLEMS... THIS HAS TO DO WITH WHAT HALFEDGE IS
      // CONSIDERED INFINITE AND WHICH NOT...
      if (  !dual_.is_infinite( it->dual_edge().first ) &&
	    !dual_.is_infinite( it->opposite()->dual_edge().first )  ) {
	valid = valid && it->is_valid();
      }
    }
    return valid;
  }


private:
  Dual_graph  dual_;
  Voronoi_traits tr_;
  Bounded_face_tester bf_tester_;
};


CGAL_END_NAMESPACE

// TO-DO-LIST
// ----------
// * write code for vertices iterator; they need to return Vertex as type
// * decide whether infinite halfedges are indeed halfedges...

#endif // CGAL_VORONOI_DIAGRAM_ADAPTOR_2_H
