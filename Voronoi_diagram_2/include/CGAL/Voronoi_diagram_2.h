// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_H
#define CGAL_VORONOI_DIAGRAM_2_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/iterator.h>
#include <CGAL/circulator.h>
#include <CGAL/tags.h>
#include <CGAL/use.h>
#include <CGAL/assertions.h>

#include <iostream>
#include <iterator>

#include <CGAL/Voronoi_diagram_2/Halfedge.h>
#include <CGAL/Voronoi_diagram_2/Face.h>
#include <CGAL/Voronoi_diagram_2/Vertex.h>
#include <CGAL/Voronoi_diagram_2/Circulator_adaptors.h>
#include <CGAL/Voronoi_diagram_2/Iterator_adaptors.h>
#include <CGAL/Voronoi_diagram_2/Handle_adaptor.h>
#include <CGAL/Voronoi_diagram_2/Validity_testers.h>
#include <CGAL/Voronoi_diagram_2/Unbounded_faces.h>
#include <CGAL/Voronoi_diagram_2/Unbounded_edges.h>
#include <CGAL/Voronoi_diagram_2/Degeneracy_tester_binders.h>
#include <CGAL/Voronoi_diagram_2/Connected_components.h>
#include <CGAL/Voronoi_diagram_2/Accessor.h>

#include <CGAL/Identity_policy_2.h>

#include <boost/variant.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>

namespace CGAL {

//=========================================================================
//=========================================================================
//=========================================================================

template<class DG, class AT, class AP = Identity_policy_2<DG,AT> >
class Voronoi_diagram_2
{
 private:
  typedef Voronoi_diagram_2<DG,AT,AP>        Self;
  typedef Triangulation_cw_ccw_2             CW_CCW_2;

  friend class CGAL_VORONOI_DIAGRAM_2_INS::Accessor<Self>;
 public:
  //-------
  // TYPES
  //-------

  // TYPES FOR THE DUAL GRAPH

  // the (triangulated) dual graph
  typedef DG                                          Delaunay_graph;
  typedef AT                                          Adaptation_traits;
  typedef AP                                          Adaptation_policy;

  typedef typename Delaunay_graph::Geom_traits        Delaunay_geom_traits;

  typedef typename Delaunay_graph::size_type          size_type;

  typedef typename Delaunay_graph::Vertex_handle      Delaunay_vertex_handle;
  typedef typename Delaunay_graph::Face_handle        Delaunay_face_handle;
  typedef typename Delaunay_graph::Edge               Delaunay_edge;

 protected:
  typedef typename Delaunay_graph::Edge_circulator    Dual_edge_circulator;
  typedef typename Delaunay_graph::Vertex_circulator  Dual_vertex_circulator;
  typedef typename Delaunay_graph::Face_circulator    Dual_face_circulator;

  typedef typename Delaunay_graph::Finite_vertices_iterator 
  Dual_vertices_iterator;

  typedef typename Delaunay_graph::Finite_faces_iterator
  Dual_faces_iterator;

  typedef typename Delaunay_graph::Finite_edges_iterator
  Dual_edges_iterator;
  typedef typename Delaunay_graph::All_edges_iterator
  All_dual_edges_iterator;

 protected:
  // TYPES FOR THE DEGENERACY TESTERS
  typedef typename Adaptation_policy::Has_site_inserter  Has_site_inserter;
  typedef typename Adaptation_policy::Has_site_remover   Has_site_remover;

  typedef typename Adaptation_policy::Edge_rejector
  Edge_rejector;

  typedef typename Adaptation_policy::Face_rejector
  Face_rejector;

 protected:
  // DEGENERACY TESTER BINDERS
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Edge_rejector_binder<Self>
  Edge_rejector_binder;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Face_rejector_binder<Self>
  Face_rejector_binder;

  // ITERATORS FOR EDGES
  typedef Filter_iterator<Dual_edges_iterator,Edge_rejector_binder>
  Non_degenerate_edges_iterator;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Edge_iterator_adaptor
  <Self,Non_degenerate_edges_iterator>
  Edge_iterator_base;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Edge_validity_tester
  <Self,Edge_iterator_base>
  Edge_validity_tester;

  typedef Filter_iterator<Edge_iterator_base,Edge_validity_tester>
  Valid_edges_iterator;

 public:
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Edge_iterator_adaptor
  <Self,Valid_edges_iterator,Tag_false>
  Edge_iterator;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Halfedge_iterator_adaptor
  <Self,Edge_iterator>
  Halfedge_iterator;

  // THE HALFEDGE
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Halfedge<Self>     Halfedge;

 protected:
  // ITERATORS FOR FACES
  typedef Filter_iterator<Dual_vertices_iterator,Face_rejector_binder>
  Non_degenerate_faces_iterator;

 public:
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Face_iterator_adaptor
  <Self,Non_degenerate_faces_iterator>
  Face_iterator;

  // THE FACE
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Face<Self>         Face;

 protected:
  // ITERATORS FOR VERTICES
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Vertex_validity_tester<Self>
  Vertex_validity_tester;

  typedef Filter_iterator<Dual_faces_iterator,Vertex_validity_tester>
  Non_degenerate_vertices_iterator;

 public:
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Vertex_iterator_adaptor
  <Self,Non_degenerate_vertices_iterator>
  Vertex_iterator;

  // THE VERTEX
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Vertex<Self>       Vertex;

 public:
  // HANDLES
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Handle_adaptor<Halfedge>  Halfedge_handle;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Handle_adaptor<Vertex>    Vertex_handle;
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Handle_adaptor<Face>      Face_handle;

  // CIRCULATORS
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Halfedge_around_vertex_circulator_adaptor
  <Halfedge>
  Halfedge_around_vertex_circulator;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Ccb_halfedge_circulator_adaptor
  <Halfedge>
  Ccb_halfedge_circulator;

  // THE BOUNDED AND UNBOUNDED FACES ITERATOR
 protected:
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Bounded_face_tester
  <Self,Non_degenerate_faces_iterator>
  Bounded_face_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Unbounded_face_tester
  <Self,Non_degenerate_faces_iterator>
  Unbounded_face_tester;

 protected:
  typedef
  Filter_iterator<Non_degenerate_faces_iterator,Bounded_face_tester>
  Unbounded_faces_iterator_base;

  typedef
  Filter_iterator<Non_degenerate_faces_iterator,Unbounded_face_tester>
  Bounded_faces_iterator_base;

 public:
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Face_iterator_adaptor
  <Self,Unbounded_faces_iterator_base>
  Unbounded_faces_iterator;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Face_iterator_adaptor
  <Self,Bounded_faces_iterator_base>
  Bounded_faces_iterator;

  // THE BOUNDED AND UNBOUNDED HALFEDGES ITERATOR
 protected:
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Bounded_edge_tester
  <Self,Edge_iterator>
  Bounded_edge_tester;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Unbounded_edge_tester
  <Self,Edge_iterator>
  Unbounded_edge_tester;

 protected:
  typedef
  Filter_iterator<Edge_iterator,Bounded_edge_tester>
  Unbounded_edges_iterator_base;

  typedef
  Filter_iterator<Edge_iterator,Unbounded_edge_tester>
  Bounded_edges_iterator_base;

 public:
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Halfedge_iterator_adaptor
  <Self,Unbounded_edges_iterator_base>
  Unbounded_halfedges_iterator;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Halfedge_iterator_adaptor
  <Self,Bounded_edges_iterator_base>
  Bounded_halfedges_iterator;

  // GENERATOR ITERATOR
 protected:
  struct Project_site_2
  {
    typedef typename Adaptation_traits::Site_2  Site_2;
    typedef Face                                argument_type;
    typedef Site_2                              result_type;

    Site_2 operator()(const Face& f) const {
      // here we construct an adaptation traits; ideally we should get
      // the adaptation traits from the outer class
      return Adaptation_traits().access_site_2_object()(f.dual());
    }
  };

 public:

  typedef boost::transform_iterator<Project_site_2, Face_iterator> Site_iterator;

  // ACCESSOR
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Accessor<Self>  Accessor;

protected:
  // POINT LOCATION RELATED TYPES
  typedef typename Adaptation_traits::Has_nearest_site_2  Has_nearest_site_2;
public:
  typedef typename Adaptation_traits::Point_2             Point_2;

  typedef boost::variant<Face_handle,Halfedge_handle,Vertex_handle>
  Locate_result;

private:
  typedef CGAL_VORONOI_DIAGRAM_2_INS::Find_valid_vertex<Self>
  Find_valid_vertex;

public:
  //--------------
  // CONSTRUCTORS
  //--------------
  Voronoi_diagram_2(const Adaptation_traits& at = Adaptation_traits(),
		    const Adaptation_policy& ap = Adaptation_policy(),
		    const Delaunay_geom_traits& gt = Delaunay_geom_traits())
    : dual_(gt), ap_(ap), at_(at) {}

  Voronoi_diagram_2(const Delaunay_graph& dg, bool swap_dg = false,
		    const Adaptation_traits& at = Adaptation_traits(),
		    const Adaptation_policy& ap = Adaptation_policy())
    : dual_(), ap_(ap), at_(at) {
    if ( swap_dg ) {
      dual_.swap(const_cast<Delaunay_graph&>(dg));
    } else {
      dual_ = dg;
    }
  }

  template<class Iterator>
  Voronoi_diagram_2(Iterator first, Iterator beyond,
		    const Adaptation_traits& at = Adaptation_traits(),
		    const Adaptation_policy& ap = Adaptation_policy(),
		    const Delaunay_geom_traits& gt = Delaunay_geom_traits())
    : dual_(first, beyond, gt), ap_(ap), at_(at) {}

  Voronoi_diagram_2(const Voronoi_diagram_2& other)
    : dual_(other.dual_), ap_(other.ap_), at_(other.at_) {}

  Self& operator=(const Self& other) {
    dual_ = other.dual_;
    ap_ = other.ap_;
    at_ = other.at_;
    return *this;
  }

public:
  //------------------
  // ACCESS FUNCTIONS
  //------------------

  // VORONOI FEATURES FROM DELAUNAY FEATURES
  Halfedge_handle dual(const Delaunay_edge& e) const {
    return Halfedge_handle( Halfedge(this, e.first, e.second) );
  }

  Face_handle dual(Delaunay_vertex_handle v) const {
    return Face_handle( Face(this, v) );
  }

  Vertex_handle dual(Delaunay_face_handle f) const {
    return Vertex_handle( Vertex(this, f) );
  }

  // DUAL
  const Delaunay_graph& dual() const { return dual_; }

  // VORONOI TRAITS
  const Adaptation_traits& adaptation_traits() const { return at_; }

  // ADAPTATION POLICY
  const Adaptation_policy& adaptation_policy() const { return ap_; }

  // SIZE RELATED METHODS
  size_type number_of_vertices() const {
    size_type num_v = 0;
    for (Vertex_iterator it = vertices_begin();	it != vertices_end();
	 ++it, ++num_v) {}
    return num_v;
  }

  size_type number_of_faces() const {
    size_type num_f = 0;
    for (Face_iterator it = faces_begin(); it != faces_end();
	 ++it, ++num_f) {}
    return num_f;
  }

  size_type number_of_halfedges() const {
    size_type num_h = 0;
    for (Halfedge_iterator it = halfedges_begin(); it != halfedges_end();
	 ++it, ++num_h) {}
    return num_h;
  }

  size_type number_of_connected_components() const {
    return CGAL_VORONOI_DIAGRAM_2_INS::Connected_components<Self>()(*this);
  }

  // DEGENERACY TESTERS -- THESE ARE UNDOCUMENTED

  // MAYBE THE FOLLOWING TWO METHODS SHOULD BE PRIVATE AND ACCESSED
  // ONLY THROUGH THE ACCESSOR
  const Edge_rejector& edge_rejector() const {
    return ap_.edge_rejector_object();
  }

  const Face_rejector& face_rejector() const {
    return ap_.face_rejector_object();
  }

  // UNBOUNDED/BOUNDED FACE
  Face_handle unbounded_face() const {
    if ( unbounded_faces_begin() != unbounded_faces_end() ) {
      return unbounded_faces_begin();
    }
    return Face_handle();
  }

  Face_handle bounded_face() const {
    if ( bounded_faces_begin() != bounded_faces_end() ) {
      return bounded_faces_begin();
    }
    return Face_handle();
  }

  // UNBOUNDED/BOUNDED EDGE
  Halfedge_handle unbounded_halfedge() const {
    if ( unbounded_halfedges_begin() != unbounded_halfedges_end() ) {
      return unbounded_halfedges_begin();
    }
    return Halfedge_handle();
  }

  Halfedge_handle bounded_halfedge() const {
    if ( bounded_halfedges_begin() != bounded_halfedges_end() ) {
      return bounded_halfedges_begin();
    }
    return Halfedge_handle();
  }

  // FACE ITERATORS
 private:
  Non_degenerate_faces_iterator non_degenerate_faces_begin() const {
    return CGAL::filter_iterator( dual_.finite_vertices_end(),
				  Face_rejector_binder(this),
				  dual_.finite_vertices_begin() );
  }

  Non_degenerate_faces_iterator non_degenerate_faces_end() const {
    return CGAL::filter_iterator( dual_.finite_vertices_end(),
				  Face_rejector_binder(this) );
  }

 public:
  Face_iterator faces_begin() const {
    return Face_iterator(this, non_degenerate_faces_begin());
  }

  Face_iterator faces_end() const {
    return Face_iterator(this, non_degenerate_faces_end());
  }

 private:
  Unbounded_faces_iterator_base unbounded_faces_base_begin() const {
    return CGAL::filter_iterator( non_degenerate_faces_end(),
				  Bounded_face_tester(this),
				  non_degenerate_faces_begin() );
  }

  Unbounded_faces_iterator_base unbounded_faces_base_end() const {
    return CGAL::filter_iterator( non_degenerate_faces_end(),
				  Bounded_face_tester(this) );
  }

  Bounded_faces_iterator_base bounded_faces_base_begin() const {
    return CGAL::filter_iterator( non_degenerate_faces_end(),
				  Unbounded_face_tester(this),
				  non_degenerate_faces_begin() );
  }

  Bounded_faces_iterator_base bounded_faces_base_end() const {
    return CGAL::filter_iterator( non_degenerate_faces_end(),
				  Unbounded_face_tester(this) );
  }

 public:
  Unbounded_faces_iterator unbounded_faces_begin() const {
    return Unbounded_faces_iterator(this, unbounded_faces_base_begin());
  }

  Unbounded_faces_iterator unbounded_faces_end() const {
    return Unbounded_faces_iterator(this, unbounded_faces_base_end());
  }

  Bounded_faces_iterator bounded_faces_begin() const {
    return Bounded_faces_iterator(this, bounded_faces_base_begin());
  }

  Bounded_faces_iterator bounded_faces_end() const {
    return Bounded_faces_iterator(this, bounded_faces_base_end());
  }

  // EDGE ITERATORS
 private:
  Non_degenerate_edges_iterator non_degenerate_edges_begin() const {
    return CGAL::filter_iterator( dual_.finite_edges_end(),
				  Edge_rejector_binder(this),
				  dual_.finite_edges_begin() );
  }

  Non_degenerate_edges_iterator non_degenerate_edges_end() const {
    return CGAL::filter_iterator( dual_.finite_edges_end(),
				  Edge_rejector_binder(this) );
  }


  Edge_iterator_base edges_base_begin() const {
    return Edge_iterator_base(this, non_degenerate_edges_begin());
  }

  Edge_iterator_base edges_base_end() const {
    return Edge_iterator_base(this, non_degenerate_edges_end());
  }

  Valid_edges_iterator valid_edges_begin() const {
    return CGAL::filter_iterator( edges_base_end(),
				  Edge_validity_tester(this),
				  edges_base_begin() );
  }

  Valid_edges_iterator valid_edges_end() const {
    return CGAL::filter_iterator( edges_base_end(),
				  Edge_validity_tester(this) );
  }

 public:
  Edge_iterator edges_begin() const {
    return Edge_iterator(this, valid_edges_begin());
  }

  Edge_iterator edges_end() const {
    return Edge_iterator(this, valid_edges_end());
  }

  Halfedge_iterator halfedges_begin() const {
    return Halfedge_iterator(this, edges_begin());
  }

  Halfedge_iterator halfedges_end() const {
    return Halfedge_iterator(this, edges_end());
  }

 protected:
  Unbounded_edges_iterator_base unbounded_edges_base_begin() const {
    return CGAL::filter_iterator( edges_end(),
				  Bounded_edge_tester(this),
				  edges_begin() );
  }

  Unbounded_edges_iterator_base unbounded_edges_base_end() const {
    return CGAL::filter_iterator( edges_end(),
				  Bounded_edge_tester(this) );
  }

  Bounded_edges_iterator_base bounded_edges_base_begin() const {
    return CGAL::filter_iterator( edges_end(),
				  Unbounded_edge_tester(this),
				  edges_begin() );
  }

  Bounded_edges_iterator_base bounded_edges_base_end() const {
    return CGAL::filter_iterator( edges_end(),
				  Unbounded_edge_tester(this) );
  }

 public:
  Unbounded_halfedges_iterator unbounded_halfedges_begin() const {
    return Unbounded_halfedges_iterator(this, unbounded_edges_base_begin());
  }

  Unbounded_halfedges_iterator unbounded_halfedges_end() const {
    return Unbounded_halfedges_iterator(this, unbounded_edges_base_end());
  }

  Bounded_halfedges_iterator bounded_halfedges_begin() const {
    return Bounded_halfedges_iterator(this, bounded_edges_base_begin());
  }

  Bounded_halfedges_iterator bounded_halfedges_end() const {
    return Bounded_halfedges_iterator(this, bounded_edges_base_end());
  }

  // VERTEX ITERATORS
 private:
  Non_degenerate_vertices_iterator non_degenerate_vertices_begin() const {
    return CGAL::filter_iterator( dual_.finite_faces_end(),
				  Vertex_validity_tester(this),
				  dual_.finite_faces_begin() );
  }

  Non_degenerate_vertices_iterator non_degenerate_vertices_end() const {
    return CGAL::filter_iterator( dual_.finite_faces_end(),
				  Vertex_validity_tester(this) );
  }

 public:
  Vertex_iterator vertices_begin() const {
    return Vertex_iterator(this, non_degenerate_vertices_begin());
  }

  Vertex_iterator vertices_end() const {
    return Vertex_iterator(this, non_degenerate_vertices_end());
  }

  // SITE ITERATOR
  Site_iterator sites_begin() const {
    return Site_iterator(faces_begin());    
  }

  Site_iterator sites_end() const {
    return Site_iterator(faces_end());
  }

  // CIRCULATORS
  Ccb_halfedge_circulator ccb_halfedges(const Face_handle& f) const {
    return Ccb_halfedge_circulator(*f->halfedge());
  }

  Ccb_halfedge_circulator ccb_halfedges(const Face_handle& f,
					const Halfedge_handle& he) const {
    CGAL_precondition( he->face() == f );
    CGAL_USE(f);
    return Ccb_halfedge_circulator(*he);
  }


  Halfedge_around_vertex_circulator
  incident_halfedges(const Vertex_handle& v) const {
    return incident_halfedges(v, v->halfedge());
  }

  Halfedge_around_vertex_circulator
  incident_halfedges(const Vertex_handle& v, const Halfedge_handle& he) const {
    CGAL_USE(v);
    CGAL_precondition( he->target() == v );
    return Halfedge_around_vertex_circulator(*he);
  }

  // POINT LOCATION
  //---------------
 private:
  Locate_result locate(const Point_2& , const Tag_false&) const {
    CGAL_assertion_msg(false, "Point location is not supported");
    return Locate_result();
  }

  Locate_result locate(const Point_2& p, const Tag_true&) const
  {
    CGAL_precondition( dual_.number_of_vertices() > 0 );

    typedef typename Adaptation_traits::Nearest_site_2  Nearest_site_2;
    typedef typename Nearest_site_2::result_type        Query_result;

    Nearest_site_2 nearest_site = at_.nearest_site_2_object();
    Query_result ns_qr = nearest_site(dual_, p);

    if ( const Delaunay_vertex_handle* dv =
	 boost::get<Delaunay_vertex_handle>(&ns_qr) ) {
      return Face_handle( Face(this, *dv) );
    } else if ( const Delaunay_face_handle *df =
		boost::get<Delaunay_face_handle>(&ns_qr) ) {
      Find_valid_vertex vertex_finder;
      Delaunay_face_handle dfvalid = vertex_finder(this, *df);
      return Vertex_handle( Vertex(this, dfvalid) );
    } else if ( const Delaunay_edge* de =
		boost::get<Delaunay_edge>(&ns_qr) ) {
      CGAL_assertion(  !edge_rejector()(dual_, *de)  );
      if ( dual_.dimension() == 1 ) {
	Delaunay_vertex_handle v1 =
	  de->first->vertex(CW_CCW_2::ccw(de->second));
	Delaunay_vertex_handle v2 =
	  de->first->vertex(CW_CCW_2::cw(de->second) );
	return Halfedge_handle( Halfedge(this, v1, v2) );
      }
      return Halfedge_handle( Halfedge(this, de->first, de->second) );
    }

    // I should never have reached this line;
    CGAL_error();
    return Locate_result();
  }

 public:
  Locate_result locate(const Point_2& p) const {
    return locate(p, Has_nearest_site_2());
  }


  // VALIDITY TESTING
  //-----------------
  bool is_valid() const {
    bool valid = dual_.is_valid();
    valid = valid && ap_.is_valid();
    valid = valid && ap_.is_valid(dual_);

    for (Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it) {
      valid = valid && it->is_valid();
    }

    for (Face_iterator it = faces_begin(); it != faces_end(); ++it) {
      valid = valid && it->is_valid();
    }

    for (Halfedge_iterator it = halfedges_begin(); it != halfedges_end();
	 ++it) {
      valid = valid && it->is_valid();
    }

    return valid;
  }
  
  // I/O
  //----
 public:
  void file_output(std::ostream& os) const { os << dual_; }
  void file_input(std::istream& is) { is >> dual_; }

  // MISCALLANEOUS
  //--------------
 public:
  void clear() {
    dual_.clear();
    ap_.clear();
  }

  void swap(Self& other) {
    dual_.swap(other.dual_);
    ap_.swap(other.ap_);
  }

  // ACCESSOR
  //---------
  Accessor accessor() const { return Accessor(this); }
  Accessor accessor() { return Accessor(this); }

 private:
  Delaunay_graph  dual_;
  Adaptation_policy ap_;
  Adaptation_traits at_;

 protected:
  Delaunay_edge opposite(const Delaunay_edge& e) const {
    int i_mirror = dual_.tds().mirror_index(e.first, e.second);
    return Dual_edge( e.first->neighbor(e.second), i_mirror );
  }

 public:
  typedef typename Adaptation_traits::Site_2     Site_2;

 protected:
  // insert is supported...
  inline Face_handle insert(const Site_2& t, const Tag_true&) {
    Delaunay_vertex_handle v = ap_.site_inserter_object()(dual_, t);
    if ( v == Delaunay_vertex_handle() ) { return Face_handle(); }
    return Face_handle( Face(this, v) );
  }

  // insert is not really supported...
  inline Face_handle insert(const Site_2& t, const Tag_false&) {
    INSERT_IS_NOT_SUPPORTED(t);
    return Face_handle();
  }

 public:
  inline Face_handle insert(const Site_2& t) {
#if 0
    // THE FOLLOWING LINE MAY BE ADDED FOR DEBUGGING PURPOSES
    for (Halfedge_iterator it=halfedges_begin();it!=halfedges_end();++it) ;
#endif
    return insert(t, Has_site_inserter());
  }

  template<class Iterator>
  inline size_type insert(Iterator first, Iterator beyond) {
    size_type counter = 0;
    for (Iterator it = first; it != beyond; ++it, ++counter) {
      insert(*it);
    }
    return counter;
  }

#if 0
 // REMOVAL IS NOT READY YET
 protected:
  void update_cached_testers(const Face_handle& f, const Tag_true&)
  {
    // HERE WE ALSO NEED TO ACCOUNT FOR THE DELETED VERTEX. THE
    // DELETED VERTEX CAN AFFECT THE CACHED FACE DEGENERACY TESTER
    Dual_vertex_handle v = f->dual();
    Dual_edge_circulator ec_start = dual_.incident_edges(v);
    Dual_edge_circulator ec = ec_start;
    do {
      cached_e_tester_.erase(*ec);
    } while ( ++ec != ec_start );
  }

  void update_cached_testers(const Face_handle& f, const Tag_false&) {}

  void remove(const Face_handle& f, const Tag_true&) {
    if ( f == Face_handle() ) { return; }
    update_cached_testers(f, !Has_insert() && Has_get_conflicts());
    dual_.remove(f->dual());
  }

  void remove(const Face_handle& f, const Tag_false&) {
    REMOVE_IS_NOT_SUPPORTED(f);
  }

 public:
  inline void remove(const Face_handle& f) {
    return remove(f, Has_remove());
  }
#endif
};


// I/O OPERATORS
//--------------
template<class DG, class AT, class AP>
std::ostream& operator<<(std::ostream& os,
			 const Voronoi_diagram_2<DG,AT,AP>& vd)
{
  vd.file_output(os);
  return os;
}


template<class DG, class AT, class AP>
std::istream& operator>>(std::istream& is,
			 Voronoi_diagram_2<DG,AT,AP>& vd)
{
  vd.file_input(is);
  return is;
}


} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_H
