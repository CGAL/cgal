// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
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

#ifndef VDA_TEST_VDA_H
#define VDA_TEST_VDA_H 1

#include <CGAL/basic.h>
#include <CGAL/Voronoi_diagram_adaptor_2/Accessor.h>
#include <cassert>
#include "helper_functions.h"

//==========================================================================
//==========================================================================
//==========================================================================

template<class VDA>
void test_vd_face_concept(const typename VDA::Face_handle& f)
{
  typedef typename VDA::Face                Face;

  // types  
  typedef typename Face::Halfedge            Halfedge;
  typedef typename Face::Vertex              Vertex;
  typedef typename Face::Halfedge_handle     Halfedge_handle;
  typedef typename Face::Vertex_handle       Vertex_handle;
  typedef typename Face::Face_handle         Face_handle;

  typedef typename Face::size_type           size_type;
  typedef typename Face::Point_2             Point_2;

  typedef typename Face::Dual_graph          Dual_graph;
  typedef typename Face::Dual_vertex_handle  Dual_vertex_handle;

  typedef typename Face::Holes_iterator      Holes_iterator;

  typedef typename Face::Ccb_halfedge_circulator CCBHC;


  // access methods
  Halfedge_handle e = f->halfedge();
  CGAL_assertion( e->face() == f );

  CCBHC hc_start = f->outer_ccb();
  CCBHC hc = hc_start;
  test_circulator(hc);

  hc = hc_start;
  do {
    bool b1 = f->is_halfedge_on_outer_ccb(hc);
    bool b2 = f->is_halfedge_on_inner_ccb(hc);
    CGAL_assertion( b1 && !b2 );
    hc++;
  } while ( hc != hc_start );

  bool b = f->is_unbounded();
  kill_warning(b);

  test_iterator(f->holes_begin(), f->holes_end());
  test_is_convertible_to<Halfedge_handle>(f->holes_begin());
  test_is_convertible_to<Halfedge>(*f->holes_begin());

  Dual_vertex_handle v = f->dual_vertex();
  kill_warning(v);

  CGAL_assertion( f->is_valid() );
}

//==========================================================================
//==========================================================================
//==========================================================================

template<class VDA>
void test_vd_vertex_concept(const typename VDA::Vertex_handle& v)
{
  typedef typename VDA::Vertex                Vertex;

  // types  
  typedef typename Vertex::Halfedge         Halfedge;
  typedef typename Vertex::Face             Face;
  typedef typename Vertex::Halfedge_handle  Halfedge_handle;
  typedef typename Vertex::Vertex_handle    Vertex_handle;
  typedef typename Vertex::Face_handle      Face_handle;

  typedef typename Vertex::size_type        size_type;
  typedef typename Vertex::Point_2          Point_2;

  typedef typename Vertex::Dual_graph       Dual_graph;
  typedef typename Vertex::Dual_face_handle Dual_face_handle;

  typedef typename Vertex::Halfedge_around_vertex_circulator HAVC;

  // access methods
  Halfedge_handle e = v->halfedge();
  CGAL_assertion( e->target() == v );

  HAVC hc_start = v->incident_halfedges();
  HAVC hc = hc_start;
  test_circulator(hc);
  hc = hc_start;
  size_type deg = 0;
  do {
    deg++;
    bool b1 = v->is_incident_edge(hc);
    bool b2 = v->is_incident_face(hc->face());
    bool b3 = v->is_incident_face(hc->opposite()->face());
    CGAL_assertion( b1 && b2 && b3 );
    kill_warning( b1 && b2 && b3 );
    ++hc;
  } while ( hc != hc_start );

  size_type d = v->degree();
  CGAL_assertion( d == deg );
  kill_warning( d );

  Point_2 p = v->point();
  Dual_face_handle f = v->dual_face();
  kill_warning(p);
  kill_warning(f);

  CGAL_assertion( v->is_valid() );
}

//==========================================================================
//==========================================================================
//==========================================================================

template<class VDA>
void test_vd_halfedge_concept(const typename VDA::Halfedge_handle& e)
{
  typedef typename VDA::Halfedge              Halfedge;

  // types  
  typedef typename Halfedge::Vertex           Vertex;
  typedef typename Halfedge::Face             Face;
  typedef typename Halfedge::Halfedge_handle  Halfedge_handle;
  typedef typename Halfedge::Vertex_handle    Vertex_handle;
  typedef typename Halfedge::Face_handle      Face_handle;

  typedef typename Halfedge::Curve            Curve;

  typedef typename Halfedge::Dual_graph       Dual_graph;
  typedef typename Halfedge::Dual_edge        Dual_edge;

  typedef typename Halfedge::Ccb_halfedge_circulator CCBHC;

  // access methods
  Halfedge_handle e_opp = e->opposite();
  e_opp = e->twin();
  CGAL_assertion( e_opp->opposite() == e );

  Halfedge_handle e_next = e->next();
  Halfedge_handle e_prev = e->previous();

  CGAL_assertion( e_next->previous() == e );
  CGAL_assertion( e_prev->next() == e );

  Face_handle f = e->face();
  CGAL_assertion( f->is_halfedge_on_outer_ccb(e) );

  Vertex_handle v_src, v_trg;
  if ( e->has_source() ) {
    v_src = e->source();
    CGAL_assertion( v_src == e_opp->target() );
  }
  if ( e->has_target() ) {
    v_trg = e->target();
    CGAL_assertion( v_trg == e_opp->source() );
  }

  CCBHC hc = e->ccb();
  test_circulator(hc);

  Curve c = e->curve();
  Dual_edge de = e->dual_edge();
  kill_warning(c);
  kill_warning(de);

  bool b1 = e->has_source();
  bool b2 = e->has_target();
  bool bu = e->is_unbounded();
  CGAL_assertion( (bu && (!b1 || !b2)) || (!bu && (b1 && b2)) );

  CGAL_assertion( e->is_valid() );
  kill_warning( b1 && b2 && bu );
}

//==========================================================================
//==========================================================================
//==========================================================================

template<class VDA>
void test_vda(const VDA& vda)
{
  typedef CGAL_VORONOI_DIAGRAM_2_NS::Accessor<VDA>    Accessor;

  // testing types
  //--------------
  typedef typename VDA::Dual_graph                    DG;
  typedef typename VDA::Voronoi_traits                VT;

  typedef typename VDA::size_type                     size_type;

  typedef typename VDA::Point_2                       Point_2;
  typedef typename VDA::Locate_result                 Locate_result;

  typedef typename VDA::Halfedge                      Halfedge;
  typedef typename VDA::Face                          Face;
  typedef typename VDA::Vertex                        Vertex;

  typedef typename VDA::Halfedge_handle               Halfedge_handle;
  typedef typename VDA::Face_handle                   Face_handle;
  typedef typename VDA::Vertex_handle                 Vertex_handle;

  typedef typename VDA::Edge_iterator                 Edge_iterator;
  typedef typename VDA::Halfedge_iterator             Halfedge_iterator;
  typedef typename VDA::Face_iterator                 Face_iterator;
  typedef typename VDA::Vertex_iterator               Vertex_iterator;
  typedef typename VDA::Holes_iterator                Holes_iterator;
  typedef typename VDA::Unbounded_faces_iterator      Unbounded_faces_iterator;
  typedef typename VDA::Bounded_faces_iterator        Bounded_faces_iterator;
  typedef typename VDA::Unbounded_halfedges_iterator  UH_iterator;
  typedef typename VDA::Bounded_halfedges_iterator    BH_iterator;

  typedef typename VDA::Halfedge_around_vertex_circulator  HAVC;
  typedef typename VDA::Ccb_halfedge_circulator            CCBHC;

  
  // testing access methods
  //-----------------------
  const DG& dg = vda.dual();
  const VT& vt = vda.voronoi_traits();
  kill_warning(dg);
  kill_warning(vt);

  size_type t;
  t = vda.number_of_vertices();
  t = vda.number_of_faces();
  t = vda.number_of_halfedges();
  kill_warning(t);

  {
    Face_handle f;
    f = vda.unbounded_face();
    if ( vda.unbounded_faces_begin() != vda.unbounded_faces_end() ) {
      CGAL_assertion( f != Face_handle() );
    } else {
      CGAL_assertion( f == Face_handle() );
    }

    f = vda.bounded_face();
    if ( vda.bounded_faces_begin() != vda.bounded_faces_end() ) {
      CGAL_assertion( f != Face_handle() );
    } else {
      CGAL_assertion( f == Face_handle() );
    }

    kill_warning(f);
  }
  {
    Halfedge_handle e;
    e = vda.unbounded_halfedge();
    if ( vda.unbounded_halfedges_begin() != vda.unbounded_halfedges_end() ) {
      CGAL_assertion( e != Halfedge_handle() );
    } else {
      CGAL_assertion( e == Halfedge_handle() );
    }

    e = vda.bounded_halfedge();
    if ( vda.bounded_halfedges_begin() != vda.bounded_halfedges_end() ) {
      CGAL_assertion( e != Halfedge_handle() );
    } else {
      CGAL_assertion( e == Halfedge_handle() );
    }

    kill_warning(e);
  }

  // testing traversal methods
  //--------------------------

  // various face iterators
  test_iterator(vda.faces_begin(), vda.faces_end());
  test_iterator(vda.unbounded_faces_begin(), vda.unbounded_faces_end());
  test_iterator(vda.bounded_faces_begin(), vda.bounded_faces_end());
  if ( vda.faces_begin() != vda.faces_end() ) {
    test_is_convertible_to<Face_handle>(vda.faces_begin());
    test_is_convertible_to<Face>(*vda.faces_begin());
  }
  if ( vda.unbounded_faces_begin() != vda.unbounded_faces_end() ) {
    test_is_convertible_to<Face_handle>(vda.unbounded_faces_begin());
    test_is_convertible_to<Face>(*vda.unbounded_faces_begin());
  }
  if ( vda.bounded_faces_begin() != vda.bounded_faces_end() ) {
    test_is_convertible_to<Face_handle>(vda.bounded_faces_begin());
    test_is_convertible_to<Face>(*vda.bounded_faces_begin());
  }

  // various edge iterators
  test_iterator(vda.edges_begin(), vda.edges_end());
  test_iterator(vda.halfedges_begin(), vda.halfedges_end());
  test_iterator(vda.unbounded_halfedges_begin(),
		vda.unbounded_halfedges_end());
  test_iterator(vda.bounded_halfedges_begin(), vda.bounded_halfedges_end());
  if ( vda.edges_begin() != vda.edges_end() ) {
    test_is_convertible_to<Halfedge_handle>(vda.edges_begin());
    test_is_convertible_to<Halfedge>(*vda.edges_begin());
  }
  if ( vda.halfedges_begin() != vda.halfedges_end() ) {
    test_is_convertible_to<Halfedge_handle>(vda.halfedges_begin());
    test_is_convertible_to<Halfedge>(*vda.halfedges_begin());
  }
  if ( vda.unbounded_halfedges_begin() != vda.unbounded_halfedges_end() ) {
    test_is_convertible_to<Halfedge_handle>(vda.unbounded_halfedges_begin());
    test_is_convertible_to<Halfedge>(*vda.unbounded_halfedges_begin());
  }
  if ( vda.bounded_halfedges_begin() != vda.bounded_halfedges_end() ) {
    test_is_convertible_to<Halfedge_handle>(vda.bounded_halfedges_begin());
    test_is_convertible_to<Halfedge>(*vda.bounded_halfedges_begin());
  }

  // vertex iterators
  test_iterator(vda.vertices_begin(), vda.vertices_end());
  if ( vda.vertices_begin() != vda.vertices_end() ) {
    test_is_convertible_to<Vertex_handle>(vda.vertices_begin());
    test_is_convertible_to<Vertex>(*vda.vertices_begin());
  }

  // holes iterator
  //  test_iterator(vda.holes_begin(), vda.holes_end());
  //  test_is_convertible_to<Halfedge_handle>(vda.holes_begin());
  //  test_is_convertible_to<Halfedge>(*vda.holes_begin());

  // Ccb_halfedge_circulator
  if ( vda.faces_begin() != vda.faces_end()
       && vda.dual().dimension() > 0 ) {
    Face_handle f = vda.faces_begin();
    Halfedge_handle e = f->halfedge();
    CCBHC ccbhc = vda.ccb_halfedges(f);
    test_circulator( ccbhc );
    test_is_convertible_to<Halfedge_handle>(ccbhc);
    test_is_convertible_to<Halfedge>(*ccbhc);

    CCBHC ccbhc_e = vda.ccb_halfedges(f, e);
    // test that initialization was correct
    CGAL_assertion( *ccbhc_e == *e );
    test_circulator( ccbhc_e );
  }

  // Halfedge_around_vertex_circulator
  if ( vda.vertices_begin() != vda.vertices_end() ) {
    Vertex_handle v = vda.vertices_begin();
    Halfedge_handle e = v->halfedge();
    HAVC havc = vda.incident_halfedges(v);
    test_circulator( havc );
    test_is_convertible_to<Halfedge_handle>(havc);
    test_is_convertible_to<Halfedge>(*havc);

    HAVC havc_e = vda.incident_halfedges(v, e);
    // test that initialization was correct
    CGAL_assertion( *havc_e == *e );
    test_circulator( havc_e );
  }

  test_vdlr_concept(vda);

  // testing validity
  assert( vda.is_valid() );

  // sanity tests
  //-------------
  // testing halfedge concept
  for (Halfedge_iterator it = vda.halfedges_begin();
       it != vda.halfedges_end(); ++it) {
    test_vd_halfedge_concept<VDA>(it);
  }

  // testing vertex concept
  for (Vertex_iterator it = vda.vertices_begin();
       it != vda.vertices_end(); ++it) {
    test_vd_vertex_concept<VDA>(it);
  }

  // testing vertex concept
  for (Face_iterator it = vda.faces_begin();
       it != vda.faces_end(); ++it) {
    //    test_vd_face_concept<VDA>(it);
  }

  // testing bounded/unbounded iterators
  // testing unbounded faces iterator
  Unbounded_faces_iterator ufit;
  for (ufit = vda.unbounded_faces_begin();
       ufit != vda.unbounded_faces_end(); ++ufit) {
    CGAL_assertion( ufit->is_unbounded() );
  }

  if ( vda.unbounded_faces_begin() != vda.unbounded_faces_end() ) {
    for (ufit = --vda.unbounded_faces_end();
	 ufit != vda.unbounded_faces_begin(); --ufit) {
      CGAL_assertion( ufit->is_unbounded() );
    }
    CGAL_assertion( ufit->is_unbounded() );
  }

  Bounded_faces_iterator bfit;
  for (bfit = vda.bounded_faces_begin();
       bfit != vda.bounded_faces_end(); ++bfit) {
    CGAL_assertion( !bfit->is_unbounded() );
  }

  if ( vda.bounded_faces_begin() != vda.bounded_faces_end() ) {
    for (bfit = --vda.bounded_faces_end();
	 bfit != vda.bounded_faces_begin(); --bfit) {
      CGAL_assertion( !bfit->is_unbounded() );
    }
    CGAL_assertion( !bfit->is_unbounded() );
  }

  // testing unbounded halfedges iterator
  UH_iterator ueit;
  for (ueit = vda.unbounded_halfedges_begin();
       ueit != vda.unbounded_halfedges_end(); ++ueit) {
    CGAL_assertion( !ueit->has_source() || !ueit->has_target() );
  }

  if ( vda.unbounded_halfedges_begin() != vda.unbounded_halfedges_end() ) {
    for (ueit = --vda.unbounded_halfedges_end();
	 ueit != vda.unbounded_halfedges_begin(); --ueit) {
      CGAL_assertion( !ueit->has_source() || !ueit->has_target() );
    }
    CGAL_assertion( !ueit->has_source() || !ueit->has_target() );
  }

  // testing bounded halfedges iterator
  BH_iterator beit;
  for (beit = vda.bounded_halfedges_begin();
       beit != vda.bounded_halfedges_end(); ++beit) {
    CGAL_assertion( beit->has_source() && beit->has_target() );
  }

  if ( vda.bounded_halfedges_begin() != vda.bounded_halfedges_end() ) {
    for (beit = --vda.bounded_halfedges_end();
	 beit != vda.bounded_halfedges_begin(); --beit) {
      CGAL_assertion( beit->has_source() && beit->has_target() );
    }
    CGAL_assertion( beit->has_source() && beit->has_target() );
  }
}


template<class VDA>
void test_vdlr_concept(const VDA& vda)
{
  test_vdlr_concept(vda, typename VDA::Voronoi_traits::Has_point_locator());
}

template<class VDA>
void test_vdlr_concept(const VDA&, const CGAL::Tag_false&) {}

template<class VDA>
void test_vdlr_concept(const VDA& vda, const CGAL::Tag_true&)
{
  typedef typename VDA::Locate_result              Locate_result;
  typedef typename VDA::Point_2                    Point_2;
  typedef typename Locate_result::Vertex_handle    Vertex_handle;
  typedef typename Locate_result::Halfedge_handle  Halfedge_handle;
  typedef typename Locate_result::Face_handle      Face_handle;

  // nothing to do if the Voronoi diagram is empty
  if ( vda.dual().number_of_vertices() == 0 ) { return; }

  Point_2 p(0,0);
  Locate_result lr = vda.locate(p);

  if ( lr.is_vertex() ) {
    Vertex_handle v = lr;
    kill_warning(v);
  } else if ( lr.is_edge() ) {
    Halfedge_handle e = lr;
    kill_warning(e);    
  } else if ( lr.is_face() ) {
    Face_handle f = lr;
    kill_warning(f);
  }
}


#endif // VDA_TEST_VDA_H
