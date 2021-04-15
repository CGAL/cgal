// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef VDA_TEST_VDA_H
#define VDA_TEST_VDA_H 1

#include <CGAL/use.h>
#include <CGAL/Voronoi_diagram_2/Accessor.h>
#include <cassert>
#include "helper_functions.h"
#include <CGAL/Testsuite/use.h>

//==========================================================================
//==========================================================================
//==========================================================================

template<class VDA>
void test_vd_face_concept(const typename VDA::Face_handle& f, int dim)
{
  typedef typename VDA::Face                 Face;

  // types
  CGAL_USE_TYPE(typename Face::Halfedge);
  CGAL_USE_TYPE(typename Face::Vertex);
  typedef typename Face::Halfedge_handle     Halfedge_handle;
  CGAL_USE_TYPE(typename Face::Vertex_handle);
  CGAL_USE_TYPE(typename Face::Face_handle);

  CGAL_USE_TYPE(typename Face::Delaunay_graph);
  typedef typename Face::Delaunay_vertex_handle  Delaunay_vertex_handle;

  typedef typename Face::Ccb_halfedge_circulator CCBHC;

  if ( dim > 0 ) {
    // access methods
    Halfedge_handle e = f->halfedge();
    assert( e->face() == f );

    CCBHC hc_start = f->ccb();
    CCBHC hc = hc_start;
    test_circulator(hc);

    hc = hc_start;
    do {
      bool b1 = f->is_halfedge_on_ccb(hc);
      assert( b1 );
      ++hc;
    } while ( hc != hc_start );
  }

  bool b = f->is_unbounded();
  kill_warning(b);

  Delaunay_vertex_handle v = f->dual();
  kill_warning(v);

  assert( f->is_valid() );
}

//==========================================================================
//==========================================================================
//==========================================================================

template<class VDA>
void test_vd_vertex_concept(const typename VDA::Vertex_handle& v)
{
  typedef typename VDA::Vertex              Vertex;

  // types
  CGAL_USE_TYPE(typename Vertex::Halfedge);
  CGAL_USE_TYPE(typename Vertex::Face);
  typedef typename Vertex::Halfedge_handle  Halfedge_handle;
  CGAL_USE_TYPE(typename Vertex::Vertex_handle);
  CGAL_USE_TYPE(typename Vertex::Face_handle);

  typedef typename Vertex::size_type        size_type;
  typedef typename Vertex::Point_2          Point_2;

  CGAL_USE_TYPE(typename Vertex::Delaunay_graph);
  typedef typename Vertex::Delaunay_face_handle   Delaunay_face_handle;
  typedef typename Vertex::Delaunay_vertex_handle Delaunay_vertex_handle;

  typedef typename Vertex::Halfedge_around_vertex_circulator HAVC;

  // access methods
  Halfedge_handle e = v->halfedge();
  assert( e->target() == v );

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
    assert( b1 && b2 && b3 );
    kill_warning( b1 && b2 && b3 );
    ++hc;
  } while ( hc != hc_start );

  size_type d = v->degree();
  assert( d == deg );
  kill_warning( d );

  Delaunay_vertex_handle vv;
  Vertex vertex = *v;
  for (unsigned int i = 0; i < 3; i++) { vv = vertex.site(i); }
  kill_warning(vv);
  kill_warning(vertex);

  Point_2 p = v->point();
  Delaunay_face_handle f = v->dual();
  kill_warning(p);
  kill_warning(f);

  assert( v->is_valid() );
}

//==========================================================================
//==========================================================================
//==========================================================================

template<class VDA>
void test_vd_halfedge_concept(const typename VDA::Halfedge_handle& e)
{
  typedef typename VDA::Halfedge              Halfedge;

  // types
  typedef typename Halfedge::Halfedge_handle  Halfedge_handle;
  typedef typename Halfedge::Vertex_handle    Vertex_handle;
  typedef typename Halfedge::Face_handle      Face_handle;

  CGAL_USE_TYPE(typename Halfedge::Delaunay_graph);
  typedef typename Halfedge::Delaunay_edge    Delaunay_edge;

  typedef typename Halfedge::Delaunay_vertex_handle    Delaunay_vertex_handle;
  typedef typename Halfedge::Ccb_halfedge_circulator   CCBHC;

  // access methods
  Halfedge_handle e_opp = e->opposite();
  e_opp = e->twin();
  assert( e_opp->opposite() == e );

  Halfedge_handle e_next = e->next();
  Halfedge_handle e_prev = e->previous();

  assert( e_next->previous() == e );
  assert( e_prev->next() == e );

  Face_handle f = e->face();
  assert( f->is_halfedge_on_ccb(e) );

  Vertex_handle v_src, v_trg;
  if ( e->has_source() ) {
    v_src = e->source();
    assert( v_src == e_opp->target() );
  }
  if ( e->has_target() ) {
    v_trg = e->target();
    assert( v_trg == e_opp->source() );
  }

  CCBHC hc = e->ccb();
  test_circulator(hc);

  Delaunay_vertex_handle v;
  v = e->up();
  v = e->down();
  if ( e->has_source() ) { v = e->left(); }
  if ( e->has_target() ) { v = e->right(); }
  Delaunay_edge de = e->dual();
  kill_warning(v);
  kill_warning(de);

  bool b1 = e->has_source();
  bool b2 = e->has_target();
  bool bu = e->is_unbounded();
  assert( (bu && (!b1 || !b2)) || (!bu && (b1 && b2)) );

  bool bs = e->is_segment();
  bool br = e->is_ray();
  bool bb = e->is_bisector();
  assert( (bs && (b1 && b2)) || (bb && (!b1 && !b2)) ||
                  (br && ((!b1 && b2) || (b1 && !b2))) );

  assert( e->is_valid() );
  kill_warning( b1 && b2 && bu && bs && br && bb );
}

//==========================================================================
//==========================================================================
//==========================================================================

template<class VDA>
void test_vdqr_concept(const VDA&, const CGAL::Tag_false&) {}

template<class VDA>
void test_vdqr_concept(const VDA& vda, const CGAL::Tag_true&)
{
  typedef typename VDA::Locate_result              Locate_result;
  typedef typename VDA::Point_2                    Point_2;
  typedef typename VDA::Vertex_handle              Vertex_handle;
  typedef typename VDA::Halfedge_handle            Halfedge_handle;
  typedef typename VDA::Face_handle                Face_handle;

  // nothing to do if the Voronoi diagram is empty
  if ( vda.dual().number_of_vertices() == 0 ) { return; }

  Point_2 p(0,0);
  Locate_result lr = vda.locate(p);

  if ( Vertex_handle* vv = boost::get<Vertex_handle>(&lr) ) {
    Vertex_handle v = *vv;
    kill_warning(v);
  } else if ( Halfedge_handle* ee = boost::get<Halfedge_handle>(&lr) ) {
    Halfedge_handle e = *ee;
    kill_warning(e);
  } else if ( Face_handle* ff = boost::get<Face_handle>(&lr) ) {
    Face_handle f = *ff;
    kill_warning(f);
  }
}

template<class VDA>
void test_vdqr_concept(const VDA& vda)
{
  test_vdqr_concept(vda, typename VDA::Adaptation_traits::Has_nearest_site_2());
}

//==========================================================================
//==========================================================================
//==========================================================================



template<class VDA>
void test_vda(const VDA& vda)
{
  CGAL_USE_TYPE(typename VDA::Accessor);

  // testing types
  //--------------
  typedef typename VDA::Delaunay_graph                DG;
  typedef typename VDA::Adaptation_traits             AT;
  typedef typename VDA::Adaptation_policy             AP;

  typedef typename VDA::size_type                     size_type;

  CGAL_USE_TYPE(typename VDA::Point_2);
  CGAL_USE_TYPE(typename VDA::Site_2);
  CGAL_USE_TYPE(typename VDA::Locate_result);

  typedef typename VDA::Halfedge                      Halfedge;
  typedef typename VDA::Face                          Face;
  typedef typename VDA::Vertex                        Vertex;

  typedef typename VDA::Halfedge_handle               Halfedge_handle;
  typedef typename VDA::Face_handle                   Face_handle;
  typedef typename VDA::Vertex_handle                 Vertex_handle;

  CGAL_USE_TYPE(typename VDA::Delaunay_geom_traits);
  CGAL_USE_TYPE(typename VDA::Delaunay_edge);
  CGAL_USE_TYPE(typename VDA::Delaunay_face_handle);
  CGAL_USE_TYPE(typename VDA::Delaunay_vertex_handle);

  CGAL_USE_TYPE(typename VDA::Edge_iterator);
  typedef typename VDA::Halfedge_iterator             Halfedge_iterator;
  typedef typename VDA::Face_iterator                 Face_iterator;
  typedef typename VDA::Vertex_iterator               Vertex_iterator;
  typedef typename VDA::Unbounded_faces_iterator      Unbounded_faces_iterator;
  typedef typename VDA::Bounded_faces_iterator        Bounded_faces_iterator;
  typedef typename VDA::Unbounded_halfedges_iterator  UH_iterator;
  typedef typename VDA::Bounded_halfedges_iterator    BH_iterator;

  typedef typename VDA::Halfedge_around_vertex_circulator  HAVC;
  typedef typename VDA::Ccb_halfedge_circulator            CCBHC;

  typedef typename VDA::Site_iterator                 Site_iterator;

  // testing access methods
  //-----------------------
  const DG& dg = vda.dual();
  const AT& at = vda.adaptation_traits();
  const AP& ap = vda.adaptation_policy();
  kill_warning(dg);
  kill_warning(at);

  size_type t;
  t = vda.number_of_vertices();
  t = vda.number_of_faces();
  t = vda.number_of_halfedges();
  t = vda.number_of_connected_components();
  kill_warning(t);

  {
    Face_handle f;
    f = vda.unbounded_face();
    if ( vda.unbounded_faces_begin() != vda.unbounded_faces_end() ) {
      assert( f != Face_handle() );
    } else {
      assert( f == Face_handle() );
    }

    f = vda.bounded_face();
    if ( vda.bounded_faces_begin() != vda.bounded_faces_end() ) {
      assert( f != Face_handle() );
    } else {
      assert( f == Face_handle() );
    }

    kill_warning(f);
  }
  {
    Halfedge_handle e;
    e = vda.unbounded_halfedge();
    if ( vda.unbounded_halfedges_begin() != vda.unbounded_halfedges_end() ) {
      assert( e != Halfedge_handle() );
    } else {
      assert( e == Halfedge_handle() );
    }

    e = vda.bounded_halfedge();
    if ( vda.bounded_halfedges_begin() != vda.bounded_halfedges_end() ) {
      assert( e != Halfedge_handle() );
    } else {
      assert( e == Halfedge_handle() );
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
    assert( *ccbhc_e == *e );
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
    assert( *havc_e == *e );
    test_circulator( havc_e );
  }

  test_vdqr_concept(vda);

  // testing validity
  assert( vda.is_valid() );

  // testing existence of non-default constructors
  {
    VDA vda2(at);
    VDA vda3(at, ap);
    VDA vda4(at, ap, vda.dual().geom_traits());

    DG dg2(dg);

    VDA vda5(dg2);
    VDA vda6(dg2, false);
    VDA vda7(dg2, false, at);
    VDA vda8(dg2, false, at, ap);
    std::size_t n_dg2 = dg2.number_of_vertices();
    VDA vda9(dg2, true);
    if ( n_dg2 != 0 ) {
      assert( dg2.number_of_vertices() == 0 );
    }

    std::vector<typename AT::Site_2> vs;
    for (Site_iterator sit = vda.sites_begin();
         sit != vda.sites_end(); ++sit) {
      vs.push_back(*sit);
    }

    VDA vda10(vs.begin(), vs.end());
    VDA vda11(vs.begin(), vs.end(), at);
    VDA vda12(vs.begin(), vs.end(), at, ap);
    VDA vda13(vs.begin(), vs.end(), at, ap, vda.dual().geom_traits());
  }

  // testing copy constructor
  VDA vda_copy(vda);
  size_type nv = vda.number_of_vertices();
  size_type ne = vda.number_of_halfedges();
  size_type nf = vda.number_of_faces();
  assert( vda_copy.number_of_vertices() == nv );
  assert( vda_copy.number_of_halfedges() == ne );
  assert( vda_copy.number_of_faces() == nf );
  assert( vda_copy.is_valid() );
  kill_warning( nv + ne + nf );

  // testing assignment operator
  vda_copy.clear();
  vda_copy = vda;
  VDA vda_copy_2 = vda;
  assert( vda_copy.number_of_vertices() == nv );
  assert( vda_copy.number_of_halfedges() == ne );
  assert( vda_copy.number_of_faces() == nf );
  assert( vda_copy.is_valid() );

  // testing clear
  vda_copy_2.clear();
  assert( vda_copy_2.number_of_vertices() == 0 );
  assert( vda_copy_2.number_of_halfedges() == 0 );
  assert( vda_copy_2.number_of_faces() == 0 );
  assert( vda_copy_2.is_valid() );

  // testing swap
  vda_copy_2.swap(vda_copy);
  assert( vda_copy_2.number_of_vertices() == nv );
  assert( vda_copy_2.number_of_halfedges() == ne );
  assert( vda_copy_2.number_of_faces() == nf );

  assert( vda_copy.number_of_vertices() == 0 );
  assert( vda_copy.number_of_halfedges() == 0 );
  assert( vda_copy.number_of_faces() == 0 );

  assert( vda_copy.is_valid() );
  assert( vda_copy_2.is_valid() );

#ifndef VDA_TEST_RT
  // testing file I/O
  std::ofstream ofs("tmp.vd.cgal");
  assert( ofs );
  ofs << vda;
  ofs.close();

  std::ifstream ifs("tmp.vd.cgal");
  assert( ifs );
  ifs >> vda_copy;
  ifs.close();

  assert( vda_copy.number_of_vertices() == nv );
  assert( vda_copy.number_of_halfedges() == ne );
  assert( vda_copy.number_of_faces() == nf );

  assert( vda_copy.is_valid() );
#endif // VDA_TEST_RT

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

  // testing face concept
  for (Face_iterator it = vda.faces_begin();
       it != vda.faces_end(); ++it) {
    test_vd_face_concept<VDA>(it, vda.dual().dimension());
  }

  // testing bounded/unbounded iterators
  // testing unbounded faces iterator
  Unbounded_faces_iterator ufit;
  for (ufit = vda.unbounded_faces_begin();
       ufit != vda.unbounded_faces_end(); ++ufit) {
    assert( ufit->is_unbounded() );
  }

  if ( vda.unbounded_faces_begin() != vda.unbounded_faces_end() ) {
    for (ufit = --vda.unbounded_faces_end();
         ufit != vda.unbounded_faces_begin(); --ufit) {
      assert( ufit->is_unbounded() );
    }
    assert( ufit->is_unbounded() );
  }

  Bounded_faces_iterator bfit;
  for (bfit = vda.bounded_faces_begin();
       bfit != vda.bounded_faces_end(); ++bfit) {
    assert( !bfit->is_unbounded() );
  }

  if ( vda.bounded_faces_begin() != vda.bounded_faces_end() ) {
    for (bfit = --vda.bounded_faces_end();
         bfit != vda.bounded_faces_begin(); --bfit) {
      assert( !bfit->is_unbounded() );
    }
    assert( !bfit->is_unbounded() );
  }

  // testing unbounded halfedges iterator
  UH_iterator ueit;
  for (ueit = vda.unbounded_halfedges_begin();
       ueit != vda.unbounded_halfedges_end(); ++ueit) {
    assert( !ueit->has_source() || !ueit->has_target() );
  }

  if ( vda.unbounded_halfedges_begin() != vda.unbounded_halfedges_end() ) {
    for (ueit = --vda.unbounded_halfedges_end();
         ueit != vda.unbounded_halfedges_begin(); --ueit) {
      assert( !ueit->has_source() || !ueit->has_target() );
    }
    assert( !ueit->has_source() || !ueit->has_target() );
  }

  // testing bounded halfedges iterator
  BH_iterator beit;
  for (beit = vda.bounded_halfedges_begin();
       beit != vda.bounded_halfedges_end(); ++beit) {
    assert( beit->has_source() && beit->has_target() );
  }

  if ( vda.bounded_halfedges_begin() != vda.bounded_halfedges_end() ) {
    for (beit = --vda.bounded_halfedges_end();
         beit != vda.bounded_halfedges_begin(); --beit) {
      assert( beit->has_source() && beit->has_target() );
    }
    assert( beit->has_source() && beit->has_target() );
  }

  // testing generator iterator
  Site_iterator sit;
  for (sit = vda.sites_begin(); sit != vda.sites_end(); ++sit) {
    typename AT::Site_2 s = *sit;
    CGAL_USE(s);
  }

  if ( vda.sites_begin() != vda.sites_end() ) {
    for (sit = --vda.sites_end(); sit != vda.sites_begin(); --sit) {
      typename AT::Site_2 s = *sit;
      CGAL_USE(s);
    }
  }
}

#endif // VDA_TEST_VDA_H
