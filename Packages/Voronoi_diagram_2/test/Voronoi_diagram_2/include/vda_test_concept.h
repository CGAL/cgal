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

#ifndef CGAL_VDA_TEST_CONCEPT_H
#define CGAL_VDA_TEST_CONCEPT_H 1

#include <CGAL/Triangulation_utils_2.h>
#include "helper_functions.h"

//============================================================================
//============================================================================

template<class E, class Other>
void test_edge_as_pair(const E& e, const Other& o)
{
  dont_compile(e, o);
}

template<class E, class F_handle>
void test_edge_as_pair(const E&, const std::pair<F_handle,int>&)
{}

//============================================================================
//============================================================================

template<class DG>
void test_dual_graph_concept(const DG& dg)
{
  // testing types
  typedef typename DG::size_type                      size_type;
  typedef typename DG::Geom_traits                    Geom_traits;
  typedef typename DG::Triangulation_data_structure   Tds;

  typedef typename DG::Vertex                         Vertex;
  typedef typename DG::Face                           Face;
  typedef typename DG::Edge                           Edge;

  typedef typename DG::Vertex_handle                  Vertex_handle;
  typedef typename DG::Face_handle                    Face_handle;

  // testing if edges are defined as pairs of face handles and integers
  test_edge_as_pair(Edge(), std::pair<Face_handle,int>());

  typedef typename DG::All_edges_iterator             All_edges_iterator;
  typedef typename DG::Finite_edges_iterator          Finite_edges_iterator;

  typedef typename DG::All_vertices_iterator          All_vertices_iterator;
  typedef typename DG::Finite_vertices_iterator       Finite_vertices_iterator;

  typedef typename DG::All_faces_iterator             All_faces_iterator;
  typedef typename DG::Finite_faces_iterator          Finite_faces_iterator;

  typedef typename DG::Face_circulator                Face_circulator;
  typedef typename DG::Edge_circulator                Edge_circulator;
  typedef typename DG::Vertex_circulator              Vertex_circulator;

  // testing convertibility of iterators and circulators to handles
  if ( dg.number_of_vertices() > 0 ) {
    test_is_convertible_to<Vertex_handle>(dg.all_vertices_begin());
    test_is_convertible_to<Vertex_handle>(dg.finite_vertices_begin());
  }

  if ( dg.number_of_faces() > 0 ) {
    test_is_convertible_to<Face_handle>(dg.all_faces_begin());
    test_is_convertible_to<Face_handle>(dg.finite_faces_begin());
  }

  if ( dg.dimension() == 2 ) {
    Vertex_circulator vc = dg.incident_vertices(dg.finite_vertex());
    test_is_convertible_to<Vertex_handle>(vc);

    Face_circulator fc = dg.incident_faces(dg.finite_vertex());
    test_is_convertible_to<Face_handle>(fc);
  }

  // testing existence of access methods
  const Tds& tds = dg.tds();
  kill_warning(tds);

  const Geom_traits& gt = dg.geom_traits();
  kill_warning(gt);

  Vertex_handle v_inf = dg.infinite_vertex();
  kill_warning(v_inf);

  if ( dg.number_of_vertices() > 0 ) {
    Vertex_handle v_fin = dg.finite_vertex();
    kill_warning(v_fin);
  }

  if ( dg.dimension() > 1 ) {
    Face_handle f_inf = dg.infinite_face();
    kill_warning(f_inf);
  }

  int d = dg.dimension();
  kill_warning(d);

  size_type s = dg.number_of_vertices();
  s = dg.number_of_faces();
  kill_warning(s);

  // testing traversal with iterators
  test_iterator(dg.finite_vertices_begin(), dg.finite_vertices_end());
  test_iterator(dg.finite_edges_begin(), dg.finite_edges_end());
  test_iterator(dg.finite_faces_begin(), dg.finite_faces_end());
  test_iterator(dg.all_vertices_begin(), dg.all_vertices_end());
  test_iterator(dg.all_edges_begin(), dg.all_edges_end());
  test_iterator(dg.all_faces_begin(), dg.all_faces_end());

  // testing traversal with circulators
  if ( dg.dimension() == 2 ){
    Vertex_handle v_fin = dg.finite_vertex();
    Face_handle f_adj = dg.incident_faces(v_fin);

    test_circulator(dg.incident_vertices(v_fin));
    test_circulator(dg.incident_vertices(v_fin, f_adj));

    test_circulator(dg.incident_edges(v_fin));
    test_circulator(dg.incident_edges(v_fin, f_adj));

    test_circulator(dg.incident_faces(v_fin));
    test_circulator(dg.incident_faces(v_fin, f_adj));
  }

  // testing predicate methods
  if ( dg.dimension() == 2 ) {
    dg.is_infinite( Vertex_handle(dg.finite_vertex()) );
    dg.is_infinite( Face_handle(dg.infinite_face()) );
    dg.is_infinite( Face_handle(dg.infinite_face()), 0 );
    dg.is_infinite( Edge(dg.infinite_face(), 0) );
    Edge_circulator ec = dg.incident_edges(dg.infinite_vertex());
    dg.is_infinite( ec );
  }

  // testing validity
  assert( dg.is_valid() );
}

//============================================================================
//============================================================================

template<class DG, class VT>
void test_voronoi_traits_concept(const DG& dg, const VT& vt)
{
  typedef typename VT::Point_2                 Point_2;
  typedef typename VT::Site_2                  Site_2;

  typedef typename VT::Dual_graph              Dual_graph;
  typedef typename VT::Voronoi_vertex_2        VV2;
  typedef typename VT::Voronoi_edge_2          VE2;

  typedef typename VT::Edge_degeneracy_tester  EDT;
  typedef typename VT::Face_degeneracy_tester  FDT;
  typedef typename VT::Has_point_locator       Has_pl;

  typedef typename VT::Vertex_handle           Vertex_handle;

  typedef typename VT::Curve                   Curve;

  // test creation of Voronoi edges and vertices
  Vertex_handle v[4];
  switch ( dg.number_of_vertices() ) {
  case 4:
    v[3] = ++(++(++dg.finite_vertices_begin()));
  case 3:
    v[2] = ++(++dg.finite_vertices_begin());
  case 2:
    v[0] = dg.finite_vertices_begin();
    v[1] = ++dg.finite_vertices_begin();
  }

  switch ( dg.number_of_vertices() ) {
  case 4:
    vt.make_edge(v[0], v[1], v[2], v[3]);
  case 3:
    vt.make_vertex(v[0], v[1], v[2]);
    vt.make_edge(v[0], v[1], v[2], true);
    vt.make_edge(v[0], v[1], v[2], false);
  case 2:
    vt.make_edge(v[0], v[1]);
  }

  // test nested concepts
  test_edt_concept( dg, vt.edge_degeneracy_tester_object() );
  test_fdt_concept( dg, vt.face_degeneracy_tester_object() );
  test_pl_concept( dg, vt, Has_pl() );

  if ( dg.dimension() < 2 ) { return; }

  test_segment_ve(dg, vt);
  test_ray_ve(dg, vt);
  test_bisector_ve(dg, vt);

  test_vv_concept( dg, vt );
}


//============================================================================
//============================================================================

template<class DG, class VT>
void test_segment_ve(const DG& dg, const VT& vt)
{
  typedef typename DG::Vertex_handle    Vertex_handle;
  typedef typename DG::Edge             Edge;
  typedef typename VT::Voronoi_edge_2   VE2;
  typedef CGAL::Triangulation_cw_ccw_2  CW_CCW_2;

  CGAL_precondition( dg.dimension() == 2 );

  typename DG::Finite_edges_iterator it = dg.finite_edges_begin();
  bool found = false;
  Edge e_seg;
  do {
    if ( !dg.is_infinite(it->first) &&
	 !dg.is_infinite(it->first->neighbor(it->second)) ) {
      e_seg = *it;
      found = true;
      break;
    }
    ++it;
  } while ( it != dg.finite_edges_end() );

  if ( !found ) { return; }
  
  Vertex_handle v1 = e_seg.first->vertex( CW_CCW_2::cw(e_seg.second) );
  Vertex_handle v2 = e_seg.first->vertex( CW_CCW_2::ccw(e_seg.second) );
  Vertex_handle v3 = e_seg.first->vertex( e_seg.second );
  Vertex_handle v4 = dg.tds().mirror_vertex( e_seg.first, e_seg.second );

  test_ve_concept( vt.make_edge(v1, v2, v3, v4) );
}

//============================================================================

template<class DG, class VT>
void test_ray_ve(const DG& dg, const VT& vt)
{
  typedef typename DG::Vertex_handle    Vertex_handle;
  typedef typename DG::Edge             Edge;
  typedef typename VT::Voronoi_edge_2   VE2;
  typedef CGAL::Triangulation_cw_ccw_2  CW_CCW_2;

  CGAL_precondition( dg.dimension() == 2 );

  typename DG::All_edges_iterator it = dg.all_edges_begin();
  bool found = false;
  Edge e_ray;
  do {
    Vertex_handle v1 = it->first->vertex( it->second );
    Vertex_handle v2 = dg.tds().mirror_vertex( it->first, it->second );
    
    if ( (dg.is_infinite(v1) && !dg.is_infinite(v2)) ||
	 (!dg.is_infinite(v1) && dg.is_infinite(v2)) ) {
      e_ray = *it;
      found = true;
      break;
    }
    ++it;
  } while ( it != dg.all_edges_end() );

  if ( !found ) { return; }

  Vertex_handle v1 = e_ray.first->vertex( CW_CCW_2::cw(e_ray.second) );
  Vertex_handle v2 = e_ray.first->vertex( CW_CCW_2::ccw(e_ray.second) );
  Vertex_handle v3 = e_ray.first->vertex( e_ray.second );
  Vertex_handle v4 = dg.tds().mirror_vertex( e_ray.first, e_ray.second );

  if ( dg.is_infinite(v3) ) {
    test_ve_concept( vt.make_edge(v1, v2, v4, false) );
  } else {
    CGAL_assertion( !dg.is_infinite(v3) );
    test_ve_concept( vt.make_edge(v1, v2, v3, true) );
  }
}

//============================================================================

template<class DG, class VT>
void test_bisector_ve(const DG& dg, const VT& vt)
{
  typedef typename DG::Vertex_handle    Vertex_handle;
  typedef typename DG::Edge             Edge;
  typedef typename VT::Voronoi_edge_2   VE2;
  typedef CGAL::Triangulation_cw_ccw_2  CW_CCW_2;

  CGAL_precondition( dg.dimension() == 2 );

  typename DG::All_edges_iterator it = dg.all_edges_begin();
  bool found = false;
  Edge e_bis;
  do {
    Vertex_handle v1 = it->first->vertex( it->second );
    Vertex_handle v2 = dg.tds().mirror_vertex( it->first, it->second );
    
    if ( (dg.is_infinite(v1) && dg.is_infinite(v2)) ) {
      e_bis = *it;
      found = true;
      break;
    }
    ++it;
  } while ( it != dg.all_edges_end() );

  if ( !found ) { return; }

  Vertex_handle v1 = e_bis.first->vertex( CW_CCW_2::cw(e_bis.second) );
  Vertex_handle v2 = e_bis.first->vertex( CW_CCW_2::ccw(e_bis.second) );
  Vertex_handle v3 = e_bis.first->vertex( e_bis.second );
  Vertex_handle v4 = dg.tds().mirror_vertex( e_bis.first, e_bis.second );

  CGAL_assertion( dg.is_infinite(v3) && dg.is_infinite(v4) );
  test_ve_concept( vt.make_edge(v1, v2) );
}

//============================================================================

template<class VE>
void test_ve_concept(const VE& ve)
{
  typedef typename VE::Dual_graph         Dual_graph;
  typedef typename VE::Site_2             Site_2;
  typedef typename VE::Voronoi_vertex_2   VV2;

  bool b;
  b = ve.is_bisector();
  b = ve.is_ray();
  b = ve.is_segment();
  b = ve.has_source();
  b = ve.has_target();
  kill_warning(b);

  if ( ve.has_source() ) {
    VV2 vv_src = ve.source();
    kill_warning(vv_src);
  }

  if ( ve.has_target() ) {
    VV2 vv_trg = ve.target();
    kill_warning(vv_trg);
  }

  VE ve_opp = ve.opposite();

  Site_2 t = ve.north();
  t = ve.south();
  if ( ve.has_source() ) { t = ve.west(); }
  if ( ve.has_target() ) { t = ve.east(); }
  kill_warning(t);
}

//============================================================================
//============================================================================

template<class DG, class VT>
void test_vv_concept(const DG& dg, const VT& vt)
{
  if ( dg.finite_faces_begin() == dg.finite_faces_end() ) { return; }

  typename DG::Face_handle f = dg.finite_faces_begin();
  typename DG::Vertex_handle v1 = f->vertex(0);
  typename DG::Vertex_handle v2 = f->vertex(1);
  typename DG::Vertex_handle v3 = f->vertex(2);
  test_vv_concept( vt.make_vertex(v1, v2, v3) );
}

//============================================================================

template<class VV>
void test_vv_concept(const VV& vv)
{
  typedef typename VV::Dual_graph         Dual_graph;
  typedef typename VV::Site_2             Site_2;

  typedef typename Dual_graph::Geom_traits::Point_2   Point_2;

  Site_2 t;
  t = vv.first();
  t = vv.second();
  t = vv.third();

  for (unsigned int i = 0; i < 3; i++) {
    t = vv.site(i);
  }

  // testing convertibility to point
  Point_2 p = vv;

  kill_warning(p);
}

//============================================================================
//============================================================================

template<class DG, class EDT>
void test_edt_concept(const DG& dg, const EDT& edt)
{
  // types for AdaptableFunctor concept
  typedef typename EDT::result_type               result_type;
  typedef typename EDT::Arity                     Arity;


  typedef typename EDT::Dual_graph                Dual_graph;
  typedef typename EDT::Edge                      Edge;
  typedef typename EDT::Face_handle               Face_handle;
  typedef typename EDT::Edge_circulator           Edge_circulator;
  typedef typename EDT::All_edges_iterator        All_edges_iterator;
  typedef typename EDT::Finite_edges_iterator     Finite_edges_iterator;

  if ( dg.dimension() < 1 ) { return; }

  bool b;

  if ( dg.all_edges_begin() != dg.all_edges_end() ) {
    b = edt(dg, dg.all_edges_begin());
  }

  if ( dg.finite_edges_begin() != dg.finite_edges_end() ) {
    b = edt(dg, dg.finite_edges_begin());

    Edge e = *dg.finite_edges_begin();
    b = edt(dg, e);
    b = edt(dg, e.first, e.second);
  }

  Edge_circulator ec = dg.incident_edges(dg.finite_vertex());
  b = edt(dg, ec);
  kill_warning(b);
}

//============================================================================
//============================================================================

template<class DG, class FDT>
void test_fdt_concept(const DG& dg, const FDT& fdt)
{
  // types for the AdaptableFunctor concept
  typedef typename FDT::result_type               result_type;
  typedef typename FDT::Arity                     Arity;

  typedef typename FDT::Dual_graph                Dual_graph;
  typedef typename FDT::Vertex_handle             Vertex_handle;

  if ( dg.dimension() < 1 ) { return; }

  Vertex_handle v = dg.finite_vertex();
  bool b = fdt(dg, v);
  kill_warning(b);
}

//============================================================================
//============================================================================

template<class DG, class VT>
void test_pl_concept(const DG&, const VT&, CGAL::Tag_false) {}

template<class DG, class VT>
void test_pl_concept(const DG& dg, const VT& vt, CGAL::Tag_true)
{
  typedef typename VT::Point_locator             Point_locator;
  typedef typename Point_locator::Dual_graph     Dual_graph;
  typedef typename Point_locator::Locate_type    Locate_type;
  typedef typename Point_locator::Vertex_handle  Vertex_handle;
  typedef typename Point_locator::Face_handle    Face_handle;
  typedef typename Point_locator::Edge           Edge;
  typedef typename Point_locator::Point_2        Point_2;

  typedef typename Locate_type::Dual_graph       LT_Dual_graph;
  typedef typename Locate_type::Vertex_handle    LT_Vertex_handle;
  typedef typename Locate_type::Face_handle      LT_Face_handle;
  typedef typename Locate_type::Edge             LT_Edge;

  if ( dg.dimension() < 0 ) { return; }

  Point_locator pl = vt.point_locator_object();
  Point_2 p(0,0);

  Locate_type lt = pl(dg, p);

  if ( lt.is_face() ) {
    LT_Face_handle lt_f = lt.face();
    Face_handle f = lt_f;
    kill_warning(f);
  } else if ( lt.is_edge() ) {
    LT_Edge lt_e = lt.edge();
    Edge e = lt_e;
    kill_warning(e);
  } else if ( lt.is_vertex() ) {
    LT_Vertex_handle lt_v = lt.vertex();
    Vertex_handle v = lt_v;
    kill_warning(v);
  }
}

//============================================================================
//============================================================================

#endif // CGAL_VDA_TEST_CONCEPT_H
