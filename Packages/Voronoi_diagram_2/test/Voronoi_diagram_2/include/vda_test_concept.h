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

template<class DG, class VT>
void test_dual_graph_concept(const DG& dg, const VT& vt)
{
  // testing types
  typedef typename DG::size_type                      size_type;
  typedef typename DG::Geom_traits                    Geom_traits;
  typedef typename DG::Triangulation_data_structure   Tds;

  //  typedef typename DG::Site_2                         Site_2;

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

  // testing constructors
  {
    // constructor taking geom traits
    DG dg0(Geom_traits());

    // copy constructor
    DG dg2(dg);

    // assignment operator
    DG dg3;
    dg3 = dg2;

    // constructors that take an iterator range
    std::vector<typename VT::Site_2> v;
    for (Finite_vertices_iterator vit = dg.finite_vertices_begin();
	 vit != dg.finite_vertices_end(); ++vit) {
      v.push_back(vt.get_site_2_object()(vit));
    }
#ifndef CGAL_TRIANGULATION_HIERARCHY_2_H
    // HACK UNTIL Triangulation_hierarchy_2 conforms with the concept
    DG dg4(v.begin(), v.end());
    DG dg5(v.begin(), v.end(), Geom_traits());
#endif
  }

  // testing convertibility of iterators and circulators to handles;
  // testing value type coherence
  if ( dg.number_of_vertices() > 0 ) {
    test_is_convertible_to<Vertex_handle>(dg.all_vertices_begin());
    test_is_convertible_to<Vertex_handle>(dg.finite_vertices_begin());
    test_is_convertible_to<Vertex>(*dg.all_vertices_begin());
    test_is_convertible_to<Vertex>(*dg.finite_vertices_begin());
  }

  if ( dg.number_of_faces() > 0 ) {
    test_is_convertible_to<Face_handle>(dg.all_faces_begin());
    test_is_convertible_to<Face_handle>(dg.finite_faces_begin());
    test_is_convertible_to<Face>(*dg.all_faces_begin());
    test_is_convertible_to<Face>(*dg.finite_faces_begin());
  }

  if ( dg.tds().number_of_edges() > 0 ) {
    test_is_convertible_to<Edge>(*dg.all_edges_begin());
    test_is_convertible_to<Edge>(*dg.finite_edges_begin());    
  }

  if ( dg.dimension() == 2 ) {
    Vertex_circulator vc = dg.incident_vertices(dg.finite_vertex());
    test_is_convertible_to<Vertex_handle>(vc);
    test_is_convertible_to<Vertex>(*vc);

    Face_circulator fc = dg.incident_faces(dg.finite_vertex());
    test_is_convertible_to<Face_handle>(fc);
    test_is_convertible_to<Face>(*fc);

    Edge_circulator ec = dg.incident_edges(dg.finite_vertex());
    test_is_convertible_to<Edge>(*ec);
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

  // testing insertion
  if ( dg.number_of_vertices() > 0 ) {
    DG dg1;
    std::vector<typename VT::Site_2> v;
    for (Finite_vertices_iterator vit = dg.finite_vertices_begin();
	 vit != dg.finite_vertices_end(); ++vit) {
      v.push_back(vt.get_site_2_object()(vit));
    }
    CGAL_assertion( v.size() > 0 );

    // one-site insertion
    typename VT::Site_2 s = v[0];
    dg1.insert(s);

    // insertion that take an iterator range
    dg1.insert(v.begin(), v.end());
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

  typedef typename VT::Delaunay_graph          Delaunay_graph;

  typedef typename VT::Edge_degeneracy_tester  EDT;
  typedef typename VT::Face_degeneracy_tester  FDT;
  typedef typename VT::Has_nearest_site_2      Has_ns;

  typedef typename VT::Vertex_handle           Vertex_handle;

  // test nested concepts
  test_edt_concept( dg, vt.edge_degeneracy_tester_object() );
  test_fdt_concept( dg, vt.face_degeneracy_tester_object() );
  test_ns_concept( dg, vt, Has_ns() );
}


//============================================================================
//============================================================================

template<class DG, class EDT>
void test_edt_concept(const DG& dg, const EDT& edt)
{
  // types for AdaptableFunctor concept
  typedef typename EDT::result_type               result_type;
  typedef typename EDT::Arity                     Arity;


  typedef typename EDT::Delaunay_graph            Delaunay_graph;
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

  typedef typename FDT::Delaunay_graph            Delaunay_graph;
  typedef typename FDT::Vertex_handle             Vertex_handle;

  if ( dg.dimension() < 1 ) { return; }

  Vertex_handle v = dg.finite_vertex();
  bool b = fdt(dg, v);
  kill_warning(b);
}

//============================================================================
//============================================================================

template<class DG, class VT>
void test_ns_concept(const DG&, const VT&, CGAL::Tag_false) {}

template<class DG, class VT>
void test_ns_concept(const DG& dg, const VT& vt, CGAL::Tag_true)
{
  typedef typename VT::Nearest_site_2             Nearest_site_2;
  typedef typename Nearest_site_2::Delaunay_graph Delaunay_graph;
  typedef typename Nearest_site_2::Query_result   Query_result;
  typedef typename Nearest_site_2::Point_2        Point_2;

  typedef typename Query_result::Delaunay_graph   QR_Delaunay_graph;
  typedef typename Query_result::Vertex_handle    QR_Vertex_handle;
  typedef typename Query_result::Face_handle      QR_Face_handle;
  typedef typename Query_result::Edge             QR_Edge;

  if ( dg.dimension() < 0 ) { return; }

  Nearest_site_2 ns = vt.nearest_site_2_object();
  Point_2 p(0,0);

  Query_result qr = ns(dg, p);

  if ( qr.is_face() ) {
    QR_Face_handle f = qr;
    kill_warning(f);
  } else if ( qr.is_edge() ) {
    QR_Edge e = qr;
    kill_warning(e);
  } else if ( qr.is_vertex() ) {
    QR_Vertex_handle v = qr;
    kill_warning(v);
  }

  Query_result qr1 = ns(dg, p);
  bool b = (qr == qr1);
  CGAL_assertion(b);
  kill_warning(b);
}

//============================================================================
//============================================================================

#endif // CGAL_VDA_TEST_CONCEPT_H
