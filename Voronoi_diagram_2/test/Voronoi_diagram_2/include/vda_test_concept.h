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

#ifndef CGAL_VDA_TEST_CONCEPT_H
#define CGAL_VDA_TEST_CONCEPT_H 1

#include <CGAL/Triangulation_utils_2.h>
#include <cassert>
#include "helper_functions.h"
#include <list>
#include <boost/variant.hpp>


template<class DG, class AS>
void test_as_concept(const DG& dg, const AS& as);

template<class DG, class CVP>
void test_cvp_concept(const DG& dg, const CVP& cvp);

template<class DG, class AT>
void test_ns_concept(const DG&, const AT&, CGAL::Tag_false);

template<class DG, class AT>
void test_ns_concept(const DG& dg, const AT& at, CGAL::Tag_true);

template<class DG, class ER>
void test_er_concept(const DG& dg, const ER& er);

template<class DG, class FR>
void test_fr_concept(const DG& dg, const FR& fr);

template<class DG, class AT, class AP>
  void test_si_concept(const DG&, const AT&, const AP&, CGAL::Tag_false);

template<class DG, class AT, class AP>
void test_si_concept(const DG& dg, const AT& at, const AP& ap, CGAL::Tag_true);


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

template<class DG, class AT>
void test_dual_graph_concept(const DG& dg, const AT& at)
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

  CGAL_USE_TYPE(typename DG::All_edges_iterator);
  CGAL_USE_TYPE(typename DG::Finite_edges_iterator);

  CGAL_USE_TYPE(typename DG::All_vertices_iterator);
  typedef typename DG::Finite_vertices_iterator       Finite_vertices_iterator;

  CGAL_USE_TYPE(typename DG::All_faces_iterator);
  CGAL_USE_TYPE(typename DG::Finite_faces_iterator);

  typedef typename DG::Face_circulator                Face_circulator;
  typedef typename DG::Edge_circulator                Edge_circulator;
  typedef typename DG::Vertex_circulator              Vertex_circulator;

  // testing constructors
  {
    // constructor taking geom traits
    Geom_traits gt;
    DG dg0(gt);

    // copy constructor
    DG dg2(dg);

    // assignment operator
    DG dg3;
    dg3 = dg2;

    // constructors that take an iterator range
    std::vector<typename AT::Site_2> v;
    for (Finite_vertices_iterator vit = dg.finite_vertices_begin();
         vit != dg.finite_vertices_end(); ++vit) {
      v.push_back(at.access_site_2_object()(vit));
    }
    DG dg4(v.begin(), v.end());
    DG dg5(v.begin(), v.end(), Geom_traits());
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
    std::vector<typename AT::Site_2> v;
    for (Finite_vertices_iterator vit = dg.finite_vertices_begin();
         vit != dg.finite_vertices_end(); ++vit) {
      v.push_back(at.access_site_2_object()(vit));
    }
    assert( v.size() > 0 );

    // one-site insertion
    typename AT::Site_2 s = v[0];
    dg1.insert(s);

    // insertion that take an iterator range
    dg1.insert(v.begin(), v.end());
  }

  // testing validity
  assert( dg.is_valid() );
}

//============================================================================
//============================================================================

template<class DG, class AT>
void test_adaptation_traits_concept(const DG& dg, const AT& at)
{
  CGAL_USE_TYPE(typename AT::Point_2);
  CGAL_USE_TYPE(typename AT::Site_2);
  CGAL_USE_TYPE(typename AT::Delaunay_graph);
  CGAL_USE_TYPE(typename AT::Delaunay_vertex_handle);
  CGAL_USE_TYPE(typename AT::Delaunay_face_handle);
  CGAL_USE_TYPE(typename AT::Delaunay_edge);
  //  typedef typename AT::Edge_degeneracy_tester     EDT;
  //  typedef typename AT::Face_degeneracy_tester     FDT;
  CGAL_USE_TYPE(typename AT::Access_site_2);
  CGAL_USE_TYPE(typename AT::Construct_Voronoi_point_2);
  typedef typename AT::Has_nearest_site_2         Has_ns;
  //  typedef typename AT::Has_site_inserter          Has_si;

  // testing copy constructor and assignment operator
  {
    AT at2(at);
    kill_warning(at2);
    //    assert( at2.is_valid() );

    AT at3;
    at3 = at;
    kill_warning(at3);
    //    assert( at3.is_valid() );
  }

  // testing clear and swap methods
#if 0
  {
    AT at2(at);
    at2.clear();
    assert( at2.is_valid() );

    AT at3(at);
    at2.swap(at3);

    assert( at2.is_valid() );
    assert( at3.is_valid() );
  }

  // testing validity method
  bool b = at.is_valid();
  kill_warning( b );
  assert( b );
#endif

  // test nested concepts
  //  test_edt_concept( dg, at.edge_degeneracy_tester_object() );
  //  test_fdt_concept( dg, at.face_degeneracy_tester_object() );
  test_as_concept( dg, at.access_site_2_object() );
  test_cvp_concept( dg, at.construct_Voronoi_point_2_object() );
  test_ns_concept( dg, at, Has_ns() );
  //  test_si_concept( dg, at, Has_si() );
}


//============================================================================
//============================================================================

template<class DG, class AT, class AP>
void test_adaptation_policy_concept(const DG& dg, const AT& at, const AP& ap)
{
  CGAL_USE_TYPE(typename AP::Site_2);
  CGAL_USE_TYPE(typename AP::Delaunay_graph);
  CGAL_USE_TYPE(typename AP::Delaunay_vertex_handle);
  CGAL_USE_TYPE(typename AP::Delaunay_face_handle);
  CGAL_USE_TYPE(typename AP::Delaunay_edge);
  CGAL_USE_TYPE(typename AP::All_Delaunay_edges_iterator);
  CGAL_USE_TYPE(typename AP::Finite_Delaunay_edges_iterator);
  CGAL_USE_TYPE(typename AP::Delaunay_edge_circulator);
  CGAL_USE_TYPE(typename AP::Edge_rejector);
  CGAL_USE_TYPE(typename AP::Face_rejector);
  typedef typename AP::Has_site_inserter               Has_si;

  // testing copy constructor and assignment operator
  {
    AP ap2(ap);
    assert( ap2.is_valid() );

    AP ap3;
    ap3 = ap;
    assert( ap3.is_valid() );
  }

  // testing clear and swap methods
  {
    AP ap2(ap);
    ap2.clear();
    assert( ap2.is_valid() );
    assert( ap2.is_valid(dg) );

    AP ap3(ap);
    ap2.swap(ap3);

    assert( ap2.is_valid() );
    assert( ap2.is_valid(dg) );
    assert( ap3.is_valid() );
    assert( ap3.is_valid(dg) );
  }

  // testing validity method
  bool b = ap.is_valid();
  b = ap.is_valid(dg);
  kill_warning( b );
  assert( b );

  // test nested concepts
  test_er_concept( dg, ap.edge_rejector_object() );
  test_fr_concept( dg, ap.face_rejector_object() );
  //  test_as_concept( dg, at.access_site_2_object() );
  //  test_cvp_concept( dg, at.construct_Voronoi_point_2_object() );
  //  test_ns_concept( dg, at, Has_ns() );
  test_si_concept( dg, at, ap, Has_si() );
}


//============================================================================
//============================================================================

template<class DG, class ER>
void test_er_concept(const DG& dg, const ER& er)
{
  // types for AdaptableFunctor concept
  CGAL_USE_TYPE(typename ER::result_type);

  CGAL_USE_TYPE(typename ER::Delaunay_graph);
  typedef typename ER::Edge                      Edge;
  CGAL_USE_TYPE(typename ER::Face_handle);
  typedef typename ER::Edge_circulator           Edge_circulator;
  CGAL_USE_TYPE(typename ER::All_edges_iterator);
  CGAL_USE_TYPE(typename ER::Finite_edges_iterator);

  if ( dg.dimension() < 1 ) { return; }

  bool b;

  if ( dg.all_edges_begin() != dg.all_edges_end() ) {
    b = er(dg, dg.all_edges_begin());
  }

  if ( dg.finite_edges_begin() != dg.finite_edges_end() ) {
    b = er(dg, dg.finite_edges_begin());

    Edge e = *dg.finite_edges_begin();
    b = er(dg, e);
    b = er(dg, e.first, e.second);
  }

  Edge_circulator ec = dg.incident_edges(dg.finite_vertex());
  b = er(dg, ec);
  kill_warning(b);

  // testing copy constructor and assignment operator
  {
    ER er2(er);
    ER er3;
    er3 = er;
    kill_warning(er2);
    kill_warning(er3);
  }

  // testing clear and swap
  {
    ER er2(er);
    ER er3(er);
    er2.clear();
    er2.swap(er3);
    kill_warning(er2);
    kill_warning(er3);
  }

  // validity testing
  b = er.is_valid();
  kill_warning(b);
}

//============================================================================
//============================================================================

template<class DG, class FR>
void test_fr_concept(const DG& dg, const FR& fr)
{
  // types for the AdaptableFunctor concept
  CGAL_USE_TYPE(typename FR::result_type);

  CGAL_USE_TYPE(typename FR::Delaunay_graph);
  typedef typename FR::Vertex_handle             Vertex_handle;

  if ( dg.dimension() < 1 ) { return; }

  Vertex_handle v = dg.finite_vertex();
  bool b = fr(dg, v);
  kill_warning(b);

  // testing copy constructor and assignment operator
  {
    FR fr2(fr);
    FR fr3;
    fr3 = fr;
    kill_warning(fr2);
    kill_warning(fr3);
  }

  // testing clear and swap
  {
    FR fr2(fr);
    FR fr3(fr);
    fr2.clear();
    fr2.swap(fr3);
    kill_warning(fr2);
    kill_warning(fr3);
  }

  // validity testing
  b = fr.is_valid();
  kill_warning(b);
}

//============================================================================
//============================================================================

template<class DG, class AS>
void test_as_concept(const DG& dg, const AS& as)
{
  // types for the AdaptableFunctor concept
  typedef typename AS::result_type               result_type;

  CGAL_USE_TYPE(typename AS::Vertex_handle);

  if ( dg.number_of_vertices() > 0 ) {
    result_type site = as(dg.finite_vertices_begin());
    kill_warning( site );
  }
}

//============================================================================
//============================================================================

template<class DG, class CVP>
void test_cvp_concept(const DG& dg, const CVP& cvp)
{
  // types for the AdaptableFunctor concept
  typedef typename CVP::result_type               result_type;

  typedef typename CVP::Face_handle               Face_handle;

  if ( dg.dimension() > 1 &&
       dg.finite_faces_begin() != dg.finite_faces_end() ) {
    Face_handle f = dg.finite_faces_begin();
    result_type point = cvp(f);
    kill_warning( point );
  }
}

//============================================================================
//============================================================================

template<class DG, class AT, class AP>
void test_si_concept(const DG&, const AT&, const AP&, CGAL::Tag_false) {}

template<class DG, class AT, class AP>
void test_si_concept(const DG& dg, const AT& at, const AP& ap, CGAL::Tag_true)
{
  typedef typename AP::Site_inserter               Site_inserter;
  CGAL_USE_TYPE(typename Site_inserter::Delaunay_graph);
  typedef typename Site_inserter::Site_2           Site_2;

  // types for the AdaptableFunctor concept
  typedef typename Site_inserter::result_type      result_type;

  Site_inserter si = ap.site_inserter_object();

  if ( dg.number_of_vertices() == 0 ) { return; }

  DG dg2;

  assert( dg2.number_of_vertices() == 0 );

  Site_2 t = at.access_site_2_object()(dg.finite_vertices_begin());
  result_type v = si(dg2, t);
  kill_warning( v );

  assert( dg2.number_of_vertices() == 1 );

  dg2.clear();
  assert( dg2.number_of_vertices() == 0 );

  typedef std::list<result_type>               res_list;
  typedef std::back_insert_iterator<res_list>  output_iterator;

  res_list v_list;

  assert( v_list.size() == 0 );

  output_iterator oit(v_list);
  oit = si(dg2, t, oit);

  assert( v_list.size() == 1 );
  assert( dg2.number_of_vertices() == 1 );
}

//============================================================================
//============================================================================

template<class DG, class AT>
void test_ns_concept(const DG&, const AT&, CGAL::Tag_false) {}

template<class DG, class AT>
void test_ns_concept(const DG& dg, const AT& at, CGAL::Tag_true)
{
  typedef typename AT::Nearest_site_2             Nearest_site_2;
  CGAL_USE_TYPE(typename Nearest_site_2::Delaunay_graph);
  typedef typename Nearest_site_2::Point_2        Point_2;

  // types for the AdaptableFunctor concept
  typedef typename Nearest_site_2::result_type    result_type;

  typedef typename AT::Delaunay_vertex_handle     Vertex_handle;
  typedef typename AT::Delaunay_face_handle       Face_handle;
  typedef typename AT::Delaunay_edge              Edge;

  if ( dg.dimension() < 0 ) { return; }

  Nearest_site_2 ns = at.nearest_site_2_object();
  Point_2 p(0,0);

  result_type qr = ns(dg, p);

  if ( Face_handle* f = boost::get<Face_handle>(&qr) ) {
    kill_warning(f);
  } else if ( Edge* e = boost::get<Edge>(&qr) ) {
    kill_warning(e);
  } else if ( Vertex_handle* v = boost::get<Vertex_handle>(&qr) ) {
    kill_warning(v);
  } else {
    // we should have reached this line
    CGAL_error();
  }

  result_type qr1 = ns(dg, p);
  bool b = (qr == qr1);
  assert(b);
  kill_warning(b);
}

//============================================================================
//============================================================================

#endif // CGAL_VDA_TEST_CONCEPT_H
