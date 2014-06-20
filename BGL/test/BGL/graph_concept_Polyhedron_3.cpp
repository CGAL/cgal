// Copyright (c) 2012 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Philipp Möller

#define BOOST_TEST_MAIN 1
#include <boost/test/unit_test.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_concepts.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;
typedef boost::graph_traits< Polyhedron > Traits;
typedef Traits::edge_descriptor edge_descriptor;
typedef Traits::halfedge_descriptor halfedge_descriptor;
typedef Traits::vertex_descriptor vertex_descriptor;
typedef Traits::face_descriptor face_descriptor;
typedef Kernel::Point_3 Point_3;

template<typename Graph> 
void concept_check_polyhedron() {
  boost::function_requires< boost::GraphConcept<Polyhedron> >();
  boost::function_requires< boost::VertexListGraphConcept<Polyhedron> >();
  boost::function_requires< boost::EdgeListGraphConcept<Polyhedron> >();
  boost::function_requires< boost::IncidenceGraphConcept<Polyhedron> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<Polyhedron> >();
  boost::function_requires< boost::BidirectionalGraphConcept<Polyhedron> >();
  // no longer the case
  // boost::function_requires< boost::MutableGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::HalfedgeGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::FaceGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::FaceListGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::MutableHalfedgeGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::MutableFaceGraphConcept<Polyhedron> >();

  boost::function_requires< boost::concepts::PropertyGraph<
    Polyhedron, halfedge_descriptor, CGAL::halfedge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Polyhedron, halfedge_descriptor, CGAL::halfedge_is_border_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Polyhedron, edge_descriptor, boost::edge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Polyhedron, edge_descriptor, boost::edge_weight_t> >();
  boost::function_requires< boost::PropertyGraphConcept<
    Polyhedron, vertex_descriptor, CGAL::vertex_point_t> >();
  boost::function_requires< boost::concepts::PropertyGraph<
    Polyhedron, vertex_descriptor, boost::vertex_index_t> >();
  // boost::function_requires< boost::concepts::ReadablePropertyGraph<
  //   Polyhedron, vertex_descriptor, boost::vertex_is_border_t> >();
  boost::function_requires< boost::concepts::PropertyGraph<
    Polyhedron, face_descriptor, CGAL::face_index_t> >();

  // external indices
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Polyhedron, face_descriptor, CGAL::face_external_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Polyhedron, halfedge_descriptor, CGAL::halfedge_external_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Polyhedron, vertex_descriptor, CGAL::vertex_external_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Polyhedron, edge_descriptor, CGAL::edge_external_index_t> >();

  // null
  boost::graph_traits<Polyhedron>::null_vertex();
  boost::graph_traits<Polyhedron>::null_face();
};

template<typename Polyhedron>
void runtime_check_halfedgegraph()
{
  // u        v
  // +--------+
  // |\      /|
  // | \ f2 / |
  // |  \y /  |
  // | f3\/ f1|
  // |   /\   |
  // |  /  \  |
  // | / f4 \ |
  // |/      \|
  // +--------+
  // w        x
  Polyhedron p;
  vertex_descriptor 
    u = add_vertex(Point_3(0,2,0), p),
    v = add_vertex(Point_3(2,2,0), p),
    w = add_vertex(Point_3(0,0,0), p),
    x = add_vertex(Point_3(2,0,0), p),
    y = add_vertex(Point_3(1,1,0), p);
  
  add_edge(v, u, p);
  add_edge(u, y, p);
  add_edge(y, v, p);
  add_face(p);

  add_edge(u, w, p);
  add_edge(w, y, p);
  add_face(p);

  add_edge(w, x, p);
  add_edge(x, y, p);
  add_face(p);

  add_edge(x, v, p);
  add_face(p);


  BOOST_CHECK_EQUAL(num_edges(p), 8);
  BOOST_CHECK_EQUAL(num_halfedges(p), 16);
  BOOST_CHECK_EQUAL(num_faces(p), 4);
}




boost::unit_test::test_suite*
init_unit_test_suite( int, char** const)
{
  boost::unit_test::framework::master_test_suite().
    add( BOOST_TEST_CASE( &concept_check_polyhedron<Polyhedron> ) );

  boost::unit_test::framework::master_test_suite().
    add( BOOST_TEST_CASE( &runtime_check_halfedgegraph<Polyhedron> ) );

    return 0;
}



// int main()
// {
//   return 0;
// }
