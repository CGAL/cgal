// Copyright (c) 2015 GeometryFactory (France). All rights reserved.
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

#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_concepts.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/boost/graph/Dual.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

template<typename Primal> 
void concept_check_dual() {
  typedef CGAL::Dual<Primal>          Graph;
  typedef boost::graph_traits<Graph>  Traits;
  typedef typename Traits::edge_descriptor     edge_descriptor;
  typedef typename Traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Traits::vertex_descriptor   vertex_descriptor;
  typedef typename Traits::face_descriptor     face_descriptor;

  boost::function_requires< boost::GraphConcept<Graph> >();
  boost::function_requires< boost::VertexListGraphConcept<Graph> >();
  boost::function_requires< boost::EdgeListGraphConcept<Graph> >();
  boost::function_requires< boost::IncidenceGraphConcept<Graph> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<Graph> >();
  // TODO fixme
  // boost::function_requires< boost::BidirectionalGraphConcept<Graph> >();

  boost::function_requires< CGAL::HalfedgeGraphConcept<Graph> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<Graph> >();
  boost::function_requires< CGAL::FaceGraphConcept<Graph> >();
  boost::function_requires< CGAL::FaceListGraphConcept<Graph> >();

  boost::function_requires< boost::concepts::PropertyGraph<
    Graph, halfedge_descriptor, CGAL::halfedge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Graph, edge_descriptor, boost::edge_index_t> >();
  // boost::function_requires< boost::concepts::ReadablePropertyGraph<
  //   Graph, edge_descriptor, boost::edge_weight_t> >();
  // boost::function_requires< boost::PropertyGraphConcept<
  //   Graph, vertex_descriptor, CGAL::vertex_point_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Graph, vertex_descriptor, boost::vertex_index_t> >();
  // boost::function_requires< boost::concepts::ReadablePropertyGraph<
  //   Graph, vertex_descriptor, boost::vertex_is_border_t> >();
  // boost::function_requires< boost::concepts::PropertyGraph<
  //   Graph, face_descriptor, CGAL::face_index_t> >();

  // null
  boost::graph_traits<Graph>::null_vertex();
  boost::graph_traits<Graph>::null_face();
  boost::graph_traits<Graph>::null_halfedge();
}

int
main()
{
  concept_check_dual<Polyhedron>();
  std::cerr << "done\n";
  return 0;
}
