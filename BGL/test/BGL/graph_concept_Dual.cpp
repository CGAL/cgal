#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
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

  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Graph, edge_descriptor, boost::edge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Graph, vertex_descriptor, boost::vertex_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Graph, face_descriptor, boost::face_index_t> >();

  // edge properties should be forwarded
  //boost::function_requires< boost::concepts::ReadablePropertyGraph<
  //  Graph, edge_descriptor, boost::edge_weight_t> >();

  // boost::function_requires< boost::concepts::PropertyGraph<
  //   Graph, halfedge_descriptor, CGAL::halfedge_index_t> >();
  // boost::function_requires< boost::concepts::ReadablePropertyGraph<
  //   Graph, edge_descriptor, boost::edge_weight_t> >();
  // boost::function_requires< boost::PropertyGraphConcept<
  //   Graph, vertex_descriptor, CGAL::vertex_point_t> >();
  // boost::function_requires< boost::concepts::ReadablePropertyGraph<
  //   Graph, vertex_descriptor, boost::vertex_is_border_t> >();

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
