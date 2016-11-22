#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_concepts.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

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
  boost::function_requires< CGAL::HalfedgeGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::FaceGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::FaceListGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::MutableHalfedgeGraphConcept<Polyhedron> >();
  boost::function_requires< CGAL::MutableFaceGraphConcept<Polyhedron> >();

  boost::function_requires< boost::concepts::PropertyGraph<
    Polyhedron, halfedge_descriptor, CGAL::halfedge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Polyhedron, edge_descriptor, boost::edge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Polyhedron, edge_descriptor, boost::edge_weight_t> >();
  boost::function_requires< boost::PropertyGraphConcept<
    Polyhedron, vertex_descriptor, CGAL::vertex_point_t> >();
  boost::function_requires< boost::concepts::PropertyGraph<
    Polyhedron, vertex_descriptor, boost::vertex_index_t> >();
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
}




int
main()
{
  concept_check_polyhedron<Polyhedron>();

  std::cerr << "done\n";
  return 0;
}
