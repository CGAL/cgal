#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/Graph_with_descriptor_with_graph.h>

#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/graph/graph_concepts.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> SM;
typedef CGAL::Graph_with_descriptor_with_graph<SM> Surface_mesh;
typedef boost::graph_traits< Surface_mesh > Traits;
typedef Traits::edge_descriptor edge_descriptor;
typedef Traits::halfedge_descriptor halfedge_descriptor;
typedef Traits::vertex_descriptor vertex_descriptor;
typedef Traits::face_descriptor face_descriptor;

void concept_check_surface_mesh()
{
  boost::function_requires< boost::GraphConcept<Surface_mesh> >();
  boost::function_requires< boost::VertexListGraphConcept<Surface_mesh> >();
  boost::function_requires< boost::EdgeListGraphConcept<Surface_mesh> >();
  boost::function_requires< boost::IncidenceGraphConcept<Surface_mesh> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<Surface_mesh> >();
  boost::function_requires< boost::BidirectionalGraphConcept<Surface_mesh> >();
  boost::function_requires< CGAL::HalfedgeGraphConcept<Surface_mesh> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<Surface_mesh> >();
  boost::function_requires< CGAL::FaceGraphConcept<Surface_mesh> >();
  boost::function_requires< CGAL::FaceListGraphConcept<Surface_mesh> >();
  boost::function_requires< CGAL::MutableHalfedgeGraphConcept<Surface_mesh> >();
  boost::function_requires< CGAL::MutableFaceGraphConcept<Surface_mesh> >();

  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Surface_mesh, halfedge_descriptor, CGAL::halfedge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Surface_mesh, edge_descriptor, boost::edge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Surface_mesh, edge_descriptor, boost::edge_weight_t> >();
  boost::function_requires< boost::concepts::PropertyGraph<
    Surface_mesh, vertex_descriptor, CGAL::vertex_point_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Surface_mesh, vertex_descriptor, boost::vertex_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Surface_mesh, face_descriptor, CGAL::face_index_t> >();

  // null
  boost::graph_traits<Surface_mesh>::null_vertex();
  boost::graph_traits<Surface_mesh>::null_halfedge();
  boost::graph_traits<Surface_mesh>::null_face();
}

int main()
{
  concept_check_surface_mesh();
  return 0;
}
