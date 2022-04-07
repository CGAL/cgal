
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Simple_cartesian.h>

#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/graph/graph_concepts.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> SM;
typedef CGAL::Face_filtered_graph<SM, SM::Property_map<boost::graph_traits<SM>::face_descriptor , std::size_t> > Adapter;
typedef boost::graph_traits< Adapter > Traits;
typedef Traits::edge_descriptor edge_descriptor;
typedef Traits::halfedge_descriptor halfedge_descriptor;
typedef Traits::vertex_descriptor vertex_descriptor;
typedef Traits::face_descriptor face_descriptor;

void concept_check_adapter()
{
  boost::function_requires< boost::GraphConcept<Adapter> >();
  boost::function_requires< boost::VertexListGraphConcept<Adapter> >();
  boost::function_requires< boost::EdgeListGraphConcept<Adapter> >();
  boost::function_requires< boost::IncidenceGraphConcept<Adapter> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<Adapter> >();
  boost::function_requires< boost::BidirectionalGraphConcept<Adapter> >();
  boost::function_requires< CGAL::HalfedgeGraphConcept<Adapter> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<Adapter> >();
  boost::function_requires< CGAL::FaceGraphConcept<Adapter> >();
  boost::function_requires< CGAL::FaceListGraphConcept<Adapter> >();

  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Adapter, halfedge_descriptor, CGAL::halfedge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Adapter, edge_descriptor, boost::edge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Adapter, edge_descriptor, boost::edge_weight_t> >();
  boost::function_requires< boost::concepts::PropertyGraph<
    Adapter, vertex_descriptor, CGAL::vertex_point_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Adapter, vertex_descriptor, boost::vertex_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Adapter, face_descriptor, CGAL::face_index_t> >();

  // null
  boost::graph_traits<Adapter>::null_vertex();
  boost::graph_traits<Adapter>::null_halfedge();
  boost::graph_traits<Adapter>::null_face();
}

int main()
{
  concept_check_adapter();
  return 0;
}
