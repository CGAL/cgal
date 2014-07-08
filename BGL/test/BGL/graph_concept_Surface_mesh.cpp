#define BOOST_TEST_MODULE graph_traits_concepts
#include <boost/test/included/unit_test.hpp>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>

#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/graph/graph_concepts.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Sm;

BOOST_AUTO_TEST_CASE( concepts )
{
  typedef boost::graph_traits< Sm >                        Traits;
  typedef Traits::vertex_descriptor               vertex_descriptor;
  typedef Traits::edge_descriptor                 edge_descriptor;
  typedef Traits::halfedge_descriptor             halfedge_descriptor;
  typedef Traits::face_descriptor                 face_descriptor;
  typedef Traits::face_iterator                   face_iterator;
  typedef Traits::faces_size_type                 faces_size_type;
  typedef Traits::out_edge_iterator               out_edge_iterator;
  typedef Traits::in_edge_iterator                in_edge_iterator;
  typedef Traits::vertex_iterator                 vertex_iterator;
  typedef Traits::vertices_size_type              vertices_size_type;
  typedef Traits::edge_iterator                   edge_iterator;
  typedef Traits::edges_size_type                 edges_size_type;
  typedef Traits::halfedge_iterator               halfedge_iterator;
  typedef Traits::halfedges_size_type             halfedges_size_type;

  boost::function_requires< boost::GraphConcept<Sm> >();
  boost::function_requires< boost::VertexListGraphConcept<Sm> >();
  boost::function_requires< boost::EdgeListGraphConcept<Sm> >();
  boost::function_requires< boost::IncidenceGraphConcept<Sm> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<Sm> >();
  boost::function_requires< boost::BidirectionalGraphConcept<Sm> >();
  boost::function_requires< boost::MutableGraphConcept<Sm> >();

  boost::function_requires< CGAL::HalfedgeGraphConcept<Sm> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<Sm> >();
  boost::function_requires< CGAL::FaceGraphConcept<Sm> >();
  boost::function_requires< CGAL::FaceListGraphConcept<Sm> >();
  boost::function_requires< CGAL::MutableFaceGraphConcept<Sm> >();
  boost::function_requires< CGAL::MutableHalfedgeGraphConcept<Sm> >();

  // this uses some internal boost concepts, better than nothing
  boost::function_requires< boost::concepts::ReadablePropertyGraph<Sm, edge_descriptor, boost::edge_weight_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<Sm, vertex_descriptor, boost::vertex_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<Sm, edge_descriptor, boost::edge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<Sm, face_descriptor, CGAL::face_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<Sm, halfedge_descriptor,
                                                                   CGAL::halfedge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<Sm, vertex_descriptor,
                                                                   CGAL::vertex_is_border_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<Sm, halfedge_descriptor,
                                                                   CGAL::halfedge_is_border_t> >();
  boost::function_requires< boost::concepts::PropertyGraphConcept<Sm, vertex_descriptor,
                                                                  CGAL::vertex_point_t> >();

  // null
  boost::graph_traits<Sm>::null_vertex();
  boost::graph_traits<Sm>::null_face();
}

// trick cgal_test_with_cmake
// int main()
// {
//   return 0;
// }
