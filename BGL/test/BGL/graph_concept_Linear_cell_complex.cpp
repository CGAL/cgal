#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/properties_Linear_cell_complex_for_combinatorial_map.h>

#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/graph/graph_concepts.h>
#include <CGAL/boost/graph/Euler_operations.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
            <2, 3, MyTraits>::type LCC;

typedef boost::graph_traits< LCC > Traits;
typedef Traits::edge_descriptor edge_descriptor;
typedef Traits::halfedge_descriptor halfedge_descriptor;
typedef Traits::vertex_descriptor vertex_descriptor;
typedef Traits::face_descriptor face_descriptor;
typedef LCC::Traits::Point Point_3;

template<typename Graph>
void concept_check_polyhedron() {
  boost::function_requires< boost::GraphConcept<LCC> >();
  boost::function_requires< boost::VertexListGraphConcept<LCC> >();
  boost::function_requires< boost::EdgeListGraphConcept<LCC> >();
  boost::function_requires< boost::IncidenceGraphConcept<LCC> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<LCC> >();
  boost::function_requires< boost::BidirectionalGraphConcept<LCC> >();
  boost::function_requires< CGAL::HalfedgeGraphConcept<LCC> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<LCC> >();
  boost::function_requires< CGAL::FaceGraphConcept<LCC> >();
  boost::function_requires< CGAL::FaceListGraphConcept<LCC> >();
  boost::function_requires< CGAL::MutableHalfedgeGraphConcept<LCC> >();
  boost::function_requires< CGAL::MutableFaceGraphConcept<LCC> >();

  boost::function_requires< boost::concepts::PropertyGraph<
    LCC, halfedge_descriptor, boost::halfedge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    LCC, edge_descriptor, boost::edge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
     LCC, edge_descriptor, boost::edge_weight_t> >();
  boost::function_requires< boost::PropertyGraphConcept<
    LCC, vertex_descriptor, boost::vertex_point_t> >();
  boost::function_requires< boost::concepts::PropertyGraph<
    LCC, vertex_descriptor, boost::vertex_index_t> >();
  boost::function_requires< boost::concepts::PropertyGraph<
    LCC, face_descriptor, boost::face_index_t> >();

  // null
  boost::graph_traits<LCC>::null_vertex();
  boost::graph_traits<LCC>::null_halfedge();
  boost::graph_traits<LCC>::null_face();
}

template<typename LCC>
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
  LCC p;
  vertex_descriptor
    u = add_vertex(Point_3(0,2,0), p),
    v = add_vertex(Point_3(2,2,0), p),
    w = add_vertex(Point_3(0,0,0), p),
    x = add_vertex(Point_3(2,0,0), p),
    y = add_vertex(Point_3(1,1,0), p);

  std::vector<vertex_descriptor> face;
  face.push_back(v); face.push_back(u); face.push_back(y);
  CGAL::Euler::add_face(face, p);

  face.clear();
  face.push_back(v); face.push_back(y); face.push_back(x);
  CGAL::Euler::add_face(face, p);

  face.clear();
  face.push_back(x); face.push_back(y); face.push_back(w);
  CGAL::Euler::add_face(face, p);

  face.clear();
  face.push_back(w); face.push_back(y); face.push_back(u);
  CGAL::Euler::add_face(face, p);

  std::cout<<"num_edges(p)"<<num_edges(p)<<std::endl;
  std::cout<<"num_halfedges(p)"<<num_halfedges(p)<<std::endl;
  std::cout<<"num_faces(p)"<<num_faces(p)<<std::endl;

  assert(num_edges(p) ==  8);
  assert(num_halfedges(p) == 16);
  assert(num_faces(p) == 4);
}

int main()
{
  concept_check_polyhedron<LCC>();
  runtime_check_halfedgegraph<LCC>();
  std::cerr << "done\n";
  return 0;
}
