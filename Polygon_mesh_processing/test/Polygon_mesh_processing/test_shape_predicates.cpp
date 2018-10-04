#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

#include <CGAL/number_type_config.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::FT                                                       FT;
typedef K::Point_3                                                  Point_3;
typedef CGAL::Surface_mesh<Point_3>                                 Surface_mesh;

void check_edge_degeneracy(const char* fname)
{
  std::cout << "test edge degeneracy...";

  typedef boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;

  std::ifstream input(fname);
  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }
  std::vector<edge_descriptor> all_edges(edges(mesh).begin(), edges(mesh).end());

  assert(!CGAL::Polygon_mesh_processing::is_degenerate_edge(all_edges[0], mesh));
  assert(!CGAL::Polygon_mesh_processing::is_degenerate_edge(all_edges[1], mesh));
  assert(CGAL::Polygon_mesh_processing::is_degenerate_edge(all_edges[2], mesh));
  std::cout << "done" << std::endl;
}

void check_triangle_face_degeneracy(const char* fname)
{
  std::cout << "test face degeneracy...";

  typedef boost::graph_traits<Surface_mesh>::face_descriptor   face_descriptor;

  std::ifstream input(fname);
  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  std::vector<face_descriptor> all_faces(faces(mesh).begin(), faces(mesh).end());
  assert(CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[0], mesh));
  assert(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[1], mesh));
  assert(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[2], mesh));
  assert(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[3], mesh));

  std::cout << "done" << std::endl;
}

// tests merge_and_duplication
template <typename PolygonMesh>
void merge_identical_points(typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_keep,
                            typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_rm,
                            PolygonMesh& mesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h = halfedge(v_rm, mesh);
  halfedge_descriptor start = h;

  do
  {
    set_target(h, v_keep, mesh);
    h = opposite(next(h, mesh), mesh);
  }
  while( h != start );

  remove_vertex(v_rm, mesh);
}

void test_vertex_non_manifoldness(const char* fname)
{
  std::cout << "test vertex non manifoldness...";

  typedef boost::graph_traits<Surface_mesh>::vertex_descriptor   vertex_descriptor;
  typedef boost::graph_traits<Surface_mesh>::vertices_size_type  size_type;

  std::ifstream input(fname);
  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  size_type ini_nv = num_vertices(mesh);

  // create non-manifold vertex
  Surface_mesh::Vertex_index vertex_to_merge_onto(1);
  Surface_mesh::Vertex_index vertex_to_merge(7);
  merge_identical_points(vertex_to_merge_onto, vertex_to_merge, mesh);
  mesh.collect_garbage();

  assert(num_vertices(mesh) == ini_nv - 1);

  BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
  {
    if(v == vertex_to_merge_onto)
      assert(CGAL::Polygon_mesh_processing::is_non_manifold_vertex(v, mesh));
    else
      assert(!CGAL::Polygon_mesh_processing::is_non_manifold_vertex(v, mesh));
  }

  std::cout << "done" << std::endl;
}

void test_vertices_merge_and_duplication(const char* fname)
{
  std::cout << "test non manifold vertex duplication...";

  typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;

  std::ifstream input(fname);
  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }
  const std::size_t initial_vertices = num_vertices(mesh);

  // create non-manifold vertex
  Surface_mesh::Vertex_index vertex_to_merge_onto(1);
  Surface_mesh::Vertex_index vertex_to_merge(7);
  Surface_mesh::Vertex_index vertex_to_merge_2(14);
  Surface_mesh::Vertex_index vertex_to_merge_3(21);

  Surface_mesh::Vertex_index vertex_to_merge_onto_2(2);
  Surface_mesh::Vertex_index vertex_to_merge_4(8);

  merge_identical_points(vertex_to_merge_onto, vertex_to_merge, mesh);
  merge_identical_points(vertex_to_merge_onto, vertex_to_merge_2, mesh);
  merge_identical_points(vertex_to_merge_onto, vertex_to_merge_3, mesh);
  merge_identical_points(vertex_to_merge_onto_2, vertex_to_merge_4, mesh);
  mesh.collect_garbage();

  const std::size_t vertices_after_merge = num_vertices(mesh);
  assert(vertices_after_merge == initial_vertices - 4);

  std::vector<std::vector<vertex_descriptor> > duplicated_vertices;
  std::size_t new_vertices_nb =
    CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(mesh,
      CGAL::parameters::output_iterator(std::back_inserter(duplicated_vertices)));

  const std::size_t final_vertices_size = vertices(mesh).size();
  assert(final_vertices_size == initial_vertices);
  assert(new_vertices_nb == 4);
  assert(duplicated_vertices.size() == 2); // two non-manifold vertex
  assert(duplicated_vertices.front().size() == 4);
  assert(duplicated_vertices.back().size() == 2);

  std::cout << "done" << std::endl;
}

void test_needles_and_caps(const char* fname)
{
  std::cout << "test needles&caps...";

  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor     halfedge_descriptor;
  typedef boost::graph_traits<Surface_mesh>::face_iterator           face_iterator;
  typedef boost::graph_traits<Surface_mesh>::face_descriptor         face_descriptor;

  std::ifstream input(fname);
  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }

  const FT eps = std::numeric_limits<FT>::epsilon();

  face_iterator fit, fend;
  boost::tie(fit, fend) = faces(mesh);

  // (0 0 0) -- (1 0 0) -- (1 1 0) (90° cap angle)
  face_descriptor f = *fit;
  halfedge_descriptor res = PMP::is_needle_triangle_face(f, mesh, 2/*needle_threshold*/);
  assert(res == boost::graph_traits<Surface_mesh>::null_halfedge()); // not a needle
  res = PMP::is_needle_triangle_face(f, mesh, CGAL::sqrt(FT(2) - eps)/*needle_threshold*/);
  assert(res != boost::graph_traits<Surface_mesh>::null_halfedge()); // is a needle

  res = PMP::is_cap_triangle_face(f, mesh, 0./*cos(pi/2)*/);
  assert(mesh.point(target(res, mesh)) == CGAL::ORIGIN);
  res = PMP::is_cap_triangle_face(f, mesh, std::cos(91 * CGAL_PI / 180));
  assert(res == boost::graph_traits<Surface_mesh>::null_halfedge());
  res = PMP::is_cap_triangle_face(f, mesh, std::cos(2 * CGAL_PI / 3));
  assert(res == boost::graph_traits<Surface_mesh>::null_halfedge());
  ++ fit;

  // (0 0 1) -- (1 0 1) -- (10 10 1)
  f = *fit;
  res = PMP::is_needle_triangle_face(f, mesh, 20);
  assert(res == boost::graph_traits<Surface_mesh>::null_halfedge());
  res = PMP::is_needle_triangle_face(f, mesh, 10 * CGAL::sqrt(FT(2) - eps));
  assert(mesh.point(target(res, mesh)) == Point_3(1,0,1));
  res = PMP::is_needle_triangle_face(f, mesh, 1);
  assert(mesh.point(target(res, mesh)) == Point_3(1,0,1));

  res = PMP::is_cap_triangle_face(f, mesh, 0./*cos(pi/2)*/);
  assert(mesh.point(target(res, mesh)) == Point_3(0,0,1));
  res = PMP::is_cap_triangle_face(f, mesh, std::cos(2 * CGAL_PI / 3));
  assert(mesh.point(target(res, mesh)) == Point_3(0,0,1));
  res = PMP::is_cap_triangle_face(f, mesh, std::cos(0.75 * CGAL_PI));
  assert(res == boost::graph_traits<Surface_mesh>::null_halfedge());
  ++ fit;

  // (0 0 2) -- (1 0 2) -- (-0.99619469809 0.08715574274 2) (175° cap angle)
  f = *fit;
  res = PMP::is_needle_triangle_face(f, mesh, 2);
  assert(res == boost::graph_traits<Surface_mesh>::null_halfedge());
  res = PMP::is_needle_triangle_face(f, mesh, 1.9);
  assert(mesh.point(target(res, mesh)) == Point_3(0,0,2) ||
         mesh.point(target(res, mesh)) == Point_3(1,0,2));
  res = PMP::is_needle_triangle_face(f, mesh, 1);
  assert(mesh.point(target(res, mesh)) == Point_3(0,0,2) ||
         mesh.point(target(res, mesh)) == Point_3(1,0,2));

  res = PMP::is_cap_triangle_face(f, mesh, 0./*cos(pi/2)*/);
  assert(res != boost::graph_traits<Surface_mesh>::null_halfedge() &&
         mesh.point(target(res, mesh)) != Point_3(0,0,2) &&
         mesh.point(target(res, mesh)) != Point_3(1,0,2));
  res = PMP::is_cap_triangle_face(f, mesh, std::cos(2 * CGAL_PI / 3));
  assert(res != boost::graph_traits<Surface_mesh>::null_halfedge());
  res = PMP::is_cap_triangle_face(f, mesh, std::cos(175 * CGAL_PI / 180));
  assert(res != boost::graph_traits<Surface_mesh>::null_halfedge());
  res = PMP::is_cap_triangle_face(f, mesh, std::cos(176 * CGAL_PI / 180));
  assert(res == boost::graph_traits<Surface_mesh>::null_halfedge());

  std::cout << "done" << std::endl;
}

int main()
{
  check_edge_degeneracy("data_degeneracies/degtri_edge.off");
  check_triangle_face_degeneracy("data_degeneracies/degtri_four.off");
  test_vertex_non_manifoldness("data/blobby.off");
  test_vertices_merge_and_duplication("data/blobby.off");
  test_needles_and_caps("data_degeneracies/caps_and_needles.off");

  return EXIT_SUCCESS;
}
