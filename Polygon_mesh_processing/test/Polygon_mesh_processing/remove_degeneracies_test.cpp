#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/stitch_holes.h>

#include <iostream>
#include <fstream>

//note : when
//CGAL::get_default_random()::get_seed() = 1473902576
//the last test (on trihole.off) does not terminate
//

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;

void fix(const char* fname)
{
  std::ifstream input(fname);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }
  CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh);

  assert( CGAL::is_valid_polygon_mesh(mesh) );
}

void check_edge_degeneracy(const char* fname)
{
  std::ifstream input(fname);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }
  typedef typename boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;
  std::vector<edge_descriptor> all_edges(edges(mesh).begin(), edges(mesh).end());

  CGAL_assertion(!CGAL::Polygon_mesh_processing::is_degenerate_edge(all_edges[0], mesh));
  CGAL_assertion(!CGAL::Polygon_mesh_processing::is_degenerate_edge(all_edges[1], mesh));
  CGAL_assertion(CGAL::Polygon_mesh_processing::is_degenerate_edge(all_edges[2], mesh));
}

void check_triangle_face_degeneracy(const char* fname)
{
  std::ifstream input(fname);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  typedef typename boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;
  std::vector<face_descriptor> all_faces(faces(mesh).begin(), faces(mesh).end());
  CGAL_assertion(CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[0], mesh));
  CGAL_assertion(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[1], mesh));
  CGAL_assertion(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[2], mesh));
  CGAL_assertion(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[3], mesh));
}

void test_vertices_merge_and_duplication(const char* fname)
{
  std::ifstream input(fname);
  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }
  const std::size_t initial_vertices = vertices(mesh).size();

  // create non-manifold vertex
  typedef typename boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
  std::vector<vertex_descriptor> all_vertices(vertices(mesh).begin(), vertices(mesh).end());
  CGAL::Polygon_mesh_processing::merge_identical_points(mesh, all_vertices[1], all_vertices[7]);

  const std::size_t vertices_after_merge = vertices(mesh).size();
  CGAL_assertion(vertices_after_merge == initial_vertices - 1);

  CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(mesh);
  const std::size_t final_vertices = vertices(mesh).size();
  CGAL_assertion(final_vertices == vertices_after_merge + 1);
  CGAL_assertion(final_vertices == initial_vertices);
}

void test_vertex_non_manifoldness(const char* fname)
{
  std::ifstream input(fname);
  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  // create non-manifold vertex
  typedef typename boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
  std::vector<vertex_descriptor> all_vertices(vertices(mesh).begin(), vertices(mesh).end());
  CGAL::Polygon_mesh_processing::merge_identical_points(mesh, all_vertices[1], all_vertices[7]);
  std::vector<vertex_descriptor> vertices_with_non_manifold(vertices(mesh).begin(), vertices(mesh).end());
  CGAL_assertion(vertices_with_non_manifold.size() == all_vertices.size() - 1);

  BOOST_FOREACH(std::size_t iv, vertices(mesh))
  {
    vertex_descriptor v = vertices_with_non_manifold[iv];
    if(iv == 1)
      CGAL_assertion(CGAL::Polygon_mesh_processing::is_non_manifold_vertex(v, mesh));
    else
      CGAL_assertion(!CGAL::Polygon_mesh_processing::is_non_manifold_vertex(v, mesh));
  }
}

void test_needle(const char* fname)
{
  std::ifstream input(fname);
  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  const double threshold = 0.8;
  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::face_descriptor f, faces(mesh))
  {
    CGAL_assertion(CGAL::Polygon_mesh_processing::is_needle_triangle_face(f, mesh, threshold));
  }
}

void test_cap(const char* fname)
{
  std::ifstream input(fname);
  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  const double threshold = 0.8;
  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::face_descriptor f, faces(mesh))
  {
    CGAL_assertion(CGAL::Polygon_mesh_processing::is_cap_triangle_face(f, mesh, threshold));
  }
}

int main()
{
  fix("data_degeneracies/degtri_2dt_1edge_split_twice.off");
  fix("data_degeneracies/degtri_four-2.off");
  fix("data_degeneracies/degtri_four.off");
  fix("data_degeneracies/degtri_on_border.off");
  fix("data_degeneracies/degtri_three.off");
  fix("data_degeneracies/degtri_single.off");
  fix("data_degeneracies/trihole.off");
  check_edge_degeneracy("data_degeneracies/degtri_edge.off");
  check_triangle_face_degeneracy("data_degeneracies/degtri_four.off");
  test_vertices_merge_and_duplication("data_degeneracies/non_manifold_vertex_duplicated.off");
  test_vertex_non_manifoldness("data_degeneracies/non_manifold_vertex_duplicated.off");
  test_needle("data_degeneracies/needle.off");
  test_cap("data_degeneracies/cap.off");

  return 0;
}
