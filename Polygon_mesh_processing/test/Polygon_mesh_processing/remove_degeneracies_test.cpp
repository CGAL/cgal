#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

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

  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::edge_descriptor e, edges(mesh))
    CGAL::Polygon_mesh_processing::is_degenerate_edge(e, mesh);
}

void check_triangle_face_degeneracy(const char* fname)
{
  std::ifstream input(fname);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::face_descriptor f, faces(mesh))
    CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(f, mesh);
}

void test_vetices_duplication(const char* fname)
{
  std::ifstream input(fname);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(mesh);
}

void test_vertex_non_manifoldness(const char* fname)
{
  std::ifstream input(fname);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::vertex_descriptor v, vertices(mesh))
  {
    if(CGAL::Polygon_mesh_processing::is_non_manifold_vertex(v, mesh))
      std::cout << "true\n";

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

  double threshold = 0.8;

  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::face_descriptor f, faces(mesh))
  {
    if(CGAL::Polygon_mesh_processing::is_needle_triangle_face(f, mesh, threshold))
      std::cout << "needle\n";
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

  double threshold = 0.8;

  BOOST_FOREACH(typename boost::graph_traits<Surface_mesh>::face_descriptor f, faces(mesh))
  {
    if(CGAL::Polygon_mesh_processing::is_cap_triangle_face(f, mesh, threshold))
      std::cout << "cap\n";
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
  test_vetices_duplication("data_degeneracies/degtri_four.off");
  test_vertex_non_manifoldness("data/non_manifold_vertex.off");
  test_needle("data_degeneracies/needle.off");
  test_cap("data_degeneracies/cap.off");

  return 0;
}
