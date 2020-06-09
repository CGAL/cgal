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
    std::cerr << fname << " is not a valid off file." << std::endl;
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
    std::cerr << fname << " is not a valid off file." << std::endl;
    std::exit(1);
  }

  std::vector<face_descriptor> all_faces(faces(mesh).begin(), faces(mesh).end());
  assert(CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[0], mesh));
  assert(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[1], mesh));
  assert(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[2], mesh));
  assert(!CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(all_faces[3], mesh));

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
    std::cerr << fname << " is not a valid off file." << std::endl;
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

  test_needles_and_caps("data_degeneracies/caps_and_needles.off");

  return EXIT_SUCCESS;
}
