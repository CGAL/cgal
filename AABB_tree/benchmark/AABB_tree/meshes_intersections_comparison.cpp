#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/AABB_tree/internal/AABB_traversal_traits.h>

#include <CGAL/box_intersection_d.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <CGAL/AABB_tree/internal/AABB_two_tree_traversal.h>

#include <tbb/tbb.h>

#include <CGAL/Timer.h>
#include <CGAL/Real_timer.h>

#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/boost/graph/generators.h>

#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Aff_transformation_3.h>

#include <CGAL/AABB_meshes_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::face_descriptor vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;

void two_meshes_intersection(std::string fname1, std::string fname2){

  Surface_mesh tm1;
  Surface_mesh tm2;
  CGAL::IO::read_polygon_mesh(fname1, tm1);
  CGAL::IO::read_polygon_mesh(fname2, tm2);
  PMP::triangulate_faces(tm1);
  PMP::triangulate_faces(tm2);

  CGAL::Timer t;
  CGAL::Real_timer rt;
  t.start();
  rt.start();

  std::vector<std::pair<face_descriptor, face_descriptor>> out;

  out.clear();
  PMP::experimental::AABB_two_tree_meshes_intersections(tm1, tm2, std::back_inserter(out), CGAL::parameters::concurrency_tag(CGAL::Parallel_tag()));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "Two tree AABB intersecton time: " << rt.time() << "sec (" << t.time() << "s all cpu)." << std::endl;

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  out.clear();
  PMP::experimental::AABB_two_tree_meshes_intersections(tm1, tm2, std::back_inserter(out));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "Two tree AABB intersecton time (Sequential): " << rt.time() << "sec (" << t.time() << "s all cpu)." << std::endl;

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  out.clear();
  PMP::experimental::mixed_meshes_intersections(tm1, tm2, std::back_inserter(out), CGAL::parameters::concurrency_tag(CGAL::Parallel_tag()));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "Mixed intersecton time: " << rt.time() << "sec (" << t.time() << "s all cpu)." << std::endl;

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  out.clear();
  PMP::experimental::mixed_meshes_intersections(tm1, tm2, std::back_inserter(out));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "Mixed intersecton time (Sequential): " << rt.time() << "sec (" << t.time() << "s all cpu)." << std::endl;

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  out.clear();
  PMP::AABB_meshes_intersections(tm1, tm2, std::back_inserter(out), CGAL::parameters::concurrency_tag(CGAL::Parallel_tag()));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "AABB intersecton time: " << rt.time() << "sec (" << t.time() << "s all cpu)." << std::endl;
  std::sort(out.begin(), out.end(), [](const auto &a, const auto &b){ return a.first<b.first || (a.first == b.first && a.second<b.second); });

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  out.clear();
  PMP::AABB_meshes_intersections(tm1, tm2, std::back_inserter(out));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "AABB intersecton time (Sequential): " << rt.time() << "sec (" << t.time() << "s all cpu)." << std::endl;


  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  out.clear();
  PMP::experimental::box_meshes_intersections(tm1, tm2, std::back_inserter(out), CGAL::parameters::concurrency_tag(CGAL::Parallel_tag()));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "Box intersecton time: " << rt.time() << "sec (" << t.time() << "s all cpu)." << std::endl;
  std::sort(out.begin(), out.end(), [](const auto &a, const auto &b){ return a.first<b.first || (a.first == b.first && a.second<b.second); });

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  out.clear();
  PMP::experimental::box_meshes_intersections(tm1, tm2, std::back_inserter(out));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "Box intersecton time (Sequential): " << rt.time() << "sec (" << t.time() << "s all cpu)." << std::endl;

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();
}

int main(int argc, char* argv[])
{
  two_meshes_intersection((argc>1)?argv[1]:CGAL::data_file_path("meshes/tetrahedron.off"),
                          (argc>2)?argv[2]:CGAL::data_file_path("meshes/beam.off"));
  // two_meshes_intersection((argc>1)?argv[1]:CGAL::data_file_path("meshes/beam.off"),
  //                         (argc>2)?argv[2]:"beam_transformed.off");
  return EXIT_SUCCESS;
}
