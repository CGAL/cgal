#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Timer.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/boost/graph/generators.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/AABB_meshes_intersections_3.h>

#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;


int main(int argc, char* argv[])
{
  std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("meshes/refined_elephant.off");
  Surface_mesh tmesh;
  CGAL::IO::read_polygon_mesh(fname, tmesh);

  CGAL::Real_timer rt;
  CGAL::Timer t;
  t.start(); rt.start();
  std::vector<std::pair<face_descriptor, face_descriptor>> out;

  std::cout << " Compare AABB_self_intersections with self_intersections on " << fname << std::endl;

  PMP::AABB_self_intersections(tmesh, std::back_inserter(out), CGAL::parameters::concurrency_tag(CGAL::Parallel_tag()));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "AABB self intersecton time: " << rt.time() << "sec (" << t.time() << " s)." << std::endl;

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  out.clear();
  PMP::self_intersections<CGAL::Parallel_tag>(tmesh, std::back_inserter(out));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "self intersecton time: " << rt.time() << "sec (" << t.time() << " s)." << std::endl;

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  out.clear();
  PMP::AABB_self_intersections(tmesh, std::back_inserter(out));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "AABB self intersecton time (Sequential): " << rt.time() << "sec (" << t.time() << " s)." << std::endl;

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  out.clear();
  PMP::self_intersections(tmesh, std::back_inserter(out));
  std::cout << "number intersections: " << out.size() << std::endl;
  std::cout << "self intersecton time (Sequential): " << rt.time() << "sec (" << t.time() << " s)." << std::endl;

  t.stop(); rt.stop();
  t.reset(); rt.reset();
  t.start(); rt.start();

  return EXIT_SUCCESS;
}
