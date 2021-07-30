#include <cmath>
#include <vector>
#include <fstream>

#include <CGAL/Real_timer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Kd_tree.h>

using SCF   = CGAL::Simple_cartesian<float>;
using SCD   = CGAL::Simple_cartesian<double>;
using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using Timer = CGAL::Real_timer;

template<typename Kernel>
void bench_random_in_bbox(const std::string filename) {

  std::cout << std::endl;
  std::cout << "---- RANDOM IN BBOX BENCH" << std::endl;
  std::cout << "- filename " << filename << " ... " << std::endl;

  using Point_3 = typename Kernel::Point_3;
  using Surface_mesh = CGAL::Surface_mesh<Point_3>;

  // Read mesh.
  std::cout << std::endl << "* reading input mesh" << std::endl;
  Surface_mesh surface_mesh;
  CGAL::IO::read_polygon_mesh(filename, surface_mesh);
  std::cout << "- number of vertices: " << surface_mesh.number_of_vertices() << std::endl;
  std::cout << "- number of faces: " << surface_mesh.number_of_faces() << std::endl;
  assert(surface_mesh.number_of_vertices() > 0);
  assert(surface_mesh.number_of_faces() > 0);

  // Create a bbox.
  std::cout << std::endl << "* creating loose bbox" << std::endl;

  std::cout << std::endl;
}

int main(int argc, char* argv[]) {

  const std::string filename = (argc > 1 ? argv[1] : "data/pyramid.off");
  bench_random_in_bbox<SCD>(filename);
}
