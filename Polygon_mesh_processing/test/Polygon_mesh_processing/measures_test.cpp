#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Real_timer.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <list>


namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;

template <typename Mesh>
void
test(std::string fname)
{
  std::cout << typeid(Mesh).name() << std::endl;
  Mesh pmesh;
  std::ifstream in(fname);
  in >> pmesh;

  CGAL::Real_timer timer;
  timer.start();
  auto [shortest_edge_pair_par, longest_edge_pair_par] =
    PMP::minmax_edge_length<CGAL::Parallel_if_available_tag>(pmesh);
  timer.stop();
  std::cout << "minmax edge length (Parallel if available) took: " << timer.time() << std::endl;

  timer.reset();
  timer.start();
  auto [shortest_edge_pair, longest_edge_pair] =
    PMP::minmax_edge_length<CGAL::Sequential_tag>(pmesh);
  timer.stop();
  std::cout << "minmax edge length (Sequential) took: " << timer.time() << std::endl;
}

int main(int argc, char* argv[])
{
  const std::string filename_polyhedron =
    (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/mech-holes-shark.off");

  test<CGAL::Surface_mesh<Epic::Point_3>>(filename_polyhedron);
  test<CGAL::Polyhedron_3<Epic>>(filename_polyhedron);

  std::cerr << "All done." << std::endl;
  return 0;
}
