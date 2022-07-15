#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iostream>
#include <unordered_map>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <fstream>

#include <boost/graph/graph_traits.hpp>

#include <CGAL/Polygon_mesh_processing/Curvatures/interpolated_corrected_curvature_measures.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <chrono>
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <CGAL/draw_surface_mesh.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Polyhedron_3<Epic> Polyhedron;
typedef CGAL::Surface_mesh<Epic::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;


int main(int argc, char* argv[])
{
  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("small_bunny.obj");

  Mesh g1;
  if(!CGAL::IO::read_polygon_mesh(filename, g1))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  std::unordered_map<vertex_descriptor, Epic::Vector_3> vnm_vec;
  boost::associative_property_map< std::unordered_map<vertex_descriptor, Epic::Vector_3>> vnm(vnm_vec);

  PMP::compute_vertex_normals(g1, vnm);


  std::vector<Epic::FT> mu0_map, mu1_map, mu2_map;

  mu0_map = PMP::interpolated_corrected_measure_i(g1, PMP::MU0_AREA_MEASURE, CGAL::parameters::vertex_normal_map(vnm));
  mu1_map = PMP::interpolated_corrected_measure_i(g1, PMP::MU1_MEAN_CURVATURE_MEASURE, CGAL::parameters::vertex_normal_map(vnm));
  mu2_map = PMP::interpolated_corrected_measure_i(g1, PMP::MU2_GAUSSIAN_CURVATURE_MEASURE, CGAL::parameters::vertex_normal_map(vnm));

  int n = g1.faces().size();

  for (int i = 0; i < n; i++)
  {
      std::cout << mu0_map[i] << "\n";
  }

  std::cout << "\n";

  for (int i = 0; i < n; i++)
  {
      std::cout << mu1_map[i] << "\n";
  }

  std::cout << "\n";

  for (int i = 0; i < n; i++)
  {
      std::cout << mu2_map[i] << "\n";
  }


  CGAL::draw(g1);

  return EXIT_SUCCESS;
}
