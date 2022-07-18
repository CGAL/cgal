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
typedef std::unordered_map<face_descriptor, Epic::FT> FaceMeasureMap_tag;


int main(int argc, char* argv[])
{
  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/small_bunny.obj");

  Mesh g1;
  if(!CGAL::IO::read_polygon_mesh(filename, g1))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  std::unordered_map<vertex_descriptor, Epic::Vector_3> vnm_init;
  boost::associative_property_map< std::unordered_map<vertex_descriptor, Epic::Vector_3>> vnm(vnm_init);

  PMP::compute_vertex_normals(g1, vnm);

  FaceMeasureMap_tag mu0_init, mu1_init, mu2_init;
  boost::associative_property_map<FaceMeasureMap_tag> mu0_map(mu0_init), mu1_map(mu1_init), mu2_map(mu2_init);

  PMP::interpolated_corrected_measure_mesh(g1, mu0_map, PMP::MU0_AREA_MEASURE);
  PMP::interpolated_corrected_measure_mesh(g1, mu1_map, PMP::MU1_MEAN_CURVATURE_MEASURE, CGAL::parameters::vertex_normal_map(vnm));
  PMP::interpolated_corrected_measure_mesh(g1, mu2_map, PMP::MU2_GAUSSIAN_CURVATURE_MEASURE, CGAL::parameters::vertex_normal_map(vnm));

  for (face_descriptor f: g1.faces())
  {
      std::cout << f.idx() << ": " << get(mu0_map, f) << ", " << get(mu1_map, f) << ", " << get(mu2_map, f) << ", " << "\n";
  }
}
