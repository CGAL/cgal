#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Curvatures/interpolated_corrected_curvature_measures.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <unordered_map>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_Kernel;
typedef CGAL::Polyhedron_3<Epic_Kernel> Polyhedron;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  Polyhedron g1;
  const std::string filename = (argc>1) ?
      argv[1] :
      CGAL::data_file_path("meshes/small_bunny.obj");

  if(!CGAL::IO::read_polygon_mesh(filename, g1))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  std::unordered_map<vertex_descriptor, Epic_Kernel::FT> mean_curvature_map, gaussian_curvature_map;
  std::unordered_map<vertex_descriptor, std::tuple<
      Epic_Kernel::FT,
      Epic_Kernel::FT,
      Epic_Kernel::Vector_3,
      Epic_Kernel::Vector_3
      >> principal_curvature_map;

  PMP::interpolated_corrected_mean_curvature(
      g1,
      boost::make_assoc_property_map(mean_curvature_map)
  );

  // uncomment this to compute a curvature while specifying named parameters
  // Example: an expansion ball radius of 0.5 and a vertex normals map (does not have to depend on positions)

  /*std::unordered_map<vertex_descriptor, Epic_Kernel::Vector_3> vnm;

  PMP::interpolated_corrected_mean_curvature(
      g1,
      boost::make_assoc_property_map(mean_curvature_map),
      CGAL::parameters::ball_radius(0.5).vertex_normal_map(boost::make_assoc_property_map(vnm))
  );*/

  PMP::interpolated_corrected_gaussian_curvature(
      g1,
      boost::make_assoc_property_map(gaussian_curvature_map)
  );
  PMP::interpolated_corrected_principal_curvatures(
      g1,
      boost::make_assoc_property_map(principal_curvature_map)
  );

  int i = 0;
  for (vertex_descriptor v : vertices(g1))
  {
    auto PC = principal_curvature_map[v];
      std::cout << i << ": HC = " << mean_curvature_map[v]
                     << ", GC = " << gaussian_curvature_map[v] << "\n"
                     << ", PC = [ " << std::get<0>(PC) << " , " << std::get<1>(PC) << " ]\n";
    i++;
  }
}
