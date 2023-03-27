#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <unordered_map>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef CGAL::Polyhedron_3<Epic_kernel> Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  Polyhedron polyhedron;
  const std::string filename = (argc > 1) ?
    argv[1] :
    CGAL::data_file_path("meshes/sphere.off");

  if (!CGAL::IO::read_polygon_mesh(filename, polyhedron))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  boost::property_map<Polyhedron, CGAL::dynamic_vertex_property_t<Epic_kernel::FT>>::type
    mean_curvature_map = get(CGAL::dynamic_vertex_property_t<Epic_kernel::FT>(), polyhedron),
    Gaussian_curvature_map = get(CGAL::dynamic_vertex_property_t<Epic_kernel::FT>(), polyhedron);

  boost::property_map<Polyhedron, CGAL::dynamic_vertex_property_t<PMP::Principal_curvatures_and_directions<Epic_kernel>>>::type
    principal_curvatures_and_directions_map =
    get(CGAL::dynamic_vertex_property_t<PMP::Principal_curvatures_and_directions<Epic_kernel>>(), polyhedron);

  PMP::interpolated_corrected_mean_curvature(polyhedron, mean_curvature_map);

  PMP::interpolated_corrected_Gaussian_curvature(polyhedron, Gaussian_curvature_map);

  PMP::interpolated_corrected_principal_curvatures_and_directions(polyhedron, principal_curvatures_and_directions_map);

  // uncomment this to compute a curvature while specifying named parameters
  // Example: an expansion ball radius of 0.5 and a vertex normals map (does not have to depend on positions)

  /*std::unordered_map<vertex_descriptor, Epic_Kernel::Vector_3> vnm;

  PMP::interpolated_corrected_mean_curvature(
      polyhedron,
      mean_curvature_map,
      CGAL::parameters::ball_radius(0.5).vertex_normal_map(boost::make_assoc_property_map(vnm))
  );*/

  // This function can be used to compute multiple curvature types by specifiying them as named parameters
  // This is more efficient than computing each one separately (shared computations).
  PMP::interpolated_corrected_curvatures(
    polyhedron,
    CGAL::parameters::vertex_mean_curvature_map(mean_curvature_map)
    .vertex_principal_curvatures_and_directions_map(principal_curvatures_and_directions_map));

  int i = 0;
  for (vertex_descriptor v : vertices(polyhedron))
  {
    auto PC = get(principal_curvatures_and_directions_map, v);
    std::cout << i << ": HC = " << get(mean_curvature_map, v)
      << ", GC = " << get(Gaussian_curvature_map, v) << "\n"
      << ", PC = [ " << PC.min_curvature << " , " << PC.max_curvature << " ]\n";
    i++;
  }
}
