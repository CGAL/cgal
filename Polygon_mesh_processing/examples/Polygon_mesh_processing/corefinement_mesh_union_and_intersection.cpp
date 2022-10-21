#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <array>
#include <iostream>
#include <string>
#include <tuple>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>                      Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

int main(int argc, char* argv[])
{
  const std::string filename1 = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  const std::string filename2 = (argc > 2) ? argv[2] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh1, mesh2;
  if(!PMP::IO::read_polygon_mesh(filename1, mesh1) || !PMP::IO::read_polygon_mesh(filename2, mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  Mesh out_union, out_intersection;
  std::array<boost::optional<Mesh*>, 4> output;
  output[PMP::Corefinement::UNION] = &out_union;
  output[PMP::Corefinement::INTERSECTION] = &out_intersection;

  // for the example, we explicit the named parameters, this is identical to
  // PMP::corefine_and_compute_boolean_operations(mesh1, mesh2, output)
  std::array<bool, 4> res =
    PMP::corefine_and_compute_boolean_operations(
      mesh1, mesh2,
      output,
      params::default_values(), // mesh1 named parameters
      params::default_values(), // mesh2 named parameters
      std::make_tuple(
        params::vertex_point_map(get(boost::vertex_point, out_union)), // named parameters for out_union
        params::vertex_point_map(get(boost::vertex_point, out_intersection)), // named parameters for out_intersection
        params::default_values(), // named parameters for mesh1-mesh2 not used
        params::default_values() )// named parameters for mesh2-mesh1 not used)
    );

  if (res[PMP::Corefinement::UNION])
  {
    std::cout << "Union was successfully computed\n";
    CGAL::IO::write_polygon_mesh("union.off", out_union, CGAL::parameters::stream_precision(17));
  }
  else
    std::cout << "Union could not be computed\n";

  if (res[PMP::Corefinement::INTERSECTION])
  {
    std::cout << "Intersection was successfully computed\n";
    CGAL::IO::write_polygon_mesh("intersection.off", out_intersection, CGAL::parameters::stream_precision(17));
  }
  else
    std::cout << "Intersection could not be computed\n";

  return 0;
}
