#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::Polygon_mesh_processing::parameters;

int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/blobby.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/eight.off";
  std::ifstream input(filename1);

  Mesh mesh1, mesh2;
  if (!input || !(input >> mesh1))
  {
    std::cerr << "First mesh is not a valid off file." << std::endl;
    return 1;
  }
  input.close();
  input.open(filename2);
  if (!input || !(input >> mesh2))
  {
    std::cerr << "Second mesh is not a valid off file." << std::endl;
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
      params::all_default(), // mesh1 named parameters
      params::all_default(), // mesh2 named parameters
      std::make_tuple(
        params::vertex_point_map(get(boost::vertex_point, out_union)), // named parameters for out_union
        params::vertex_point_map(get(boost::vertex_point, out_intersection)), // named parameters for out_intersection
        params::all_default(), // named parameters for mesh1-mesh2 not used
        params::all_default() )// named parameters for mesh2-mesh1 not used)
    );

  if (res[PMP::Corefinement::UNION])
  {
    std::cout << "Union was successfully computed\n";
    std::ofstream output("union.off");
    output.precision(17);
    output << out_union;
  }
  else
    std::cout << "Union could not be computed\n";

  if (res[PMP::Corefinement::INTERSECTION])
  {
    std::cout << "Intersection was successfully computed\n";
    std::ofstream output("intersection.off");
    output.precision(17);
    output << out_intersection;
  }
  else
    std::cout << "Intersection could not be computed\n";

  return 0;
}
