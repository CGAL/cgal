#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/OpenVDB_mesh_domain_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

#include <openvdb/openvdb.h>
#include <openvdb/io/File.h>
#include <openvdb/io/Stream.h>
#include <openvdb/tools/Interpolation.h>

#include <iostream>
#include <string>
#include <fstream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using GT = K;
using Point = K::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;

using Mesh_domain = CGAL::OpenVDB_mesh_domain_3<K>;
using Tr = typename CGAL::Mesh_triangulation_3<Mesh_domain>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;


int main(int argc, char** argv)
{
  if (argc < 5)
  {
    std::cerr << "Usage: ./main <VDB_FILE> <angular_bound> <radius_bound> <distance_bound>\n";
    return 1;
  }

  const std::string input_filename = argv[1];
  const double angular_bound = std::stod(argv[2]);
  const double radius_bound = std::stod(argv[3]);
  const double distance_bound = std::stod(argv[4]);

  std::cout << "Meshing: " << input_filename << ", angular bound: " << angular_bound << ", radius bound: " << radius_bound << ", distance bound: " << distance_bound << '\n';

  std::ifstream ifs(input_filename, std::ios::binary);

  openvdb::initialize();
  openvdb::io::Stream stream(ifs);
  openvdb::GridPtrVecPtr grids = stream.getGrids();
  openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(grids->front());

  Mesh_domain surface(grid);

  std::cout << "OpenVdb Surface created" << std::endl;

  Mesh_criteria criteria(CGAL::parameters::facet_size = radius_bound,
                         CGAL::parameters::facet_angle = angular_bound,
                         CGAL::parameters::facet_distance = distance_bound);

  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(surface, criteria,
                                      CGAL::parameters::no_perturb(),
                                      CGAL::parameters::no_exude(),
                                      CGAL::parameters::manifold());

  std::cout << "Meshing done" << std::endl;

  Surface_mesh sm;
  CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, sm);

  std::ofstream ofs("mesh.ply");
  CGAL::IO::write_PLY(ofs, sm, "OpenVDB surface");

  return EXIT_SUCCESS;
}
