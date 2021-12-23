#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

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

  std::cout << "Number of vertices before corefinement "
            << num_vertices(mesh1) << " and "
            << num_vertices(mesh2) << "\n";

  PMP::corefine(mesh1,mesh2);

  std::cout << "Number of vertices after corefinement "
            << num_vertices(mesh1) << " and "
            << num_vertices(mesh2) << "\n";

  CGAL::IO::write_polygon_mesh("mesh1_refined.off", mesh1, CGAL::parameters::stream_precision(17));
  CGAL::IO::write_polygon_mesh("mesh2_refined.off", mesh2, CGAL::parameters::stream_precision(17));

  return 0;
}
