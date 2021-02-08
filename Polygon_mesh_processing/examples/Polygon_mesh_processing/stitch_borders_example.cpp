
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Polyhedron_3<K>                                 Polyhedron;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/quads_to_stitch.off");

  Polyhedron mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::cout << "Before stitching : " << std::endl;
  std::cout << "\t Number of vertices  :\t" << mesh.size_of_vertices() << std::endl;
  std::cout << "\t Number of halfedges :\t" << mesh.size_of_halfedges() << std::endl;
  std::cout << "\t Number of facets    :\t" << mesh.size_of_facets() << std::endl;

  PMP::stitch_borders(mesh);

  std::cout << "Stitching done : " << std::endl;
  std::cout << "\t Number of vertices  :\t" << mesh.size_of_vertices() << std::endl;
  std::cout << "\t Number of halfedges :\t" << mesh.size_of_halfedges() << std::endl;
  std::cout << "\t Number of facets    :\t" << mesh.size_of_facets() << std::endl;

  CGAL::IO::write_polygon_mesh("mesh_stitched.off", mesh, CGAL::parameters::stream_precision(17));

  return 0;
}
