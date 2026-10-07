#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <input_file>\n";
    return 1;
  }

  Surface_mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(argv[1], mesh) || mesh.is_empty())
  {
    std::cerr << "Cannot read " << argv[1] << "\n";
    return 1;
  }

  if(!CGAL::is_triangle_mesh(mesh))
    PMP::triangulate_faces(mesh);

  PMP::autorefine(mesh,
                  CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag()));

  return 0;
}
