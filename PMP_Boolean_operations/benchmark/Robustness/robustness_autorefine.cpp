#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

enum EXIT_CODES { VALID_OUTPUT=0,
                  INVALID_INPUT=1,
                  SELF_INTERSECTING_OUTPUT=3,
                  SIGSEGV=10,
                  SIGSABRT=11,
                  SIGFPE=12,
                  TIMEOUT=13 };

int main(int argc, char** argv)
{
  if(argc < 2)
    return INVALID_INPUT;

  Surface_mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(argv[1], mesh) || mesh.is_empty())
    return INVALID_INPUT;

  if(!CGAL::is_triangle_mesh(mesh))
    PMP::triangulate_faces(mesh);

  PMP::autorefine(mesh, CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag()));

  if(PMP::does_self_intersect<CGAL::Parallel_if_available_tag>(mesh))
    return SELF_INTERSECTING_OUTPUT;

  return VALID_OUTPUT;
}
