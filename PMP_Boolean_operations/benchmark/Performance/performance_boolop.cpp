#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <input_file>\n";
    return 1;
  }

  Mesh tm1;
  if(!CGAL::IO::read_polygon_mesh(argv[1], tm1) || tm1.is_empty())
  {
    std::cerr << "Cannot read " << argv[1] << "\n";
    return 1;
  }

  Mesh tm2 = tm1;
  const CGAL::Bbox_3 bb = PMP::bbox(tm1);
  const Vector_3 translation((bb.xmax()-bb.xmin()) * 0.2,
                             (bb.ymax()-bb.ymin()) * 0.2,
                             (bb.zmax()-bb.zmin()) * 0.2);
  for(auto v : tm2.vertices())
    tm2.point(v) = tm2.point(v) + translation;

  Mesh out;
  PMP::corefine_and_compute_union(
    tm1, tm2, out,
    CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag()));
  return 0;
}
