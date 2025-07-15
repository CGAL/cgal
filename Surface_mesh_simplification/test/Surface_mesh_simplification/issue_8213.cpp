#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

CGAL::Bbox_3 bbox_g;
double max_bbox_g;
#define CGAL_CHECK_EXPENSIVE

#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE 5
#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_LT_TRACE 4

void Surface_simplification_external_trace(const std::string& s)
{
  std::cout << s << std::endl;
}

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
// #include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounding_box_filter.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Face_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

using Kernel = CGAL::Simple_cartesian<double>;
// using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Surface_mesh = CGAL::Surface_mesh<Kernel::Point_3>;

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const char* filename = (argc > 1) ? argv[1] : "data/issue_8213.off";

  Surface_mesh sm;

  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Error: failed to read input data" << std::endl;
    return EXIT_FAILURE;
  }

  auto bb=PMP::bbox(sm);
  std::cout << "Bbox:" << bb << std::endl;
  std::cout << "Input mesh has " << num_vertices(sm) << " vertices" << std::endl;
  std::cout << "Input mesh has " << num_faces(sm) << " faces" << std::endl;


  SMS::Face_count_stop_predicate<Surface_mesh> stop(1);
   SMS::edge_collapse(
    sm,
    stop,
    CGAL::parameters::
      //  filter(SMS::Bounded_normal_change_filter<>())  //<SMS::Bounding_box_filter<>>(SMS::Bounding_box_filter<>()))
      get_cost(SMS::LindstromTurk_cost<Surface_mesh>())
      .get_placement(SMS::LindstromTurk_placement<Surface_mesh>())
    );
  CGAL::IO::write_OFF(std::cout, sm, CGAL::parameters::stream_precision(17));

  for(auto v : vertices(sm))
  {
    // To be within the bounding box isn't a guarantee, but here it is a sufficient test
    // to check if things went awry
    if(!CGAL::do_overlap(bb, sm.point(v).bbox()))
    {
      std::cerr << "Error: " << sm.point(v) << " is outside the initial bbox" << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
