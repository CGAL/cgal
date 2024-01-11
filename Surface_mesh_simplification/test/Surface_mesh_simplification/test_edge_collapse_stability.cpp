#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>

#include <CGAL/IO/polygon_mesh_io.h>

//bbox
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <iostream>
#include <fstream>

namespace SMS = CGAL::Surface_mesh_simplification;

typedef CGAL::Simple_cartesian<double>                        Kernel;

typedef Kernel::Point_3                                       Point_3;
typedef CGAL::Surface_mesh<Point_3>                           Surface;

typedef SMS::LindstromTurk_cost<Surface>                      Cost;
typedef SMS::LindstromTurk_placement<Surface>                 Placement;

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  std::string filename = (argc > 1) ? argv[1] : "data/far_xy.off";

  Surface mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "input has " << num_vertices(mesh) << " vertices." << std::endl;

  CGAL::Iso_cuboid_3<Kernel> bbox(CGAL::Polygon_mesh_processing::bbox(mesh));

  // scale a bit the bounding box bbox, because the kernel is SC
  bbox = { (bbox.min)() - 0.01 * ((bbox.max)() - (bbox.min)()),
           (bbox.max)() + 0.01 * ((bbox.max)() - (bbox.min)()) };

  std::cout << "Bbox: " << bbox << std::endl;

  SMS::Edge_count_stop_predicate<Surface> stop(num_halfedges(mesh)/10);
  Placement placement_ref;

  SMS::edge_collapse(mesh, stop,
                     CGAL::parameters::get_cost(Cost())
                                      .get_placement(placement_ref));

  CGAL::IO::write_polygon_mesh("out.off", mesh, CGAL::parameters::stream_precision(17));

  for(auto v : vertices(mesh))
  {
    if(bbox.has_on_unbounded_side(mesh.point(v)))
    {
      std::cerr << "Error: " << mesh.point(v) << " is outside" << std::endl;
    }
  }

  std::cout << "output has " << vertices(mesh).size() << " vertices." << std::endl;
  return EXIT_SUCCESS;
}
