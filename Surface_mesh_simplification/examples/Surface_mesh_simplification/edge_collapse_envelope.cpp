#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Timer.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Polyhedral_envelope_filter.h>

//bbox
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <iostream>
#include <fstream>

namespace SMS = CGAL::Surface_mesh_simplification;

typedef CGAL::Simple_cartesian<double>                        Kernel;

typedef Kernel::Point_3                                       Point_3;
typedef CGAL::Surface_mesh<Point_3>                           Surface;
typedef CGAL::Timer                                           Timer;

typedef SMS::LindstromTurk_cost<Surface>                      Cost;
typedef SMS::LindstromTurk_placement<Surface>                 Placement;
typedef SMS::Polyhedral_envelope_filter<Kernel,SMS::Bounded_normal_change_filter<> > Filter;


int main(int argc, char** argv)
{
  Surface mesh;
  std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/helmet.off");

  CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, mesh);

  SMS::Edge_count_stop_predicate<Surface> stop(0); // go as far as you can while in the envelope

  Timer t;
  t.start();
  CGAL::Iso_cuboid_3<Kernel> bbox(CGAL::Polygon_mesh_processing::bbox(mesh));

  Point_3 cmin = (bbox.min)();
  Point_3 cmax = (bbox.max)();
  const double diag = CGAL::approximate_sqrt(CGAL::squared_distance(cmin, cmax));

  double eps = (argc>2) ? std::stod(argv[2]) : 0.01*diag;
  std::cout << "eps = " << eps << std::endl;
  Placement placement;
  Filter filter(eps);
  std::cout << "start edge_collapse()" << std::endl;
  SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(Cost()).filter(filter).get_placement(placement));

  std::cout << t.time() << "sec." << std::endl;
  std::ofstream out("out.off");
  out << mesh << std::endl;
  out.close();

  return EXIT_SUCCESS;
}
