#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Timer.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/FastEnvelope_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>


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


typedef SMS::FastEnvelope_placement<Placement, Kernel>    Filtered_placement;

int main(int argc, char** argv)
{
  Surface ref_mesh;
  std::ifstream is(argc > 1 ? argv[1] : "data/helmet.off");
  is >> ref_mesh;

  SMS::Count_stop_predicate<Surface> stop(num_halfedges(ref_mesh)/10);

  std::cout << "Input has " << num_vertices(ref_mesh) << " vertices and " << num_edges(ref_mesh) << " edges" << std::endl;
  CGAL::Iso_cuboid_3<Kernel> bbox(CGAL::Polygon_mesh_processing::bbox(ref_mesh));

  Point_3 cmin = (bbox.min)();
  Point_3 cmax = (bbox.max)();
  const double diag = CGAL::approximate_sqrt(CGAL::squared_distance(cmin, cmax));

  Surface mesh_cpy = ref_mesh; // need a copy to keep the AABB tree valid
  Surface small_mesh = ref_mesh;
  Surface big_mesh = ref_mesh;
  Surface huge_mesh = ref_mesh;

  CGAL::Timer t;
  t.start();
  Placement placement_ref;
  SMS::edge_collapse(ref_mesh, stop, CGAL::parameters::get_cost(Cost()).get_placement(placement_ref));
  std::cout << "Output has " << vertices(ref_mesh).size() << " vertices and " << edges(ref_mesh).size() << " edges" << std::endl;
  std::cout << t.time() << "sec\n";
  t.reset();

  std::cout << "eps = " << 0.00005*diag << std::endl;
  Filtered_placement placement_small(0.00005*diag, placement_ref);
  SMS::edge_collapse(small_mesh, stop, CGAL::parameters::get_cost(Cost()).get_placement(placement_small));
  std::cout << "Output has " << vertices(small_mesh).size() << " vertices and " << edges(small_mesh).size() << " edges" << std::endl;
  std::cout << t.time() << "sec\n";
  t.reset();

  std::cout << "eps = " << 0.0001*diag << std::endl;
  Filtered_placement placement_big(0.0001, placement_ref);
  SMS::edge_collapse(big_mesh, stop, CGAL::parameters::get_cost(Cost()).get_placement(placement_big));
  std::cout << "Output has " << vertices(big_mesh).size() << " vertices and " << edges(big_mesh).size() << " edges" << std::endl;
  std::cout << t.time() << "sec\n";

 std::cout << "eps = " << 0.0002*diag << std::endl;
  Filtered_placement placement_huge(0.0002, placement_ref);
  SMS::edge_collapse(huge_mesh, stop, CGAL::parameters::get_cost(Cost()).get_placement(placement_huge));
  std::cout << "Output has " << vertices(huge_mesh).size() << " vertices and " << edges(huge_mesh).size() << " edges" << std::endl;
  std::cout << t.time() << "sec\n";

  std::ofstream out("big.off");
  out << big_mesh << std::endl;
  out.close();

  std::cout << "no filtering: " << vertices(ref_mesh).size() << " vertices left" << std::endl;
  std::cout << "huge filtering distance: " << vertices(huge_mesh).size() << " vertices left" << std::endl;
  std::cout << "large filtering distance: " << vertices(big_mesh).size() << " vertices left" << std::endl;
  std::cout << "small filtering distance: " << vertices(small_mesh).size() << " vertices left" << std::endl;


  assert(vertices(ref_mesh).size() < vertices(small_mesh).size());
  assert(vertices(huge_mesh).size() < vertices(small_mesh).size());
  assert(vertices(ref_mesh).size() < vertices(big_mesh).size());

  return EXIT_SUCCESS;
}
