#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_distance_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>

//AABB_tree
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

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

typedef CGAL::AABB_face_graph_triangle_primitive<Surface>     Primitive;
typedef CGAL::AABB_traits_3<Kernel, Primitive>                  Traits;
typedef CGAL::AABB_tree<Traits>                               Tree;

typedef SMS::Bounded_distance_placement<Placement, Kernel>    Filtered_placement;
typedef SMS::Bounded_distance_placement<Placement, Tree>      Filtered_placement_with_tree;

int main(int argc, char** argv)
{
  Surface ref_mesh;
  std::ifstream is(argc > 1 ? argv[1] : "data/helmet.off");
  is >> ref_mesh;

  SMS::Edge_count_stop_predicate<Surface> stop(num_halfedges(ref_mesh)/10);

  std::cout << "input has " << num_vertices(ref_mesh) << " vertices." << std::endl;
  CGAL::Iso_cuboid_3<Kernel> bbox(CGAL::Polygon_mesh_processing::bbox(ref_mesh));

  Point_3 cmin = (bbox.min)();
  Point_3 cmax = (bbox.max)();
  const double diag = CGAL::approximate_sqrt(CGAL::squared_distance(cmin, cmax));

  Surface mesh_cpy = ref_mesh; // need a copy to keep the AABB tree valid
  Surface small_mesh = ref_mesh;
  Surface big_mesh = ref_mesh;
  Tree tree(faces(mesh_cpy).first, faces(mesh_cpy).second, mesh_cpy);

  Placement placement_ref;
  SMS::edge_collapse(ref_mesh, stop, CGAL::parameters::get_cost(Cost()).get_placement(placement_ref));

  Filtered_placement_with_tree placement_small(0.00005*diag, tree, placement_ref);
  SMS::edge_collapse(small_mesh, stop, CGAL::parameters::get_cost(Cost()).get_placement(placement_small));

  Filtered_placement placement_big(0.005*diag, placement_ref); // lazily builds the AABB tree
  SMS::edge_collapse(big_mesh, stop, CGAL::parameters::get_cost(Cost()).get_placement(placement_big));

  std::cout << "no filtering: " << vertices(ref_mesh).size() << " vertices left" << std::endl;
  std::cout << "large filtering distance: " << vertices(big_mesh).size() << " vertices left" << std::endl;
  std::cout << "small filtering distance: " << vertices(small_mesh).size() << " vertices left" << std::endl;

  assert(vertices(ref_mesh).size() < vertices(small_mesh).size());
  assert(vertices(big_mesh).size() < vertices(small_mesh).size());
  assert(vertices(ref_mesh).size() < vertices(big_mesh).size());

  return EXIT_SUCCESS;
}
