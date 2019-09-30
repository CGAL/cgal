#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>

//AABB_tree
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

//bbox
#include <CGAL/Polygon_mesh_processing/bbox.h>


//Timer
//#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_distance_placement.h>

namespace SMS = CGAL::Surface_mesh_simplification ;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface;


typedef CGAL::AABB_face_graph_triangle_primitive<Surface> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;


int main( int argc, char** argv )
{
  Surface ref_mesh;
  std::ifstream is(argc > 1 ? argv[1] : "data/helmet.off");
  is >> ref_mesh;

  SMS::Count_stop_predicate<Surface> stop(num_halfedges(ref_mesh)/10);

  typedef SMS::Bounded_distance_placement<SMS::LindstromTurk_placement<Surface>, Tree> Filtered_placement;
  typedef SMS::LindstromTurk_placement<Surface> Placement;
  //double placement1_time, placement2_time, tree_time;

  //std::cout<<"input has "<<num_vertices(ref_mesh)<<" vertices."<<std::endl;
  CGAL::Iso_cuboid_3<Kernel> bbox(CGAL::Polygon_mesh_processing::bbox(ref_mesh));


  Kernel::Point_3 cmin = (bbox.min)();
  Kernel::Point_3 cmax = (bbox.max)();
  double diag = std::sqrt(CGAL::squared_distance(cmin,cmax));

  Surface small_mesh = ref_mesh;
  Surface big_mesh = ref_mesh;
  Tree tree( faces(ref_mesh).first, faces(ref_mesh).second, ref_mesh);

  Placement placement_ref;
  Filtered_placement placement_small(0.00005*diag, tree, placement_ref);
  Filtered_placement placement_big(50*diag, tree, placement_ref);

  SMS::edge_collapse( small_mesh,
                      stop,
                      CGAL::parameters::get_cost (SMS::LindstromTurk_cost<Surface>())
                      .get_placement(placement_small)
                      );

  SMS::edge_collapse( big_mesh,
                      stop,
                      CGAL::parameters::get_cost (SMS::LindstromTurk_cost<Surface>())
                      .get_placement(placement_big)
                      );

  SMS::edge_collapse( ref_mesh,
                      stop,
                      CGAL::parameters::get_cost (SMS::LindstromTurk_cost<Surface>())
                      .get_placement(placement_ref)
                      );

  ref_mesh.collect_garbage();
  small_mesh.collect_garbage();
  big_mesh.collect_garbage();
  //std::cout<<"Filtered placement time = "<<placement2_time<<"s."<<std::endl;
  //std::cout<<" Not filtered placement took "<<placement1_time<<"s."<<std::endl;
  //std::cout<<"There are "<<num_vertices(surface_mesh)<<" vertices left when filtered and "<<num_vertices(copy_mesh)<<" when not filtered."<<std::endl;
  assert(num_vertices(big_mesh) == num_vertices(ref_mesh));
  assert(num_vertices(small_mesh) != num_vertices(ref_mesh));

  return EXIT_SUCCESS;
}

// EOF //
