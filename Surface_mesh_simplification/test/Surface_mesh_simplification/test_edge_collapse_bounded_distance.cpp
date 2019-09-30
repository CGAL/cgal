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


//Timer
#include <CGAL/Timer.h>
//filter{
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_distance_placement.h>
//} end filter

namespace SMS = CGAL::Surface_mesh_simplification ;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface;

int main( int argc, char** argv )
{
  Surface surface_mesh;
  std::ifstream is(argc > 1 ? argv[1] : "input.off");
  is >> surface_mesh;
  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface mesh drops below the specified number (1000)
  SMS::Count_stop_predicate<Surface> stop(num_halfedges(surface_mesh)/10);

  typedef SMS::Bounded_distance_placement<SMS::LindstromTurk_placement<Surface> > Filtered_placement;
  typedef SMS::LindstromTurk_placement<Surface> Placement;
  double placement1_time, placement2_time;

  std::cout<<"input has "<<num_vertices(surface_mesh)<<" vertices."<<std::endl;
  double input_sq_dist=0.000000001;
  std::cout<<"threshold squared distance = "<<input_sq_dist<<std::endl;
  CGAL::Timer timer;
  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface mesh lack an "id()" field.
  Surface copy_mesh = surface_mesh;
  Placement placement1;
  Filtered_placement placement2(input_sq_dist, placement1);
  timer.start();

  SMS::edge_collapse( surface_mesh,
                      stop,
                      CGAL::parameters::get_cost (SMS::LindstromTurk_cost<Surface>())
                      .get_placement(placement2)
                      );
  timer.stop();
  placement2_time=timer.time();
  timer.reset();
  timer.start();
  SMS::edge_collapse( copy_mesh,
                      stop,
                      CGAL::parameters::get_cost (SMS::LindstromTurk_cost<Surface>())
                      .get_placement(placement1)
                      );
  timer.stop();
  placement1_time=timer.time();
  surface_mesh.collect_garbage();
  copy_mesh.collect_garbage();
  std::cout<<"Filtered placement time = "<<placement2_time<<"s."<<std::endl;
  std::cout<<" Not filtered placement took "<<placement1_time<<"s."<<std::endl;
  std::cout<<"There are "<<num_vertices(surface_mesh)<<" vertices left when filtered and "<<num_vertices(copy_mesh)<<" when not filtered."<<std::endl;
  std::ofstream os( "out_filtered.off" );
  os.precision(17);
  os << surface_mesh;
  os.close();
  os.open("out_clear.off");
  os << copy_mesh;
  os.close();
  return EXIT_SUCCESS;
}

// EOF //
