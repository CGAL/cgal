#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface_mesh; 

namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv ) 
{
  Surface_mesh surface_mesh;

  std::ifstream is(argc > 1 ? argv[1] : "data/fold.off") ; is >> surface_mesh ;

  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface mesh drops below the specified number (1000)
  SMS::Count_stop_predicate<Surface_mesh> stop(num_halfedges(surface_mesh)/2 - 1);
     
typedef SMS::Bounded_normal_change_placement<SMS::LindstromTurk_placement<Surface_mesh> > Placement;


  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface mesh lack an "id()" field.
  SMS::edge_collapse( surface_mesh,
                      stop,
                      CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,surface_mesh))
                        .halfedge_index_map  (get(CGAL::halfedge_external_index  ,surface_mesh))
                        .get_cost (SMS::LindstromTurk_cost<Surface_mesh>())
                        .get_placement(Placement())
                      );

  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface_mesh ;

  return 0 ;      
}
