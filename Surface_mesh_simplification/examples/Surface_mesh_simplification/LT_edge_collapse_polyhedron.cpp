#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Unique_hash_map.h>

#include <CGAL/Surface_mesh_simplification/Polyhedron.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse.h>

// === EXAMPLE SPECIFIC HEADERS BEGINS HERE ===

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_pred.h>

// === EXAMPLE SPECIFIC HEADERS ENDS HERE ===

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface; 

namespace TSMS = CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  // The following code sets up the per-edge extra pointer needed by the algorithm 
  // as an external map becasue the Polyhedron used in this example has not been 
  // enriched to support the extra pointer directly (as shown in the next example)
  CGAL::Unique_hash_map<Surface::Halfedge_handle,void*> edge2ptr ;
  for ( Surface::Halfedge_iterator hi = surface.halfedges_begin()
      ; hi != surface.halfedges_end() 
      ; ++ hi 
      )
    edge2ptr[hi] = 0 ;

  // This is a stop-condition policy (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface drops below the specified number (1000)
  TSMS::Count_stop_condition<Surface> stop_policy(1000);
     
  // This the actual call to the simplification algorithm.
  // The surface and stop conditions are mandatory arguments,
  // while third argument could have been omited, as shown in the next example.
  int r = TSMS::edge_collapse(surface
                             ,stop_policy                              // StopCondition 
                             ,boost::make_assoc_property_map(edge2ptr) // EdgeExtraPointerMap 
                             );
  
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===

  std::cout << "\nFinished...\n" << r << " edges removed.\n"  << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
