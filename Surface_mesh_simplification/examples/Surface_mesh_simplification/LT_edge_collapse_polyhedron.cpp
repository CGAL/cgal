#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Surface_mesh_simplification/polyhedron.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// === EXAMPLE SPECIFIC HEADERS BEGINS HERE ===

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

// === EXAMPLE SPECIFIC HEADERS ENDS HERE ===

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface; 

namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  // This is a stop-condition policy (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface drops below the specified number (1000)
  SMS::Count_stop_predicate<Surface> stop_policy(1000);
     
  // This the actual call to the simplification algorithm.
  // The surface and stop conditions are mandatory arguments.
  // The third argument is needed because the edges on this
  // surface lack an "id()" field.
  int r = SMS::edge_collapse(surface
                            ,stop_policy                 
                            ,CGAL::edge_index_map(surface)
                            );
  
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===

  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
