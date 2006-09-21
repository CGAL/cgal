#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Unique_hash_map.h>

#include <CGAL/Surface_mesh_simplification/Polyhedron.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse.h>

// === EXAMPLE SPECIFIC HEADERS BEGINS HERE ===

#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_pred.h>

// === EXAMPLE SPECIFIC HEADERS ENDS HERE ===

typedef CGAL::Simple_cartesian<double> Kernel;

// === EXAMPLE SPECIFIC DETAILS BEGINS HERE ===

//
// Setup an enriched polyhedron type which stores a ID in the halfedges
//
typedef CGAL::Polyhedron_3<Kernel,Polyhedron_items_with_id_3> Surface; 

// === EXAMPLE SPECIFIC DETAILS ENDS HERE ===

namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  
  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  // The halfedges in this polyhedron have an "id()" field 
  // which the default "edge_index_map()" uses to get the index
  // of an edge.
  // However, the Polyhedron_3 class doesn't assign any value to
  // this id(), so we must do it here:
  int index = 0 ;
  
  for( typename Surface::Halfedge_iterator eb = surface.halfedges_begin()
     , ee = surface.halfedges_end()
     ; eb != ee
     ; ++ eb
     ) 
    eb->id() = index++;

  SMS::Count_stop_condition<Surface> stop_policy(1000);
     
  int r = SMS::edge_collapse(surface,stop_policy);
  
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===

  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
