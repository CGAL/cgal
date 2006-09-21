#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Surface_mesh_simplification/Polyhedron.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse.h>

// === EXAMPLE SPECIFIC HEADERS BEGINS HERE ===

#include <CGAL/Unique_hash_map.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_pred.h>

// === EXAMPLE SPECIFIC HEADERS ENDS HERE ===


typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface; 

namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  
  // In this example, the flag that indicates if the vertex is fixed is set up as an external
  // map (unitrusively).
  //
  CGAL::Unique_hash_map<Surface::Vertex_const_handle,bool> is_vertex_fixed ;
  
  // This irrealistic loop just illustrates how the flag would be set.
  //
  for ( Surface::Vertex_const_iterator vi = surface.vertices_begin()
      ; vi != surface.vertices_end() 
      ; ++ vi 
      )
    is_vertex_fixed[vi] = false ; // some would be set to true
 
  int r = SMS::edge_collapse(surface
                            ,SMS::Count_ratio_stop_condition<Surface>(0.10)
                            ,CGAL::edge_index_map(surface)
                            .CGAL::vertex_is_fixed_map(boost::make_assoc_property_map(is_vertex_fixed))
                            );
  
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===

  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges." ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
