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
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_pred.h>

// Cost-strategy policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h> 

// === EXAMPLE SPECIFIC HEADERS ENDS HERE ===

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface; 

namespace TSMS = CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  typedef Surface::Halfedge_handle Halfedge_handle ;
  CGAL::Unique_hash_map<Halfedge_handle,void*> edge2ptr ;
  for ( Surface::Halfedge_iterator hi = surface.halfedges_begin()
      ; hi != surface.halfedges_end() 
      ; ++ hi
      )
    edge2ptr[hi] = 0 ;

  // In this example, wich indicates that the cost and placement values 
  // should be computed on demand with no caching.
  // That is specified by passing the "Set_empty_collapse_data" policy
  // along with the GetCost and GetPlacement policies that do the
  // actual on-deman computation.
  //   
  int r = TSMS::edge_collapse(surface
                             ,TSMS::Count_ratio_stop_condition<Surface>(0.10) // StopCondition
                             ,boost::make_assoc_property_map(edge2ptr)        // EdgeExtraPointerMap
                             
                             ,CGAL::Vertex_is_fixed_map_always_false<Surface>() // VertexIsFixedMap [default].
                             
                             ,TSMS::Set_empty_collapse_data         <Surface>() // SetCollapseData
                             ,TSMS::Edge_length_cost                <Surface>() // GetCost
                             ,TSMS::Midpoint_placement              <Surface>() // GetPlacement
                             );
      
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===
  
  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
