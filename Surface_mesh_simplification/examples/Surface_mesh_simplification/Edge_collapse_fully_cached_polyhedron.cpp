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

// Cost-strategy policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>  
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

// === EXAMPLE SPECIFIC HEADERS ENDS HERE ===


typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface; 

enum { LT, MP } ;

namespace TSMS = CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  int strategy = atoi(argv[1]);
  
  std::ifstream is(argv[2]) ; is >> surface ;

  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  CGAL::Unique_hash_map<Surface::Halfedge_handle,void*> edge2ptr ;
  for ( Surface::Halfedge_iterator hi = surface.halfedges_begin(); hi != surface.halfedges_end() ; ++ hi )
    edge2ptr[hi] = 0 ;
    
  if ( strategy == LT )
  {
    // This section shows how to use a full cache with the LindstromTurk strategy:
    //
    // The policy "Set_full_collapse_data_LindstromTurk" is the one actually 
    // computing and caching the cost and placement values. 
    //
    // Since that policy is computing the actual values, the other two policies
    // are set to "Cached_cost" and "Cached_placement", since they only need
    // to get the values already computed.
    // 
    TSMS::edge_collapse(surface
                       ,TSMS::Count_stop_condition<Surface>(1000)             // StopCondition
                       ,boost::make_assoc_property_map(edge2ptr)              // EdgeExtraPointerMap 
                       
                       ,CGAL::Vertex_is_fixed_map_always_false<Surface>()     // VertexIsFixedMap
                       
                       ,TSMS::Set_full_collapse_data_LindstromTurk<Surface>() // SetCollapseData   
                       ,TSMS::Cached_cost<Surface>()                          // GetCost
                       ,TSMS::Cached_placement<Surface>()                     // GetPlacement
                       );
  }
  else
  {
    // This section shows how to use a full cache with a non-default cost strategy.
    //
    //
    // The policy "Set_full_collapse_data<>" is the one actually computing
    // and caching the cost and placement value.
    // However, the policy object by itself it doesn't know how to do hat, 
    // so it delegates the computation to external policies passed
    // to the SetCollapseData policy object constructor.
    // 
    // In this example, the non-default cost strategy is to use 
    // the edge length for the cost and the edge midpoint for the placement, so
    // we instantiate policy objects for that and pass them to the 
    // Set_full_collapse_data<> policy object.
    // 
    // Still, the GetCost and GetPlacement policies **passed to the algorithm** are
    // Cached_cost and Cached_placement as they only need to extract the cached values 
    // from the collapse-data.
    //
    typedef TSMS::Edge_length_cost  <Surface> Compute_cost ;
    typedef TSMS::Midpoint_placement<Surface> Compute_placement ;
    
    Compute_cost      compute_cost ;
    Compute_placement compute_placement ;
    
    TSMS::Set_full_collapse_data<Surface,Compute_cost,Compute_placement> 
      set_collapse_data(compute_cost,compute_placement);
    
    TSMS::edge_collapse(surface
                       ,TSMS::Count_stop_condition<Surface>(1000)         // StopCondition
                       ,boost::make_assoc_property_map(edge2ptr)          // EdgeExtraPointerMap 
                       
                       ,CGAL::Vertex_is_fixed_map_always_false<Surface>() // VertexIsFixedMap [default]
                       
                       ,set_collapse_data                                 // SetCollapseData   
                       ,TSMS::Cached_cost<Surface>()                      // GetCost
                       ,TSMS::Cached_placement<Surface>()                 // GetPlacement
                       );
  }  
  
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===

  std::cout << "\nFinished: " << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 3 ? argv[3] : "out.off" ) ;  os << surface ;
  
  return 0 ;      
}

// EOF //
