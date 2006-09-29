#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

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

namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[2]) ; is >> surface ;

  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  int strategy = atoi(argv[1]);

  if ( strategy == LT )
  {
    // This section shows how to use a full cache with the LindstromTurk strategy:
    //
    // The SetCache policy "LindstromTurk_set_cost_and_placement_cache"
    // is the one actually computing and caching the cost and placement values. 
    //
    // Since that policy is computing the actual values, the other two policies
    // are set to "Cached_cost" and "Cached_placement", since they only need
    // to get the values already computed.
    // 
    SMS::edge_collapse(surface
                      ,SMS::Count_stop_condition<Surface>(1000)             
                       
                      ,CGAL::edge_index_map(surface)
                      
                      .CGAL::set_cache    (SMS::LindstromTurk_set_cost_and_placement_cache<Surface>())
                      .CGAL::get_cost     (SMS::Cached_cost     <Surface>())
                      .CGAL::get_placement(SMS::Cached_placement<Surface>())
                      );
  }
  else
  {
    // This section shows how to use a full cache with a non-default cost strategy.
    //
    //
    // The SetCache policy "Set_cost_and_placemet_cache<>" is the one actually
    // computing and caching the cost and placement values.
    // However, the policy object by itself doesn't know how to do that, 
    // so it delegates the computation to external policies passed
    // to its constructor.
    // 
    // In this example, the non-default cost strategy is to use 
    // the edge length for the cost and the edge midpoint for the placement,
    // so the corresponding policies are instantiated and passed to the 
    // Set_cost_and_placemet_cache constructor.
    // 
    // Still, the GetCost and GetPlacement policies **passed to the algorithm** are
    // Cached_cost and Cached_placement as they only need to extract the cached values.
    //
    typedef SMS::Edge_length_cost  <Surface> Compute_cost ;
    typedef SMS::Midpoint_placement<Surface> Compute_placement ;
    
    Compute_cost      compute_cost ;
    Compute_placement compute_placement ;
    
    SMS::Set_cost_and_placement_cache<Surface,Compute_cost,Compute_placement> 
      set_full_cache(compute_cost,compute_placement);
    
    SMS::edge_collapse(surface
                      ,SMS::Count_stop_condition<Surface>(1000)         
                       
                      ,CGAL::edge_index_map(surface)
                      
                      .CGAL::set_cache    (set_full_cache)
                      .CGAL::get_cost     (SMS::Cached_cost     <Surface>())
                      .CGAL::get_placement(SMS::Cached_placement<Surface>())
                      );
  }  
  
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===

  std::cout << "\nFinished: " << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 3 ? argv[3] : "out.off" ) ;  os << surface ;
  
  return 0 ;      
}

// EOF //
