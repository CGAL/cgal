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

// === EXAMPLE SPECIFIC HEADERS ENDS HERE ===


typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface; 

// === EXAMPLE SPECIFIC DETAILS BEGINS HERE ===

typedef Surface::Halfedge_handle Halfedge_handle ;
typedef Surface::Vertex_handle   Vertex_handle ;

// The following is the Visitor that keeps track of the simplification process.
// In this example the progress is printed real-time and a few statistics are
// recorded (and printed in the end).
//
struct Visitor
{
  Visitor() 
    : collected(0)
    , processed(0)
    , collapsed(0)
    , non_collapsable(0)
    , cost_uncomputable(0) 
  {} 

  // Called on algorithm entry  
  void OnStarted( Surface& aSurface ) {} 
  
  // Called on algorithm exit  
  void OnFinished ( Surface& aSurface ) { std::cerr << "\n" << std::flush ; } 
  
  // Called when the stop condition returned true
  void OnStopConditionReached( Surface& aSurface ) {} 
  
  // Called during the collecting phase for each edge collected.
  void OnCollected( Halfedge_handle const& aEdge
                  , bool                   aIsFixed
                  , Surface&               aSurface
                  )
  {
    ++ collected ;
    std::cerr << "\rEdges collected: " << collected << std::flush ;
 }                
  
  // Called during the processing phase for each edge processed.
  // If aVertex is a valid handle then the edge was collapsed and aVertex is the replacement.
  void OnProcessed(Halfedge_handle const&  aEdge
                  ,Surface&                aSurface
                  ,boost::optional<double> aCost
                  ,Vertex_handle const&    aVertex
                  )
  {
    ++ processed ;
    if ( aVertex == Vertex_handle() )
    {
      if ( !aCost )
           ++ cost_uncomputable ;
      else ++ non_collapsable ;
    }
    else 
    {
      ++ collapsed;
    }
  }                
  
  // Call at each step in the processing phase (when each edge is selected for processing) before the
  // stop condition is evaluated.
  void OnStep( Halfedge_handle const& aEdge
             , Surface&               aSurface
             , std::size_t            aInitial
             , std::size_t            aCurrent
             )
  {
    if ( aCurrent == aInitial )
      std::cerr << "\n" << std::flush ;
    std::cerr << "\r" << aCurrent << std::flush ;
  }                
  
  std::size_t  collected
             , processed
             , collapsed
             , non_collapsable
             , cost_uncomputable ; 
} ;

// === EXAMPLE SPECIFIC DETAILS ENDS HERE ===

namespace TSMS = CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  CGAL::Unique_hash_map<Surface::Halfedge_handle,void*> edge2ptr ;
  for ( Surface::Halfedge_iterator hi = surface.halfedges_begin()
      ; hi != surface.halfedges_end() 
      ; ++ hi 
      )
    edge2ptr[hi] = 0 ;
 
  Visitor visitor ;

  // Since the visitor is the last parameter, all the parameters, which could have been
  // ommited, must be explicitely passed (so this examples shows all the defaults).
  //  
  TSMS::LindstromTurk_params params ;
    
  int r = TSMS::edge_collapse(surface
                             ,TSMS::Count_ratio_stop_condition<Surface>(0.10)          // StopCondition
                             ,boost::make_assoc_property_map(edge2ptr)                 // EdgeExtraPointerMap
                             
                             // These are the same as the default arguments
                             ,CGAL::Vertex_is_fixed_map_always_false<Surface>()        // VertexIsFixedMap
                             ,TSMS::Set_partial_collapse_data_LindstromTurk<Surface>() // SetCollapseData
                             ,TSMS::Cached_cost<Surface>()                             // GetCost
                             ,TSMS::LindstromTurk_placement<Surface>()                 // GetPlacement
                             ,&params                                                  // GetCost parameters
                             ,&params                                                  // GetPlacement parameters
                             
                             ,&visitor
                             );

  
  std::cout << "\nEdges collected: " << visitor.collected
            << "\nEdges proccessed: " << visitor.processed
            << "\nEdges collapsed: " << visitor.collapsed
            << std::endl
            << "\nEdges not collapsed due to topological constrians: " 
            << visitor.non_collapsable
            << "\nEdge not collapsed due to computational constrians: " 
            << visitor.cost_uncomputable 
            << std::endl ; 
            
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===
  
  std::cout << "\nFinished...\n" << r << " edges removed.\n"  
            << (surface.size_of_halfedges()/2) << " final edges." ;
  
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
