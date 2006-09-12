// Lindstrom-Turk edge-collapse with ordinary surface and external per-edge pointer.
//
// Explicit arguments:
// 
//   The surface is an ordinary Polyhedron_3
//
//   The stop condition is to finish when the number of undirected edges 
//   drops below a certain absolute number.
//
//   An external map is used to store the per-edge extra pointer unintrusively.
//
//   No vertex is fixed.
//
//   The cost strategy is Lindstrom-Turk with partial cache (only cost cached)+
//
//   A visitor is used to track the simplification.
//
#include <iostream>
#include <iomanip>
#include <fstream>

#include <CGAL/Simple_cartesian.h>

// Target surface type. (this include Polyhedron_3.h itself)
#include <CGAL/Surface_mesh_simplification/Polyhedron.h>

// Policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_pred.h>

// Simplification method
#include <CGAL/Surface_mesh_simplification/Edge_collapse.h>

#include <CGAL/IO/Polyhedron_iostream.h>

using namespace std ;
using namespace boost ;
using namespace CGAL ;

typedef Simple_cartesian<double> Kernel;
typedef Kernel::Vector_3         Vector;
typedef Kernel::Point_3          Point;

typedef Polyhedron_3<Kernel> Surface; 

typedef Surface::Halfedge_handle Halfedge_handle ;
typedef Surface::Vertex_handle   Vertex_handle ;

using namespace CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;

struct Visitor
{
  Visitor() : collected(0), processed(0), collapsed(0), non_collapsable(0), cost_uncomputable(0) {} 

  // Called on algorithm entry  
  void OnStarted( Surface& aSurface ) {} 
  
  // Called on algorithm exit  
  void OnFinished ( Surface& aSurface ) { cerr << "\n" << flush ; } 
  
  // Called when the stop condition returned true
  void OnStopConditionReached( Surface& aSurface ) {} 
  
  // Called during the collecting phase for each edge collected.
  void OnCollected( Halfedge_handle const& aEdge
                  , bool                   aIsFixed
                  , Surface&               aSurface
                  )
  {
    ++ collected ;
    cerr << "\rEdges collected: " << collected << flush ;
 }                
  
  // Called during the processing phase for each edge processed.
  // If aVertex is a valid handle then the edge was collapsed and aVertex is the replacement.
  void OnProcessed(Halfedge_handle const& aEdge
                  ,Surface&               aSurface
                  ,optional<double>       aCost
                  ,Vertex_handle const&   aVertex
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
  void OnStep(Halfedge_handle const& aEdge, Surface& aSurface, size_t aInitial, size_t aCurrent)
  {
    if ( aCurrent == aInitial )
      cerr << "\n" << flush ;
    cerr << "\r" << aCurrent << flush ;
  }                
  
  size_t collected, processed, collapsed, non_collapsable, cost_uncomputable ; 
} ;


int main( int argc, char** argv ) 
{
  Surface surface; 
  
  ifstream is(argv[1]) ;
  is >> surface ;

  // Extra pointer external map
  Unique_hash_map<Halfedge_handle,void*> edge2ptr ;
  for ( Surface::Halfedge_iterator hi = surface.halfedges_begin(); hi != surface.halfedges_end() ; ++ hi )
    edge2ptr[hi] = 0 ;
 
  Visitor visitor ;

  // Since the visitor is the last parameter, all othwerwise default arguments must be explicitely passed.
  // (so this examples shows the default policies).
  //  
  LindstromTurk_params params ;
    
  int r = edge_collapse(surface
                       ,Count_ratio_stop_condition<Surface>(0.10)
                       ,boost::make_assoc_property_map(edge2ptr)
                       
                       ,Vertex_is_fixed_map_always_false<Surface>()        // Same as default
                       ,Set_partial_collapse_data_LindstromTurk<Surface>() // Same as default
                       ,Cached_cost<Surface>()                             // Same as default
                       ,LindstromTurk_placement<Surface>()                 // Same as default
                       ,&params                                            // Same as default
                       ,&params                                            // Same as default
                       
                       ,&visitor
                       );

  cout << "\nFinished...\n" << r << " edges removed.\n"  << (surface.size_of_halfedges()/2) << " final edges." ;
  
  cout << "\nEdges collected: " << visitor.collected
       << "\nEdges proccessed: " << visitor.processed
       << "\nEdges collapsed: " << visitor.collapsed
       << endl
       << "\nEdges not collapsed due to topological constrians: " << visitor.non_collapsable
       << "\nEdge not collapsed due to computational constrians: " << visitor.cost_uncomputable 
       << endl ; 
        
  ofstream os( argc > 2 ? argv[2] : "out.off" ) ;
  os << surface ;
  
  return 0 ;      
}

// EOF //
