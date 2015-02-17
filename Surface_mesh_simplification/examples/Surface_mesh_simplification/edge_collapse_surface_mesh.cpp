#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h> 

typedef CGAL::Simple_cartesian<double> Kernel;

typedef Kernel::Point_3 Point_3;

//
// Setup an enriched polyhedron type which stores an id() field in the items
//
typedef CGAL::Surface_mesh<Point_3> Surface_mesh; 

typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor ;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification ;

typedef SMS::Edge_profile<Surface_mesh> Profile ;


// The following is a Visitor that keeps track of the simplification process.
// In this example the progress is printed real-time and a few statistics are
// recorded (and printed in the end).
//
struct Stats
{
  Stats() 
    : collected(0)
    , processed(0)
    , collapsed(0)
    , non_collapsable(0)
    , cost_uncomputable(0) 
    , placement_uncomputable(0) 
  {} 
  
  std::size_t collected ;
  std::size_t processed ;
  std::size_t collapsed ;
  std::size_t non_collapsable ;
  std::size_t cost_uncomputable  ;
  std::size_t placement_uncomputable ; 
} ;

struct My_visitor : SMS::Edge_collapse_visitor_base<Surface_mesh>
{
  My_visitor( Stats* s) : stats(s){} 

  // Called during the collecting phase for each edge collected.
  void OnCollected( Profile const&, boost::optional<double> const& )
  {
    ++ stats->collected ;
    std::cerr << "\rEdges collected: " << stats->collected << std::flush ;
  }                
  
  // Called during the processing phase for each edge selected.
  // If cost is absent the edge won't be collapsed.
  void OnSelected(Profile const&          
                 ,boost::optional<double> cost
                 ,std::size_t             initial
                 ,std::size_t             current
                 )
  {
    ++ stats->processed ;
    if ( !cost )
      ++ stats->cost_uncomputable ;
      
    if ( current == initial )
      std::cerr << "\n" << std::flush ;
    std::cerr << "\r" << current << std::flush ;
  }                
  
  // Called during the processing phase for each edge being collapsed.
  // If placement is absent the edge is left uncollapsed.
  void OnCollapsing(Profile const&          
                   ,boost::optional<Point>  placement
                   )
  {
    if ( !placement )
      ++ stats->placement_uncomputable ;
  }                
  
  // Called for each edge which failed the so called link-condition,
  // that is, which cannot be collapsed because doing so would
  // turn the surface mesh into a non-manifold.
  void OnNonCollapsable( Profile const& )
  {
    ++ stats->non_collapsable;
  }                
  
  // Called AFTER each edge has been collapsed
  void OnCollapsed( Profile const&, vertex_descriptor )
  {
    ++ stats->collapsed;
  }                
  
  Stats* stats ;
} ;


int main( int argc, char** argv ) 
{
  Surface_mesh surface_mesh; 
  
  std::ifstream is(argv[1]) ; is >> surface_mesh ;

  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(0.1);
 
  Stats stats ;
  
  My_visitor vis(&stats) ;
    
  // The index maps are not explicitelty passed as in the previous
  // example because the surface mesh items have a proper id() field.
  // On the other hand, we pass here explicit cost and placement
  // function which differ from the default policies, ommited in
  // the previous example.
  int r = SMS::edge_collapse
           (surface_mesh
           ,stop
            ,CGAL::parameters::get_cost     (SMS::Edge_length_cost  <Surface_mesh>())
                              .get_placement(SMS::Midpoint_placement<Surface_mesh>())
                              .visitor      (vis)
           );
  
  std::cout << "\nEdges collected: "  << stats.collected
            << "\nEdges proccessed: " << stats.processed
            << "\nEdges collapsed: "  << stats.collapsed
            << std::endl
            << "\nEdges not collapsed due to topological constraints: "  << stats.non_collapsable
            << "\nEdge not collapsed due to cost computation constraints: "  << stats.cost_uncomputable 
            << "\nEdge not collapsed due to placement computation constraints: " << stats.placement_uncomputable 
            << std::endl ; 
            
  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << num_edges(surface_mesh) << " final edges.\n" ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface_mesh ;
  
  return 0 ;      
}

// EOF //
