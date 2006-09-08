// Edge-collapse with ordinary surface and external per-edge pointer.
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
//   Both the cost and placement are cached in the collapse_data.
//   The example shows two choices: following the Lindstrom-Turk strategy or the Edge-length/Midpoint strategy.
//
//   The on-deman cost dirtecty given by the collapse data.
// 
//   The on-deman placement is directly given by the collapse data.
//
// Implicit arguments:
//
//   No visitor is passed.
//
#include <iostream>
#include <iomanip>
#include <fstream>

#include <CGAL/Simple_cartesian.h>

// Target surface type. (this include Polyhedron_3.h itself)
#include <CGAL/Surface_mesh_simplification/Polyhedron.h>

// Policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>  //<==== NOTICE THIS 
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_pred.h>

// Simplification method
#include <CGAL/Surface_mesh_simplification/Edge_collapse.h>

#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Unique_hash_map.h>

using namespace std ;
using namespace boost ;
using namespace CGAL ;

typedef Simple_cartesian<double> Kernel;
typedef Kernel::Vector_3         Vector;
typedef Kernel::Point_3          Point;

typedef Polyhedron_3<Kernel> Surface; 
typedef Surface::Halfedge_handle Halfedge_handle ;

using namespace CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;

enum { LT, MP } ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  int strategy = atoi(argv[1]);
  
  ifstream is(argv[2]) ;
  is >> surface ;

  // Extra pointer external map
  Unique_hash_map<Halfedge_handle,void*> edge2ptr ;
  for ( Surface::Halfedge_iterator hi = surface.halfedges_begin(); hi != surface.halfedges_end() ; ++ hi )
    edge2ptr[hi] = 0 ;
    
  if ( strategy == LT )
  {
    //
    // Set_full_collapse_data_LindstromTurk is the one actually computing and caching the cost and placement
    // following the Lindstrom-Turk strategy. This, Cached_cost and Cached_placement as passed since they
    // only need to extract the cached values from the collapse-data.
    // 
    edge_collapse(surface
                 ,Count_stop_condition<Surface>(1000)
                 ,boost::make_assoc_property_map(edge2ptr)
                 
                 ,Vertex_is_fixed_map_always_false<Surface>() // Same as default
                 
                 ,Set_full_collapse_data_LindstromTurk<Surface>()
                 ,Cached_cost<Surface>()
                 ,Cached_placement<Surface>()
                 );
  }
  else
  {
    //
    // In this case, Set_full_collapse_data<> is the one actually computing and caching the cost and placement.
    // But by itself it doesn't know how to do that, so it delegates the computation to external policies.
    // 
    // That is, you need to construct GetCost and GetPlacement policies which actually compute the values
    // and passed them to the constructor of Set_full_collapse_data<>.
    //
    // In this example the strategy is Midpoint-Length.
    // 
    // Still, the GetCost and GetPlacement policies **passed to the algorithm** are
    // Cached_cost and Cached_placement as they only need to extract the cached values 
    // from the collapse-data.
    //
    typedef Edge_length_cost  <Surface> Compute_cost ;
    typedef Midpoint_placement<Surface> Compute_placement ;
    
    typedef Set_full_collapse_data<Surface,Compute_cost,Compute_placement> Set_full_collapse_data ;
    
    Compute_cost      compute_cost ;
    Compute_placement compute_placement ;
    Set_full_collapse_data set_collapse_data(compute_cost,compute_placement);
    
    edge_collapse(surface
                 ,Count_stop_condition<Surface>(1000)
                 ,boost::make_assoc_property_map(edge2ptr)
                 
                 ,Vertex_is_fixed_map_always_false<Surface>() // Same as default
                 
                 ,set_collapse_data
                 ,Cached_cost<Surface>()
                 ,Cached_placement<Surface>()
                 );
  }  

  cout << "\nFinished: " << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  ofstream os( argc > 3 ? argv[3] : "out.off" ) ;
  os << surface ;
  
  return 0 ;      
}

// EOF //
