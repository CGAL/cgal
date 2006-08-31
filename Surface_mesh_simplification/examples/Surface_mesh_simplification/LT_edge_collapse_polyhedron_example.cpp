#include <iostream>
#include <iomanip>
#include <fstream>

#include <CGAL/Simple_cartesian.h>

// Target surface type. (this include Polyhedron_3.h itself)
#include <CGAL/Surface_mesh_simplification/Polyhedron.h>

// Policies
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

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  ifstream is(argv[1]) ;
  is >> surface ;

  //
  // The simplification algorithm needs to associate some data to each edge.
  // In this example we use an external map to do that, properly initialized with null pointers
  // and passed to the edge_collapse function as a boost associative property map
  //
  Unique_hash_map<Halfedge_handle,void*> edge2ptr ;
  for ( Surface::Halfedge_iterator hi = surface.halfedges_begin(); hi != surface.halfedges_end() ; ++ hi )
    edge2ptr[hi] = 0 ;
 
  //
  // Simplify surface down to 1000 edges using an external hasmap to store the per-edge extra pointer
  //
  // Implicit policies:
  //   No fixed vertices.
  //   LindstromTurk strategy, with default parameters.
  //   No visitor.
  //
  int r = edge_collapse(surface, Count_stop_condition<Surface>(1000), boost::make_assoc_property_map(edge2ptr) );

  cout << "\nFinished...\n" << r << " edges removed.\n"  << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  ofstream os( argc > 2 ? argv[2] : "out.off" ) ;
  os << surface ;
  
  return 0 ;      
}

// EOF //
