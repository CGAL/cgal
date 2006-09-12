// Lindstrom-Turk edge-collapse with fixed vertices in an ordinary polyhedron.
//
// Explicit arguments:
// 
//   The surface is an ordinary polyhedron.
//
//   The stop condition is to finish when the number of undirected edges 
//   drops below a certain absolute number.
//
//   An external map is used to store the per-edge extra pointer unintrusively.
//
//   An external map is used to specify which vertices are fixed.
//
// Implicit arguments:
// 
//   The cost strategy is Lindstrom-Turk with partial cache.
//   No visitor is passed.
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

typedef Surface::Halfedge_handle     Halfedge_handle ;
typedef Surface::Vertex_const_handle Vertex_const_handle ;

using namespace CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;


int main( int argc, char** argv ) 
{
  Surface surface; 
  
  ifstream is(argv[1]) ;
  is >> surface ;

  // Extra pointer external map
  Unique_hash_map<Halfedge_handle,void*> edge2ptr ;
  for ( Surface::Halfedge_iterator hi = surface.halfedges_begin(); hi != surface.halfedges_end() ; ++ hi )
    edge2ptr[hi] = 0 ;
 
  //
  // The following Map specifies which vertices are fixed.  
  //
  Unique_hash_map<Vertex_const_handle,bool> is_vertex_fixed ;
  
  for ( Surface::Vertex_const_iterator vi = surface.vertices_begin(); vi != surface.vertices_end() ; ++ vi )
    is_vertex_fixed[vi] = false ; // some would be set to true
 
  int r = edge_collapse(surface
                       ,Count_ratio_stop_condition<Surface>(0.10)
                       ,boost::make_assoc_property_map(edge2ptr)
                       ,boost::make_assoc_property_map(is_vertex_fixed)
                       );

  cout << "\nFinished...\n" << r << " edges removed.\n"  << (surface.size_of_halfedges()/2) << " final edges." ;
        
  ofstream os( argc > 2 ? argv[2] : "out.off" ) ;
  os << surface ;
  
  return 0 ;      
}

// EOF //
