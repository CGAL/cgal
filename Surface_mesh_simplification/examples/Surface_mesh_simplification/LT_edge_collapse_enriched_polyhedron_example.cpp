// Lindstrom-Turk edge-collapse with enriched surface for embeeded per-edge pointer.
//
// Explicit arguments:
// 
//   The surface is an enriched Polyhedron_3 which embeeds the extra pointer in the edge.
//
//   The stop condition is to finish when the number of undirected 
//   edges drops below a certain relative fraction of the initial count.
//
// Implicit arguments:
// 
//   The per-edge extra pointer is given by the edge itself.
//   No vertex is fixed.
//   The cost strategy is Lindstrom-Turk
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

//
// Setup an enriched polyhedron type which stores in a halfedge the extra pointer needed by the algorithm
//

template <class Refs, class Traits>
struct My_halfedge : public HalfedgeDS_halfedge_base<Refs> 
{ 
  My_halfedge() : extra_ptr_(0) {} // The extra pointer MUST be initialized to NULL
 
  void*& extra_pointer() { return extra_ptr_ ; }
    
  void* extra_ptr_ ;
};

struct My_items : public Polyhedron_items_3 
{
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef CGAL::HalfedgeDS_vertex_base<Refs,Tag_true,Point> Vertex ;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper { 
        typedef My_halfedge<Refs,Traits>  Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef CGAL::HalfedgeDS_face_base<Refs,CGAL::Tag_true> Face;
    };
};

typedef Polyhedron_3<Kernel,My_items> Surface; 


using namespace CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  ifstream is(argv[1]) ;
  is >> surface ;

  int r = edge_collapse(surface, Count_ratio_stop_condition<Surface>(0.10) );

  cout << "\nFinished...\n" << r << " edges removed.\n"  << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  ofstream os( argc > 2 ? argv[2] : "out.off" ) ;
  os << surface ;
  
  return 0 ;      
}

// EOF //
