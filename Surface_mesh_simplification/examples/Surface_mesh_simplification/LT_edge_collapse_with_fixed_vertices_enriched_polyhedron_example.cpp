// Lindstrom-Turk edge-collapse with fixed vertices on an enriched polyhedron.
//
// Explicit arguments:
// 
//   The surface is an enriched Polyhedron_3 which embeeds the extra pointer in the edge and a boolean flag in the vertices.
//
//   The stop condition is to finish when the number of undirected edges 
//   drops below a certain absolute number.
//
//   The per-edge extra pointer is stored in the edge itself.
//
//   The flag determining whether the vertex is fixed is stored in the vertex itself
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
#include <CGAL/Surface_mesh_simplification/Vertex_is_fixed_map_stored.h>  //<==== NOTICE THIS 
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
// Setup an enriched polyhedron type which stores in a halfedge each extra pointer needed by the algorithm
// and in each vertex whether it is fixed or not.
//

template <class Refs, class Traits>
struct My_vertex : public HalfedgeDS_vertex_base<Refs,Tag_true,Point>
{ 
  typedef HalfedgeDS_vertex_base<Refs,Tag_true,Point> Base ;
  
  My_vertex() : is_fixed_(false) {}
  
  My_vertex( Point const& p ) : Base(p), is_fixed_(false) {}
 
  bool is_fixed() const { return is_fixed_ ; }
    
  bool is_fixed_ ;
};

template <class Refs, class Traits>
struct My_halfedge : public HalfedgeDS_halfedge_base<Refs> 
{ 
  My_halfedge() : extra_ptr_(0) {}
 
  void*& extra_pointer() { return extra_ptr_ ; }
    
  void* extra_ptr_ ;
};

struct My_items : public Polyhedron_items_3 
{
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef My_vertex<Refs,Traits>  Vertex;
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

typedef Surface::Halfedge_handle Halfedge_handle ;
typedef Surface::Vertex_handle   Vertex_handle ;

using namespace CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;


int main( int argc, char** argv ) 
{
  Surface surface; 
  
  ifstream is(argv[1]) ;
  is >> surface ;

  for ( Surface::Vertex_iterator vi = surface.vertices_begin(); vi != surface.vertices_end() ; ++ vi )
    vi->is_fixed_ = true ; // only some would be set to true, of cotrue
  
  int r = edge_collapse(surface
                       ,Count_ratio_stop_condition<Surface>(0.10)
                       ,Edge_extra_pointer_map_stored<Surface>()
                       ,Vertex_is_fixed_map_stored<Surface>()
                       );

  cout << "\nFinished...\n" << r << " edges removed.\n"  << (surface.size_of_halfedges()/2) << " final edges." ;
  
        
  ofstream os( argc > 2 ? argv[2] : "out.off" ) ;
  os << surface ;
  
  return 0 ;      
}

// EOF //
