#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Unique_hash_map.h>

#include <CGAL/Surface_mesh_simplification/Polyhedron.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse.h>

// === EXAMPLE SPECIFIC POLICIES BEGINS HERE ===
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_pred.h>
// === EXAMPLE SPECIFIC POLICIES ENDS HERE ===

typedef CGAL::Simple_cartesian<double> Kernel;


//
// Setup an enriched polyhedron type which stores in the halfedge the extra pointer needed by the algorithm
//

template <class Refs, class Traits>
struct My_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs> 
{ 
  My_halfedge() : extra_ptr_(0) {} // The extra pointer MUST be initialized to NULL
 
  void*& extra_pointer() { return extra_ptr_ ; }
    
  void* extra_ptr_ ;
};

typedef Kernel::Point_3 Point;

struct My_items : public CGAL::Polyhedron_items_3 
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

typedef CGAL::Polyhedron_3<Kernel,My_items> Surface; 

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  CGAL::Triangulated_surface_mesh
      ::Simplification
      ::Edge_collapse
      ::Count_stop_condition<Surface> stop_policy(1000);
     
  // The third argument is ommited this time because the default policy
  // is to put the extra-pointer in the edge itself.
  int r = CGAL::Triangulated_surface_mesh
              ::Simplification
              ::Edge_collapse
              ::edge_collapse(surface
                             ,stop_policy
                             );

  std::cout << "\nFinished...\n" << r << " edges removed.\n"  << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
