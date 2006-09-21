#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Surface_mesh_simplification/Polyhedron.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse.h>


// === EXAMPLE SPECIFIC HEADERS BEGINS HERE ===

#include <CGAL/Surface_mesh_simplification/Vertex_is_fixed_map_stored.h>  

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_pred.h>

// === EXAMPLE SPECIFIC HEADERS ENDS HERE ===


typedef CGAL::Simple_cartesian<double> Kernel;

// === EXAMPLE SPECIFIC DETAILS BEGINS HERE ===

typedef Kernel::Point_3 Point;

//
// Setup an enriched polyhedron type which stores in a halfedge
// each extra pointer needed by the algorithm
// and in each vertex whether it is fixed or not.
//

template <class Refs, class Traits>
struct My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs,CGAL::Tag_true,Point>
{ 
  typedef CGAL::HalfedgeDS_vertex_base<Refs,CGAL::Tag_true,Point> Base ;
  
  My_vertex() : is_fixed_(false) {}
  
  My_vertex( Point const& p ) : Base(p), is_fixed_(false) {}
 
  bool is_fixed() const { return is_fixed_ ; }
    
  bool is_fixed_ ;
};

struct My_items : public CGAL::Polyhedron_items_3 
{
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef My_vertex<Refs,Traits>  Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper { 
        typedef CGAL::HalfedgeDS_halfedge_max_base_with_id<Refs> Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef CGAL::HalfedgeDS_face_base<Refs,CGAL::Tag_true> Face;
    };
};

typedef CGAL::Polyhedron_3<Kernel,My_items> Surface; 

// === EXAMPLE SPECIFIC DETAILS ENDS HERE ===

namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  // In this example, the polyhedron has been enriched to store a flag 
  // right in the vertex that indicates if the vertex is fixed.
  //
  // This irrealistic loop just illustrates how the flag would be set.
  //
  for ( Surface::Vertex_iterator vi = surface.vertices_begin()
      ; vi != surface.vertices_end() 
      ; ++ vi 
      )
    vi->is_fixed_ = true ; // only some would be set to true, of cotrue
  
  int r = SMS::edge_collapse(surface
                             ,SMS::Count_ratio_stop_condition<Surface>(0.10) 
                             ,SMS.vertex_is_fixed_map(CGAL::Vertex_is_fixed_map_stored<Surface>())
                             // The edge_index_map parameter is ommited becasue
                             // the halfedge in this polyhedron supports the stored id().
                             );
           
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===
  
  std::cout << "\nFinished...\n" << r << " edges removed.\n"  << (surface.size_of_halfedges()/2) << " final edges." ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
