#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Unique_hash_map.h>

#include <CGAL/Surface_mesh_simplification/Polyhedron.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse.h>

// === EXAMPLE SPECIFIC HEADERS BEGINS HERE ===

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_pred.h>

// === EXAMPLE SPECIFIC HEADERS ENDS HERE ===

typedef CGAL::Simple_cartesian<double> Kernel;

// === EXAMPLE SPECIFIC DETAILS BEGINS HERE ===

//
// Setup an enriched polyhedron type which stores a ID in the halfedges
//

typedef Kernel::Point_3 Point;

struct My_items : public CGAL::Polyhedron_items_3 
{
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef CGAL::HalfedgeDS_vertex_base<Refs,CGAL::Tag_true,Point> Vertex ;
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
  
  SMS::Count_stop_condition<Surface> stop_policy(1000);
     
  // The edge_index_map() parameter is ommited becasue
  // the halfedges in this polyhedron support the stored id().
  int r = SMS::edge_collapse(surface,stop_policy);
  
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===

  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
