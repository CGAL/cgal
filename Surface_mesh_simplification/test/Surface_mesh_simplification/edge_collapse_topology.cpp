#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface; 

namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv ) 
{
  if (argc!=2){
    std::cerr << "Please provide only an off-file as input\n";
    return 1;
  }

  std::ifstream is(argv[1]);

  if (!is){
    std::cerr << "Error reading the input\n";
    return 1;
  }

  Surface surface;
  is >> surface ;

  std::size_t initial_count = (surface.size_of_halfedges()/2);
  std::cout << "Initial count " << initial_count << " edges.\n" ;

  // Contract the surface as much as possible
  SMS::Count_stop_predicate<Surface> stop(0);

  int r = SMS::edge_collapse
            (surface
            ,stop
             ,CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,surface)) 
                               .halfedge_index_map  (get(CGAL::halfedge_external_index  ,surface)) 
            );

  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;

  assert( initial_count == (surface.size_of_halfedges()/2) + r );
  // std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;

  return 0 ;      
}

// EOF //
