#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/tuple.h>
#include <fstream>
#include <iostream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;

int main(int argc, char** argv)
{
  if ( argc!=3){
    std::cerr << "Usage: " << argv[0] <<  " intput.off profile.txt\n";
    return EXIT_FAILURE;
  }
  
  std::ifstream input(argv[1]);
  std::ofstream output(argv[2]);

  Polyhedron_3 p;
  CGAL::scan_OFF( input, p,true);
  if ( !p.is_closed () ) std::cout << "WARNING the polyhedron is not closed!!!!" << std::endl;
  
  std::multimap<double,Polyhedron_3::Halfedge_handle> distance_map;
  
  for (Polyhedron_3::Edge_iterator it=p.edges_begin();it!=p.edges_end();++it)
  {
    double sq_dist = CGAL::squared_distance( it->vertex()->point(),it->opposite()->vertex()->point() );
    distance_map.insert( std::make_pair( sq_dist,Polyhedron_3::Halfedge_handle(it) ) );
  }
  
  output.precision(45);
  for(std::multimap<double,Polyhedron_3::Halfedge_handle>::iterator it=distance_map.begin();it!=distance_map.end();++it)
  {
    output << it->first << " " << it->second->vertex()->point() << " " << it->second->opposite()->vertex()->point() << "\n";
  }
  
  
}
  

