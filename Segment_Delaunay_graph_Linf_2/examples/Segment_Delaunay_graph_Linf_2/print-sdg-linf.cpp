// standard includes
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

// define the kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Kernel>  Gt;
typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt>             SDG2;

using namespace std;

int main( int argc, char *argv[] ) {
  if ( argc >= 3 ) {
    std::cout << "usage: " << argv[0] << " [filename]" << std::endl;
  }

  ifstream ifs( (argc == 1) ? "data/sites2.cin" : argv[1] );
  assert( ifs );

  SDG2          sdg;
  SDG2::Site_2  site;

  // read the sites from the stream and insert them in the diagram
  while ( ifs >> site ) {
    sdg.insert( site );
  }

  ifs.close();

  //std::cout << "About to validate diagram ..." << std::endl;

  // validate the diagram
  //assert( sdg.is_valid(true, 1) );
  //cout << endl << endl;

  //std::cout << "Diagram validated." << std::endl;
  std::cout
     << "About to print sdg for input file: "
     << ((argc == 1) ? "data/sites2.cin" : argv[1])
     << std::endl ;

  sdg.file_output_verbose(std::cout);

  return 0;
}
