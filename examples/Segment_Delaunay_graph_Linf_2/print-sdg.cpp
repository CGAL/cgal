// standard includes

//#define CGAL_SDG_VERBOSE 

#ifndef CGAL_SDG_VERBOSE
#define CGAL_SDG_DEBUG(a)
#else
#define CGAL_SDG_DEBUG(a) { a }
#endif

#include <iostream>
#include <fstream>
#include <cassert>
#include <string>

// define the kernel
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

typedef CGAL::Simple_cartesian<double>    CK;
typedef CGAL::Filtered_kernel<CK>         Kernel;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_Linf_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_hv_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2.h>

typedef CGAL::Segment_Delaunay_graph_Linf_traits_2<Kernel>  Gt;
typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt>             SDG2;

typedef CGAL::Segment_Delaunay_graph_Linf_hv_traits_2<Kernel>  Gt_hv;
typedef CGAL::Segment_Delaunay_graph_Linf_2<Gt_hv>             SDG2_HV;

using namespace std;

int main( int argc, char *argv[] ) {
  bool use_hv = false;
  int fileat = 0;
  if ( argc >= 4 ) {
    std::cout <<"usage: "<< argv[0] <<" [filename]\n" <<
    "[-h] option for sdg_hv" << std::endl;
  }
  use_hv = (argc == 1) ? false :
  (argc == 2 and (argv[1][0] == '-' and argv[1][1] == 'h') ) ? true :
  (argc == 3 and ( (argv[1][0] == '-' and argv[1][1] == 'h')
                 or(argv[2][0] == '-' and argv[2][1] == 'h') )) ? true : false;
  fileat = (argc == 1) ? 0 :
  (argc == 2) ? (use_hv == true ? 0 : 1) :
  (argc == 3 and (argv[1][0] == '-' and argv[1][1] == 'h')) ? 1 : 2;
  
  ifstream ifs( (fileat == 0) ? "data/sites2.cin" :
                (fileat == 1) ? argv[1] : argv[2] );
  assert( ifs );

  SDG2          sdg;
  SDG2_HV       sdg_hv;
  SDG2::Site_2  site;

  // read the sites from the stream and insert them in the diagram
  if( use_hv ) {
    while ( ifs >> site ) {
      sdg_hv.insert( site );
    }
  } else {
    while ( ifs >> site ) {
      sdg.insert( site );
    }
  }
  
  ifs.close();


  //std::cout << "About to validate diagram ..." << std::endl;

  // validate the diagram
  //assert( sdg.is_valid(true, 1) );
  //cout << endl << endl;

  //std::cout << "Diagram validated." << std::endl;
  (use_hv) ?
  std::cout << "About to print sdg_hv for input file: " << ((fileat == 0) ?
               "data/sites2.cin" : (fileat == 1) ? argv[1] : argv[2]) << std::endl :
  std::cout << "About to print sdg for input file: " << ((fileat == 0) ?
               "data/sites2.cin" : (fileat == 1) ? argv[1] : argv[2]) << std::endl ;
  sdg.file_output_verbose(std::cout);
  
  return 0;
}
