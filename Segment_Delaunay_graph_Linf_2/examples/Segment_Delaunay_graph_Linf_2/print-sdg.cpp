// standard includes
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

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
    std::cout << "usage: " << argv[0] << " [filename]" << std::endl
              << "       -h: option for sdg_hv" << std::endl;
  }

  use_hv =
    (argc == 1) ?
      false :
      (argc == 2 and (strlen(argv[1]) == 2) and
       (argv[1][0] == '-' and argv[1][1] == 'h') ) ?
        true :
        (argc == 3 and
         (((strlen(argv[1]) == 2) and
            argv[1][0] == '-' and argv[1][1] == 'h') or
          ((strlen(argv[2]) == 2) and
           argv[2][0] == '-' and argv[2][1] == 'h')   )) ?
          true :
          false;

  fileat = (argc == 1) ? 0 :
    (argc == 2) ? (use_hv == true ? 0 : 1) :
    (argc == 3 and (argv[1][0] == '-' and argv[1][1] == 'h')) ? 2 : 1;

  std::cout << "use_hv = " << use_hv
            << " fileat: " << fileat << std::endl;

  ifstream ifs( (fileat == 0) ? "data/sites2.cin" : argv[fileat] );
  assert( ifs );

  SDG2          sdg;
  SDG2_HV       sdg_hv;
  SDG2::Site_2  site;

  // read the sites from the stream and insert them in the diagram
  while ( ifs >> site ) {
    if( use_hv ) {
      if (site.is_segment()) {
        if ( not( site.segment().is_horizontal() or
                  site.segment().is_vertical() ) ) {
          std::cout << "input is not axis parallel " << site << endl;
          return 1;
        }
      }
      sdg_hv.insert( site );
    } else {
      sdg.insert( site );
    }
  }

  ifs.close();


  //std::cout << "About to validate diagram ..." << std::endl;

  // validate the diagram
  //assert( sdg.is_valid(true, 1) );
  //cout << endl << endl;

  //std::cout << "Diagram validated." << std::endl;
  std::cout
     << "About to print " << ((use_hv) ? "sdg_hv" : "sdg")
     << " for input file: "
     << ((fileat == 0) ? "data/sites2.cin" : argv[fileat])
     << std::endl ;

  if (use_hv) {
    sdg_hv.file_output_verbose(std::cout);
  } else {
    sdg.file_output_verbose(std::cout);
  }

  return 0;
}
