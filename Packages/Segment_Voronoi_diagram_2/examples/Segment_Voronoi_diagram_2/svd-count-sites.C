// file: svd-count-sites.C
#include <CGAL/basic.h>

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// define the number type
#  include <CGAL/Quotient.h>
#  include <CGAL/MP_Float.h>
typedef CGAL::Quotient<CGAL::MP_Float>  EFT;

// define the kernel
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_exact.h>

typedef CGAL::Filtered_exact<double,EFT>  FT;
typedef CGAL::Simple_cartesian<FT>        Kernel;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Voronoi_diagram_traits_2.h>
#include <CGAL/Segment_Voronoi_diagram_2.h>

typedef CGAL::Segment_Voronoi_diagram_traits_2<Kernel>  Gt;
typedef CGAL::Segment_Voronoi_diagram_2<Gt>             SVD2;

using namespace std;

int main() {
  ifstream ifs("data/sitesx.cin");
  assert( ifs );

  SVD2          svd;
  SVD2::Site_2  site;

  while ( ifs >> site ) { svd.insert( site ); }

  ifs.close();

  assert( svd.is_valid(true, 1) );
  cout << endl << endl;

  // print the number of input and output sites
  cout << "# of input sites : " << svd.number_of_input_sites() << endl;
  cout << "# of output sites: " << svd.number_of_output_sites() << endl;

  unsigned int n_ipt(0), n_iseg(0), n_opt(0), n_oseg(0), n_ptx(0);

  // count the number of input points and input segments
  SVD2::Input_sites_iterator iit;
  for (iit = svd.input_sites_begin(); iit != svd.input_sites_end(); ++iit)
    {
      if ( iit->is_point() ) { n_ipt++; } else { n_iseg++; }
    }

  // count the number of output points and output segments, as well
  // as the number of points that are points of intersection of pairs
  // of strongly intersecting sites
  SVD2::Output_sites_iterator oit;
  for (oit = svd.output_sites_begin(); oit != svd.output_sites_end(); ++oit)
    {
      if ( oit->is_segment() ) { n_oseg++; } else {
        n_opt++;
        if ( !oit->is_input() ) { n_ptx++; }
      }
    }

  cout << endl << "# of input segments:  " << n_iseg << endl;
  cout << "# of input points:    " << n_ipt << endl << endl;
  cout << "# of output segments: " << n_oseg << endl;
  cout << "# of output points:   " << n_opt << endl << endl;
  cout << "# of intersection points: " << n_ptx << endl;

  return 0;
}
