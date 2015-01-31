// standard includes

#include <iostream>
#include <fstream>
#include <cassert>

#include "svd-typedefs.h"

using namespace std;

int main( int argc, char *argv[] ) {
  if ( not (( argc == 1 ) or (argc == 2)) ) {
    std::cout <<"usage: "<< argv[0] <<" [filename]\n";
  }

  ifstream ifs( (argc == 1) ? "data/sitesx.cin" : argv[1] );
  assert( ifs );

  SDG_2          sdg;
  SDG_2::Site_2  site;

  while ( ifs >> site ) { sdg.insert( site ); }

  ifs.close();

  // print the number of input and output sites
  cout << "# of input sites : " << sdg.number_of_input_sites() << endl;
  cout << "# of output sites: " << sdg.number_of_output_sites() << endl;

  unsigned int n_ipt(0), n_iseg(0), n_opt(0), n_oseg(0), n_ptx(0);

  // count the number of input points and input segments
  SDG_2::Input_sites_iterator iit;
  for (iit = sdg.input_sites_begin(); iit != sdg.input_sites_end(); ++iit)
    {
      if ( iit->is_point() ) { n_ipt++; } else { n_iseg++; }
    }

  // count the number of output points and output segments, as well
  // as the number of points that are points of intersection of pairs
  // of strongly intersecting sites
  SDG_2::Output_sites_iterator oit;
  for (oit = sdg.output_sites_begin(); oit != sdg.output_sites_end(); ++oit)
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
