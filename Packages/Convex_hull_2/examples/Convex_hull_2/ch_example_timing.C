#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_real.h>
#endif // CGAL_USE_LEDA

#include <CGAL/convex_hull_traits_2.h>

#include <fstream>

#include <deque>
#include <list>
#include <vector>

#include <CGAL/ch_akl_toussaint.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/ch_eddy.h>
#include <CGAL/ch_bykat.h>
#include <CGAL/ch_jarvis.h>

#include <CGAL/ch_timing_2.h>

using namespace std;

typedef double                                      nu_type;
typedef CGAL::Cartesian< nu_type >                  TraitsCls;
typedef TraitsCls::Point_2                          Point2;

int
main( int argc, char* argv[] )
{
  if (argc != 3)   // assertion
  {
      cerr << "Usage: ch_example_timing datafilename ";
      cerr << "number_of_iterations";
      exit(1);
  }
  vector< Point2 > V;
  vector< Point2 > VE;
  ifstream F(argv[1]);
  CGAL::set_ascii_mode( F );
  istream_iterator< Point2>  in_start( F );
  istream_iterator< Point2>  in_end;
  copy( in_start, in_end , back_inserter(V) );
  copy( V.begin(), V.end(), back_inserter(VE) );
  int iterations = atoi( argv[2] );
  CGAL::ch_timing(V.begin(), V.end(), VE.begin(), iterations, TraitsCls() ); 
  return 0;
}
