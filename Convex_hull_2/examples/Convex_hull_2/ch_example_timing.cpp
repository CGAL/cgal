#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_traits_2.h>
#include <CGAL/ch_akl_toussaint.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/ch_eddy.h>
#include <CGAL/ch_bykat.h>
#include <CGAL/ch_jarvis.h>
#include <CGAL/ch_timing_2.h>

#include <fstream>
#include <vector>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

int main( int argc, char* argv[] )
{
  if (argc < 2 || argc > 3)   // assertion
  {
      std::cerr << "Usage: ch_example_timing data_file_name ";
      std::cerr << "[number_of_iterations]" << std::endl;
      std::exit(1);
  }
  std::vector< Point_2 > V;
  std::vector< Point_2 > VE;
  std::ifstream F(argv[1]);
  CGAL::set_ascii_mode( F );
  std::istream_iterator< Point_2>  in_start( F );
  std::istream_iterator< Point_2>  in_end;
  std::copy( in_start, in_end , std::back_inserter(V) );
  std::copy( V.begin(), V.end(), std::back_inserter(VE) );
  int iterations;
  if (argc == 3)
    iterations = std::atoi( argv[2] );
  else
    iterations = 1;
  CGAL::ch_timing(V.begin(), V.end(), VE.begin(), iterations, K() );
  return 0;
}
