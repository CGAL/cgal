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
  if (argc < 1 || argc > 3)
  {
      std::cerr << "Usage: " << argv[0] << " [data_file_name = files/CD500] ";
      std::cerr << "[number_of_iterations = 10]" << std::endl;
      std::exit(1);
  }

  std::ifstream F( (argc >= 2) ? argv[1] : "files/CD500");
  CGAL::set_ascii_mode( F );
  std::istream_iterator< Point_2>  in_start( F );
  std::istream_iterator< Point_2>  in_end;

  std::vector< Point_2 > V (in_start, in_end);
  std::vector< Point_2 > VE = V;

  int iterations = (argc == 3) ? std::atoi( argv[2] ) : 10;

  CGAL::ch_timing(V.begin(), V.end(), VE.begin(), iterations, K() );

  return 0;
}
