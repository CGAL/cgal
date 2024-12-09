#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/WKT.h>

#include <ostream>
#include <fstream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
//using Kernel = CGAL::Simple_cartesian<double>;
using Traits = CGAL::Frechet_distance_traits_3<Kernel>;
using Point = Traits::Point_d;

int main(int argc, char* argv[])
{
#if 0
    std::vector<Point> A, B;
    {
      std::ifstream in(CGAL::data_file_path("wkt/moebius.wkt"));
      CGAL::IO::read_linestring_WKT(in, A);
    }
    {
      std::ifstream in(CGAL::data_file_path("wkt/moebius2.wkt"));
      CGAL::IO::read_linestring_WKT(in, B);
    }
    #else

    // two identical points fails
    std::array<Point,4> A = { Point(0,0,0), Point(0,0,0), Point(1,0,1), Point(1,1,0)};
    std::array<Point,4> B = { Point(0,0,0), Point(0,0,0), Point(1,0,0), Point(1,1,0)};

#endif

    std::pair<double, double> res = CGAL::approximate_Frechet_distance<Traits>(A, B, 0.000001);
    std::cout << "The Frechet distance between the polylines is between " <<  res.first << " and " << res.second << std::endl;
    return 0;
}
