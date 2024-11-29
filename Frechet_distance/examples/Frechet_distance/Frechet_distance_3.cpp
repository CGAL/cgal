#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/WKT.h>

#include <ostream>
#include <fstream>

using Kernel = CGAL::Simple_cartesian<double>;
using Traits = CGAL::Frechet_distance_traits_3<Kernel>;
using Point = Traits::Point_d;

int main(int argc, char* argv[])
{
    std::vector<Point> A, B;
    {
      std::ifstream in(CGAL::data_file_path("wkt/moebius.wkt"));
      CGAL::IO::read_linestring_WKT(in, A);
    }
    {
      std::ifstream in(CGAL::data_file_path("wkt/moebius2.wkt"));
      CGAL::IO::read_linestring_WKT(in, B);
    }
    std::pair<double, double> res = CGAL::approximate_Frechet_distance<Traits>(A, B, 0.000001);
    std::cout << "The Frechet distance between the polylines is between " <<  res.first << " and " << res.second << std::endl;
    return 0;
}
