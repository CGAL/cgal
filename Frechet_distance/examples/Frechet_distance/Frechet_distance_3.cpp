#include <CGAL/Frechet_distance.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/WKT.h>

#include <ostream>
#include <fstream>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;

int main(int argc, char* argv[])
{
    std::vector<Point> polylineA, polylineB;
    {
      std::ifstream in((argc > 1) ? argv[1] : CGAL::data_file_path("wkt/moebius.wkt"));
      CGAL::IO::read_linestring_WKT(in, polylineA);
    }
    {
      std::ifstream in((argc > 2) ? argv[2] : CGAL::data_file_path("wkt/moebius2.wkt"));
      CGAL::IO::read_linestring_WKT(in, polylineB);
    }

    std::pair<double, double> res = CGAL::bounded_error_Frechet_distance(polylineA, polylineB, 0.000001);
    std::cout << "The Frechet distance between the polylines is between " <<  res.first << " and " << res.second << std::endl;
    return 0;
}
