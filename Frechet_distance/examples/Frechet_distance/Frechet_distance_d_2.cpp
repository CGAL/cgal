#include <CGAL/Epick_d.h>
#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_d.h>
#include <CGAL/IO/WKT.h>

#include <ostream>
#include <fstream>

using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<2>>;
using Traits = CGAL::Frechet_distance_traits_d<Kernel>;
using Point = Kernel::Point_d;

int main(int argc, char* argv[])
{
    std::vector<Point> polylineA, polylineB;
    {
      std::ifstream in((argc > 1) ? argv[1] : CGAL::data_file_path("wkt/LetterA.wkt"));
      //CGAL::IO::read_linestring_WKT(in, polylineA);
    }
    {
      std::ifstream in((argc > 1) ? argv[2] : CGAL::data_file_path("wkt/LetterAbis.wkt"));
      //CGAL::IO::read_linestring_WKT(in, polylineB);
    }
    bool res = CGAL::is_Frechet_distance_larger(polylineA, polylineB, 0.001);
    std::cout << std::boolalpha << res << std::endl;
    return 0;
}
