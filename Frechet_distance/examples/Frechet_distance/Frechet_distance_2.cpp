#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Frechet_distance.h>
#include <CGAL/IO/WKT.h>

#include <ostream>
#include <fstream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_2;

int main(int argc, char* argv[])
{
    std::vector<Point> A, B;
    {
      std::ifstream in((argc > 1) ? argv[1] : CGAL::data_file_path("wkt/LetterA.wkt"));
      CGAL::IO::read_linestring_WKT(in, A);
    }
    {
      std::ifstream in((argc > 1) ? argv[2] : CGAL::data_file_path("wkt/LetterAbis.wkt"));
      CGAL::IO::read_linestring_WKT(in, B);
    }
    bool res = CGAL::is_Frechet_distance_larger(A, B, 0.001);
    std::cout << std::boolalpha << res << std::endl;
    return 0;
}
