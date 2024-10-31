#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <ostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Traits = CGAL::Frechet_distance_traits_2<Kernel>;
using Point = Traits::Point_d;

int main(int argc, char* argv[])
{
    std::vector<Point> A, B;
    bool res = CGAL::is_Frechet_distance_larger<Traits>(A, B, 0.001);
    std::cout << std::boolalpha << res << std::endl;
    return 0;
}
