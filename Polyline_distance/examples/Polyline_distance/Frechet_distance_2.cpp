#include <CGAL/Frechet_distance.h>
#include <CGAL/Polyline_distance_traits_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <ostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Traits = CGAL::Polyline_distance_traits_2<Kernel>;
using Point = Traits::Point;

int main(int argc, char* argv[])
{
    std::vector<Point> A, B;
    bool res = CGAL::continuous_Frechet_distance_less_than<Traits>(A, B, 0.001);
    return 0;
}
