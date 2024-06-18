#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_d.h>
#include <CGAL/Epick_d.h>

#include <iostream>
#include <vector>

using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<4>>;
using Traits = CGAL::Frechet_distance_traits_d<Kernel>;
using Point = Traits::Point;

int main(int argc, char* argv[])
{
    std::vector<Point> A, B;
    bool res = CGAL::continuous_Frechet_distance_less_than<Traits>(A, B, 0.001);
    return 0;
}
