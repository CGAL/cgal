#include <CGAL/Epick_d.h>
#include <CGAL/Frechet_distance.h>
#include <CGAL/Frechet_distance_traits_d.h>

#include <iostream>
#include <vector>

using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<4>>;
using Traits = CGAL::Frechet_distance_traits_d<Kernel>;
using Point = Traits::Point_d;

int main(int , char*)
{
    std::array<Point,4> A = { Point(0,0,0,0), Point(1,0,0,0), Point(1,1,0,1),Point(1,1,1,0)};
    std::array<Point,4> B = { Point(0,0,0,0), Point(1,0,0,0), Point(1,1,0,0),Point(1,1,1,0)};

    std::pair<double, double> res = CGAL::approximate_Frechet_distance<Traits>(A, B, 0.000001);
    std::cout << "The Frechet distance between the polylines is between " <<  res.first << " and " << res.second << std::endl;
    return 0;
}
