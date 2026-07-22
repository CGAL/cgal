#include <CGAL/Epick_d.h>
#include <CGAL/Frechet_distance.h>

#include <iostream>
#include <vector>

using Kernel = CGAL::Epick_d<CGAL::Dimension_tag<4>>;
using Point = Kernel::Point_d;

int main()
{
    std::array<Point,4> polylineA = { Point(0,0,0,0), Point(1,0,0,0), Point(1,1,0,1),Point(1,1,1,0)};
    std::array<Point,4> polylineB = { Point(0,0,0,0), Point(1,0,0,0), Point(1,1,0,0),Point(1,1,1,0)};

    std::pair<double, double> res = CGAL::bounded_error_Frechet_distance(polylineA, polylineB, 0.000001);
    std::cout << "The Frechet distance between the polylines is between " <<  res.first << " and " << res.second << std::endl;
    return 0;
}
