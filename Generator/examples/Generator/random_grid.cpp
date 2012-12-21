#include <CGAL/Simple_cartesian.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

using namespace CGAL;

typedef Simple_cartesian<int>         K;
typedef K::Point_2                    Point;
typedef Creator_uniform_2<int,Point>  Creator;

int main() {
    // Create test point set. Prepare a vector for 400 points.
    std::vector<Point> points;
    points.reserve(400);

    // Create 250 points from a 16 x 16 grid. Note that the double
    // arithmetic _is_ sufficient to produce exact integer grid points.
    // The distance between neighbors is 34 pixel = 510 / 15.
    points_on_square_grid_2( 255.0, 250, std::back_inserter(points),Creator());

    // Lower, left corner.
    assert( points[0].x() == -255);
    assert( points[0].y() == -255);

    // Upper, right corner. Note that 6 points are missing to fill the grid.
    assert( points[249].x() == 255 - 6 * 34);
    assert( points[249].y() == 255);

    // Create 250 points within a disc of radius 150.
    Random_points_in_disc_2<Point,Creator> g( 150.0);
    CGAL::cpp11::copy_n( g, 250, std::back_inserter(points));

    // Check that we have really created 500 points.
    assert( points.size() == 500);
    return 0;
}
