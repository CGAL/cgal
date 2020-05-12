#include <CGAL/Simple_cartesian.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/random_selection.h>

using namespace CGAL;

typedef Simple_cartesian<double>         R;
typedef R::Point_2                       Point;
typedef Creator_uniform_2<double,Point>  Creator;
typedef std::vector<Point>               Vector;

int main() {
    // Create test point set. Prepare a vector for 1000 points.
    Vector points;
    points.reserve(1000);

    // Create 600 points within a disc of radius 150.
    Random_points_in_disc_2<Point,Creator> g( 150.0);
    std::copy_n( g, 600, std::back_inserter(points));

    // Create 200 points from a 15 x 15 grid.
    points_on_square_grid_2( 250.0, 200, std::back_inserter(points),Creator());

    // Select 100 points randomly and append them at the end of
    // the current vector of points.
    random_selection( points.begin(), points.end(), 100,
                      std::back_inserter(points));

    // Create 100 points that are collinear to two randomly chosen
    // points and append them to the current vector of points.
    random_collinear_points_2( points.begin(), points.end(), 100,
                               std::back_inserter( points));

    // Check that we have really created 1000 points.
    assert( points.size() == 1000);

    // Use a random permutation to hide the creation history
    // of the point set.
    CGAL::cpp98::random_shuffle( points.begin(), points.end());

    // Check range of values.
    for ( Vector::iterator i = points.begin(); i != points.end(); i++){
        assert( i->x() <=  251);
        assert( i->x() >= -251);
        assert( i->y() <=  251);
        assert( i->y() >= -251);
    }
    return 0;
}
