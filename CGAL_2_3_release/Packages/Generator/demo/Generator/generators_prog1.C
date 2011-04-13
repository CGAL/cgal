// demo/Generator/generators_prog1.C
// ------------------------------------------
// CGAL example program for point generators.

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/random_selection.h>
#include <CGAL/IO/Window_stream.h>  // used for visualization
#include <cassert>
#include <vector>
#include <algorithm>

typedef CGAL::Cartesian<double>                K;
typedef K::Point_2                             Point;
typedef CGAL::Creator_uniform_2<double,Point>  Creator;

typedef std::vector<Point>                     Vector;

int main() {
    // Create test point set. Prepare a vector for 1000 points.
    Vector points;
    points.reserve(1000);

    // Create 600 points within a disc of radius 0.6
    CGAL::Random_points_in_disc_2<Point,Creator> g( 0.6);
    CGAL::copy_n( g, 600, std::back_inserter(points));

    // Create 200 points from a 15 x 15 grid.
    CGAL::points_on_square_grid_2(0.95, 200, std::back_inserter(points),
                                  Creator());

    // Select 100 points randomly and append them at the end of
    // the current vector of points.
    CGAL::random_selection( points.begin(), points.end(), 100,
                            std::back_inserter(points));

    // Create 100 points that are collinear to two randomly chosen
    // points and append them to the current vector of points.
    CGAL::random_collinear_points_2( points.begin(), points.end(), 100,
                                     std::back_inserter(points));

    // Check that we have really created 1000 points.
    assert( points.size() == 1000);

    // A random permutation to hide the creation order in the point set.
    std::random_shuffle( points.begin(), points.end(), CGAL::default_random);

    // Visualize point set.
    CGAL::Window_stream* window = CGAL::create_and_display_demo_window();
    for( Vector::iterator i = points.begin(); i != points.end(); i++)
        *window << *i;

    // Wait for mouse click in window.
    (*window).read_mouse();
    delete window;
    return 0;
}
