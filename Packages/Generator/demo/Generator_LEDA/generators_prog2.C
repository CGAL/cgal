// generators_prog2.C
// ------------------------------------------------------------------
// CGAL example program for point generators creating integer points.

#include <CGAL/basic.h>
#ifndef CGAL_USE_LEDA
int main() { std::cout << "\nSorry, this demo needs LEDA\n"; return 0; }
#else
#include <cassert>
#include <vector>
#include <algorithm>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/leda_window.h>  // used for visualization

using namespace CGAL;

typedef Cartesian<double>                  R;
typedef Point_2<R>                         Point;
typedef Creator_uniform_2<double,Point>    Creator;
typedef std::vector<Point>                 Vector;

int main() {
    // Create test point set. Prepare a vector for 400 points.
    Vector points;
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
    CGAL::copy_n( g, 250, std::back_inserter(points));

    // Check that we have really created 500 points.
    assert( points.size() == 500);

    // Visualize point set.
    leda_window* window = create_and_display_demo_window();
    window->init(-262.0, 261.0, -262.0);
    for( Vector::iterator i = points.begin(); i != points.end(); i++)
        *window << *i;

    //  Wait for mouse click in window.
    (*window).read_mouse();
    delete window;
    return 0;
}
#endif // CGAL_USE_LEDA
