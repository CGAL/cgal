// demo/Generator/generators_prog3.C
// ------------------------------------------------------------------
// CGAL example program for point generators creating integer points.

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/Geomview_stream.h>  
#include <cassert>
#include <vector>
#include <algorithm>

#if !defined(__BORLANDC__) && !defined(_MSC_VER)

typedef CGAL::Cartesian<double>                  K;
typedef K::Point_3                               Point;
typedef CGAL::Creator_uniform_3<double,Point>    Creator;
typedef std::vector<Point>                       Vector;

int main() {
    // Create test point set. Prepare a vector for 400 points.
    Vector points;
    points.reserve(125);

    // Create 125 points from a 5 x 5 x 5 grid. Note that the double
    // arithmetic _is_ sufficient to produce exact integer grid points.
    CGAL::points_on_cube_grid_3(10, 125, std::back_inserter(points), 
                                Creator());

    // Check that we have really created 1000 points.
    assert( points.size() == 125);

    // Lower, left corner.
    assert( points[0].x() == -10);
    assert( points[0].y() == -10);
    assert( points[0].z() == -10);

    // Upper, right corner. 
    assert( points[124].x() == 10);
    assert( points[124].y() == 10);
    assert( points[124].z() == 10);

    // Visualize point set.
    CGAL::Geomview_stream geomview(CGAL::Bbox_3(-10,-10,-10, 10, 10, 10));
    for( Vector::iterator i = points.begin(); i != points.end(); i++)
        geomview << *i;

    char wait_ch;
    std::cout << "Press any key to end the program ";
    std::cin >> wait_ch;
    return 0;
}

#else // on windows:

int main() {
  std::cerr <<
  "This demo requires geomview, which is is not present on windows.\n";
  return 0;
}

#endif

