/*  generators_prog2.C              */
/*  ------------------------------- */
/*  CGAL example program for point generators creating integer points. */

#include <CGAL/basic.h>
#include <assert.h>
#include <vector.h>
#include <algo.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/leda_window.h>  /* only used for visualization */

typedef CGAL_Cartesian<double>                R;
typedef CGAL_Point_2<R>                       Point;
typedef CGAL_Creator_uniform_2<double,Point>  Creator;

int main()
{
    /* Create test point set. Prepare a vector for 400 points. */
    vector<Point> points;
    points.reserve(400);

    /* Create 250 points from a 16 x 16 grid. Note that the double */
    /* arithmetic _is_ sufficient to produce exact integer grid points. */
    /* The distance between neighbors is 34 pixel = 510 / 15. */
    CGAL_points_on_square_grid_2( 255.0, 250, back_inserter(points),Creator());

    /* Lower, left corner. */
    assert( points[0].x() == -255);
    assert( points[0].y() == -255);

    /* Upper, right corner. Note that 6 points are missing to fill the grid. */
    assert( points[249].x() == 255 - 6 * 34);
    assert( points[249].y() == 255);

    /* Create 250 points within a disc of radius 150. */
    CGAL_Random_points_in_disc_2<Point,Creator> g( 150.0);
    CGAL_copy_n( g, 250, back_inserter( points));

    /* Check that we have really created 500 points. */
    assert( points.size() == 500);

    /* Visualize point set. */
    leda_window* window = CGAL_create_and_display_demo_window();
    window->init(-262.0, 261.0, -262.0);
    for( vector<Point>::iterator i = points.begin(); i != points.end(); i++)
	*window << *i;

    /*  Wait for mouse click in window. */
    Point p;
    *window >> p;
    delete window;
    return 0;
}
