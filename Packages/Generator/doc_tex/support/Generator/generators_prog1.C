/*  generators_prog1.C              */
/*  ------------------------------- */
/*  CGAL example program for point generators. */

#include <CGAL/basic.h>
#include <assert.h>
#include <vector.h>
#include <algo.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/random_selection.h>
#include <CGAL/IO/leda_window.h>  /* only used for visualization */

typedef CGAL_Cartesian<double>                R;
typedef CGAL_Point_2<R>                       Point;
typedef CGAL_Creator_uniform_2<double,Point>  Creator;

int main()
{
    /* Create test point set. Prepare a vector for 1000 points. */
    vector<Point> points;
    points.reserve(1000);

    /* Create 600 points within a disc of radius 0.6 */
    CGAL_Random_points_in_disc_2<Point,Creator> g( 0.6);
    CGAL_copy_n( g, 600, back_inserter( points));

    /* Create 200 points from a 15 x 15 grid. */
    CGAL_points_on_square_grid_2(0.95, 200, back_inserter(points),Creator());

    /* Select 100 points randomly and append them at the end of */
    /* the current vector of points. */
    CGAL_random_selection( points.begin(), points.end(), 100, 
			   back_inserter( points));

    /* Create 100 points that are collinear to two randomly chosen */
    /* points and append them to the current vector of points. */
    CGAL_random_collinear_points_2( points.begin(), points.end(), 100,
				    back_inserter( points));

    /* Check that we have really created 1000 points. */
    assert( points.size() == 1000);

    /* A random permutation to hide the creation order in the point set. */
    random_shuffle( points.begin(), points.end(), CGAL_random);

    /* Visualize point set. */
    leda_window* window = CGAL_create_and_display_demo_window();
    for( vector<Point>::iterator i = points.begin(); i != points.end(); i++)
	*window << *i;

    /*  Wait for mouse click in window. */
    Point p;
    *window >> p;
    delete window;
    return 0;
}
