/*  generators_prog1.C              */
/*  ------------------------------- */
/*  CGAL example program for point generators. */

#include <assert.h>
#include <vector.h>
#include <algo.h>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/ncopy.h>
#include <CGAL/random_selection.h>
#include <CGAL/IO/Window_stream.h>  /* only for visualization used */

typedef CGAL_Cartesian<double>  R;
typedef CGAL_Point_2<R>         Point;

int main()
{
    /* Create test point set. Prepare a vector for 1000 points. */
    vector<Point> points;
    points.reserve(1000);

    /* Create 600 points within a disc of radius 150. */
    CGAL_Random_points_in_disc_2<Point> g( 150.0);
    CGAL_ncopy( g, 600, back_inserter( points));

    /* Create 200 points from a 15 x 15 grid. */
    CGAL_points_on_square_grid_2( 500.0, 200, back_inserter(points),(Point*)0);

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

    /* Use a random permutation to hide the creation history */
    /* of the point set. */
    random_shuffle( points.begin(), points.end(), CGAL_random);

    /* Visualize point set. Can be omitted, see example programs */
    /* in the CGAL source code distribution. */
    CGAL_Window_stream W(512, 512);
    W.init(-256.0, 255.0, -256.0);
    W << CGAL_BLACK;
    for( vector<Point>::iterator i = points.begin(); i != points.end(); i++)
	W << *i;

    /*  Wait for program termination. */
    char c;
    cout << " Type any character to continue: " << endl;
    cin >> c;
    cout << " done" << endl;

    return 0;
}
