// stream_iterator.C
// ----------------------------------------------------------
// CGAL example program for the CGAL stream iterator adaptor.

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/Istream_iterator.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/leda_window.h>
#include <algorithm>

using namespace CGAL;

typedef Cartesian<double>                       TutorialR;
typedef Point_2<TutorialR>                      Point;
typedef Creator_uniform_2<double,Point>         Creator;
typedef Random_points_in_disc_2<Point,Creator>  Random_points_in_disc;

void init_window( leda_window& W) {
    cgalize( W);
    W.set_fg_color( leda_green);
    W.display();
    W.init(-1.0, 1.0, -1.0);
}

int main()
{
    Point  points[100];

    // Create 100 random points uniform distributed in a disc of radius 1.
    // Use deterministic initialization for the random number generator.
    Random rnd(1);
    Random_points_in_disc  rnd_points( 1.0, rnd);
    CGAL::copy_n( rnd_points, 100, points);

    // Display points in a 512x512 pixel window.
    leda_window W(512, 512);
    init_window( W);
    std::copy( points, points+100, Ostream_iterator<Point, Window_stream>(W));

    // Wait for mouse click in window.
    Istream_iterator<Point, Window_stream>  si(W);
    Point q = *si;  // W >> q;

    return 0;
}

