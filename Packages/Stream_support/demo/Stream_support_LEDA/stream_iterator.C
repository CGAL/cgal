#ifndef CGAL_USE_LEDA
#include <iostream>
int main(){ std::cout << "This demo needs LEDA" << std::endl; return 0;}
#else


// stream_iterator.C
// ----------------------------------------------------------
// CGAL example program for the CGAL stream iterator adaptor.

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/Istream_iterator.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/Window_stream.h>
#include <algorithm>

typedef CGAL::Cartesian<double>                       TutorialR;
typedef TutorialR::Point_2                            Point;
typedef CGAL::Creator_uniform_2<double,Point>         Creator;
typedef CGAL::Random_points_in_disc_2<Point,Creator>  Random_points_in_disc;

#ifdef CGAL_USE_CGAL_WINDOW
#define leda_window CGAL::window
#define leda_green  CGAL::green
#endif

void init_window( leda_window& W)
{
    CGAL::cgalize( W);
    W.set_fg_color( leda_green);
    W.display();
    W.init(-1.0, 1.0, -1.0);
}

int main() {
    Point points[100];

    // Create 100 random points uniform distributed in a disc of radius 1.
    // Use deterministic initialization for the random number generator.
    CGAL::Random rnd(1);
    Random_points_in_disc  rnd_points( 1.0, rnd);
    CGAL::copy_n( rnd_points, 100, points);

    // Display points in a 512x512 pixel window.
    leda_window W(512, 512);
    init_window( W);
    std::copy( points, points+100,
               CGAL::Ostream_iterator<Point, CGAL::Window_stream>(W));

    // Wait for mouse click in window.
    CGAL::Istream_iterator<Point, CGAL::Window_stream>  si(W);
    Point q = *si;  // W >> q;

    return 0;
}

#endif
