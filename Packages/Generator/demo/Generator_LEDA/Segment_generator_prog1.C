// Segment_generator_prog1.C
// -------------------------------
// CGAL example program for the generic segment generator.

#include <CGAL/basic.h>
#ifndef CGAL_USE_LEDA
int main() { std::cout << "\nSorry, this demo needs LEDA\n"; return 0; }
#else
#include <cassert>
#include <vector>
#include <algorithm>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/leda_window.h>  // used for visualization

using namespace CGAL;

typedef Cartesian<double>                R;
typedef Point_2<R>                       Point;
typedef Creator_uniform_2<double,Point>  Pt_creator;
typedef Segment_2<R>                     Segment;
typedef std::vector<Segment>             Vector;

int main() {
    // Create test segment set. Prepare a vector for 200 segments.
    Vector segs;
    segs.reserve(200);

    // Prepare point generator for the horizontal segment, length 200.
    typedef  Random_points_on_segment_2<Point,Pt_creator>  P1;
    P1 p1( Point( -0.4, 0), Point( 0.4, 0));

    // Prepare point generator for random points on circle, radius 250.
    typedef  Random_points_on_circle_2<Point,Pt_creator>  P2;
    P2 p2( 1.0);

    // Create 200 segments.
    typedef Creator_uniform_2< Point, Segment> Seg_creator;
    typedef Join_input_iterator_2< P1, P2, Seg_creator> Seg_iterator;
    Seg_iterator g( p1, p2);
    CGAL::copy_n( g, 200, std::back_inserter(segs));

    // Visualize segments.
    leda_window* window = create_and_display_demo_window();
    for( Vector::iterator i = segs.begin(); i != segs.end(); i++)
        *window << *i;

    //  Wait for mouse click in window.
    (*window).read_mouse();
    delete window;
    return 0;
}
#endif // CGAL_USE_LEDA
