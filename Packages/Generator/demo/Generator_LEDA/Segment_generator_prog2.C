// Segment_generator_prog2.C
// -------------------------------
// CGAL example program generating a regular segment pattern.

#include <CGAL/basic.h>
#ifndef CGAL_USE_LEDA
int main() { std::cout << "\nSorry, this demo needs LEDA\n"; return 0; }
#else
#include <algorithm>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/Counting_iterator.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/leda_window.h>  // used for visualization

using namespace CGAL;

typedef Cartesian<double>                          R;
typedef Point_2<R>                                 Point;
typedef Segment_2<R>                               Segment;
typedef Points_on_segment_2<Point>                 PG;
typedef Creator_uniform_2< Point, Segment>         Creator;
typedef Join_input_iterator_2< PG, PG, Creator>    Segm_iterator;
typedef Counting_iterator<Segm_iterator,Segment>   Count_iterator;

int main() {
    // Open window.
    leda_window* window = create_and_display_demo_window();
    window->init(-256.0, 255.0, -256.0);

    // A horizontal like fan.
    PG p1( Point(-250, -50), Point(-250, 50),50);     // Point generator.
    PG p2( Point( 250,-250), Point( 250,250),50);
    Segm_iterator  t1( p1, p2);                       // Segment generator.
    Count_iterator t1_begin( t1);                     // Finite range.
    Count_iterator t1_end( 50);
    std::copy( t1_begin, t1_end, 
               Ostream_iterator<Segment,Window_stream>(*window));

    // A vertical like fan.
    PG p3( Point( -50,-250), Point(  50,-250),50);
    PG p4( Point(-250, 250), Point( 250, 250),50);
    Segm_iterator  t2( p3, p4);
    Count_iterator t2_begin( t2);
    Count_iterator t2_end( 50);
    std::copy( t2_begin, t2_end,
               Ostream_iterator<Segment,Window_stream>(*window));

    // Wait for mouse click in window.
    (*window).read_mouse();
    delete window;
    return 0;
}
#endif // CGAL_USE_LEDA
