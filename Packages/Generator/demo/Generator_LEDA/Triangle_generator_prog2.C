//  Triangle_generator_prog2.C
//  -------------------------------
//  CGAL example program generating a regular triangle pattern.

#include <CGAL/basic.h>
#ifndef CGAL_USE_LEDA
int main() { std::cout << "\nSorry, this demo needs LEDA\n"; return 0; }
#else
#include <algorithm>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/Counting_iterator.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/leda_window.h>  // only used for visualization

using namespace CGAL;

typedef Cartesian<double>                            R;
typedef Point_2<R>                                   Point;
typedef Triangle_2<R>                                Triangle;
typedef Points_on_segment_2<Point>                   PG;
typedef Creator_uniform_3< Point, Triangle>          Creator;
typedef Join_input_iterator_3< PG, PG, PG, Creator>  Triang_iterator;
typedef Counting_iterator<Triang_iterator,Triangle>  Count_iterator;

int main() {
    // Prepare point generator for three segments.
    PG p1( Point(  0.50,  0.90), Point( -0.50,  0.90), 50);
    PG p2( Point( -0.95,  0.00), Point( -0.50, -0.90), 50);
    PG p3( Point(  0.50, -0.90), Point(  0.95,  0.00), 50);

    // Create triangle generating iterator.
    Triang_iterator ti( p1, p2, p3);
    Count_iterator  t_begin( ti);
    Count_iterator  t_end( 50);

    // Open window and copy 50 triangles into window.
    leda_window* window = create_and_display_demo_window();
    std::copy( t_begin, t_end,
               Ostream_iterator<Triangle,Window_stream>(*window));

    // Wait for mouse click in window.
    (*window).read_mouse();
    delete window;
    return 0;
}
#endif // CGAL_USE_LEDA
