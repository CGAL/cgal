// demo/Generator/Triangle_generator_prog1.C
// -------------------------------------------------
// CGAL example program generating random triangles.

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/Window_stream.h>
#include <cassert>
#include <vector>
#include <algorithm>

typedef CGAL::Cartesian<double>                K;
typedef K::Point_2                             Point;
typedef CGAL::Creator_uniform_2<double,Point>  Pt_creator;
typedef K::Triangle_2                          Triangle;
typedef std::vector<Triangle>                  Vector;

int main() {
    // Create test triangle set. Prepare a vector for 20 triangles.
    Vector triang;
    triang.reserve(20);

    // Prepare point generator for random points in a disc.
    typedef  CGAL::Random_points_in_disc_2<Point,Pt_creator>  RP;
    RP p1( 1.0);
    RP p2( 1.0);
    RP p3( 1.0);

    // Create 20 triangles.
    typedef CGAL::Creator_uniform_3< Point, Triangle> T_creator;
    typedef CGAL::Join_input_iterator_3< RP, RP, RP, T_creator> Triang_iterator;
    Triang_iterator ti( p1, p2, p3);
    CGAL::copy_n( ti, 20, std::back_inserter(triang));

    // Visualize triangles.
    CGAL::Window_stream* window = CGAL::create_and_display_demo_window();
    for( Vector::iterator i = triang.begin(); i != triang.end(); i++)
        *window << *i;

    // Wait for mouse click in window.
    (*window).read_mouse();
    delete window;
    return 0;
}

