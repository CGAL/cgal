// Triangle_generator_prog1.C
// -------------------------------------------------
// CGAL example program generating random triangles.

#include <CGAL/basic.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/leda_window.h>  // used for visualization

using namespace CGAL;

typedef Cartesian<double>                R;
typedef Point_2<R>                       Point;
typedef Creator_uniform_2<double,Point>  Pt_creator;
typedef Triangle_2<R>                    Triangle;
typedef std::vector<Triangle>            Vector;

int main() {
    // Create test triangle set. Prepare a vector for 20 triangles.
    Vector triang;
    triang.reserve(20);

    // Prepare point generator for random points in a disc.
    typedef  Random_points_in_disc_2<Point,Pt_creator>  RP;
    RP p1( 1.0);
    RP p2( 1.0);
    RP p3( 1.0);

    // Create 20 triangles.
    typedef Creator_uniform_3< Point, Triangle> T_creator;
    typedef Join_input_iterator_3< RP, RP, RP, T_creator> Triang_iterator;
    Triang_iterator ti( p1, p2, p3);
    std::copy_n( ti, 20, std::back_inserter(triang));

    // Visualize triangles.
    leda_window* window = create_and_display_demo_window();
    for( Vector::iterator i = triang.begin(); i != triang.end(); i++)
        *window << *i;

    // Wait for mouse click in window.
    Point p;
    *window >> p;
    delete window;
    return 0;
}
