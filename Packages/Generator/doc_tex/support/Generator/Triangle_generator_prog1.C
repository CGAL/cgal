/*  Triangle_generator_prog1.C      */
/*  ------------------------------- */
/*  CGAL example program generating random triangles. */

#include <CGAL/basic.h>
#include <assert.h>
#include <vector.h>
#include <algo.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/Window_stream.h>  /* only for visualization used */

typedef CGAL_Cartesian<double>                R;
typedef CGAL_Point_2<R>                       Point;
typedef CGAL_Creator_uniform_2<double,Point>  Pt_creator;
typedef CGAL_Triangle_2<R>                    Triangle;

int main()
{
    /* Create test segment set. Prepare a vector for 20 triangles. */
    vector<Triangle> triang;
    triang.reserve(20);

    /* Prepare point generator for random points in a disc. */
    typedef  CGAL_Random_points_in_disc_2<Point,Pt_creator>  RP;
    RP p1( 250);
    RP p2( 250);
    RP p3( 250);
    
    /* Create 20 triangles. */
    typedef CGAL_Creator_uniform_3< Point, Triangle> T_creator;
    typedef CGAL_Join_input_iterator_3< RP, RP, RP, T_creator> Triang_iterator;
    Triang_iterator ti( p1, p2, p3);
    CGAL_copy_n( ti, 20, back_inserter( triang));

    /* Visualize triangles. Can be omitted, see example programs */
    /* in the CGAL source code distribution. */
    CGAL_Window_stream W(512, 512);
    W.init(-256.0, 255.0, -256.0);
    W << CGAL_BLACK;
    for( vector<Triangle>::iterator i = triang.begin(); i != triang.end(); i++)
	W << *i;

    /*  Wait for mouse click in window. */
    Point p;
    W >> p;

    return 0;
}
