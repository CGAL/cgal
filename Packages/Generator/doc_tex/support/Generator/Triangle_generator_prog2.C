/*  Triangle_generator_prog2.C      */
/*  ------------------------------- */
/*  CGAL example program generating a regular triangle pattern. */

#include <CGAL/basic.h>
#include <algo.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/Counting_iterator.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL_Cartesian<double>                            R;
typedef CGAL_Point_2<R>                                   Point;
typedef CGAL_Triangle_2<R>                                Triangle;
typedef CGAL_Points_on_segment_2<Point>                   PG;
typedef CGAL_Creator_uniform_3< Point, Triangle>          Creator;
typedef CGAL_Join_input_iterator_3< PG, PG, PG, Creator>  Triang_iterator;
typedef CGAL_Counting_iterator<Triang_iterator,Triangle>  Count_iterator;


int main()
{
    /* Prepare point generator for three segments. */
    PG p1( Point( 125, 217), Point(-125, 217), 50);
    PG p2( Point(-250,   0), Point(-125,-217), 50);
    PG p3( Point( 125,-217), Point( 250,   0), 50);

    /* Create triangle generating iterator. */
    Triang_iterator ti( p1, p2, p3);
    Count_iterator  t_begin( ti);
    Count_iterator  t_end( 50);

    /* Open window and copy 50 triangles into window. */
    CGAL_Window_stream W(512, 512);
    W.init(-256.0, 255.0, -256.0);
    W << CGAL_BLACK;
    copy( t_begin, t_end, 
	  CGAL_Ostream_iterator<Triangle,CGAL_Window_stream>(W));

    /*  Wait for mouse click in window. */
    Point p;
    W >> p;
    return 0;
}
