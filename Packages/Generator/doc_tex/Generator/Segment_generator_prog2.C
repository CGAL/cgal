/*  Segment_generator_prog2.C       */
/*  ------------------------------- */
/*  CGAL example program generating a regular segment pattern. */

#include <CGAL/basic.h>
#include <algo.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/Counting_iterator.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL_Cartesian<double>                            R;
typedef CGAL_Point_2<R>                                   Point;
typedef CGAL_Segment_2<R>                                 Segment;
typedef CGAL_Points_on_segment_2<Point>                   PG;
typedef CGAL_Creator_uniform_2< Point, Segment>           Creator;
typedef CGAL_Join_input_iterator_2< PG, PG, Creator>      Segm_iterator;
typedef CGAL_Counting_iterator<Segm_iterator,Segment>     Count_iterator;


int main()
{
    /* Open window. */
    CGAL_Window_stream W(512, 512);
    W.init(-256.0, 255.0, -256.0);
    W << CGAL_BLACK;

    /* A horizontal like fan. */
    PG p1( Point(-250, -50), Point(-250, 50),50);     /* Point generator. */
    PG p2( Point( 250,-250), Point( 250,250),50);
    Segm_iterator  t1( p1, p2);                       /* Segment generator. */
    Count_iterator t1_begin( t1);                     /* Finite range. */
    Count_iterator t1_end( 50);
    copy( t1_begin, t1_end, 
	  CGAL_Ostream_iterator<Segment,CGAL_Window_stream>(W));

    /* A vertical like fan. */
    PG p3( Point( -50,-250), Point(  50,-250),50);
    PG p4( Point(-250, 250), Point( 250, 250),50);
    Segm_iterator  t2( p3, p4);
    Count_iterator t2_begin( t2);
    Count_iterator t2_end( 50);
    copy( t2_begin, t2_end, 
	  CGAL_Ostream_iterator<Segment,CGAL_Window_stream>(W));

    /*  Wait for mouse click in window. */
    Point p;
    W >> p;
    return 0;
}
