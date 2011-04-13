// examples/Generator/Segment_generator_example2.C
// --------------------------------------------------------------
// CGAL example program for the generic segment generator
// using precomputed point locations.
// Derived from Segment_generator_prog2.C, without visualization.

#include <CGAL/basic.h>
#include <algorithm>
#include <vector>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/Counting_iterator.h>

using namespace CGAL;

typedef Cartesian<double>                            R;
typedef Point_2<R>                                   Point;
typedef Segment_2<R>                                 Segment;
typedef Points_on_segment_2<Point>                   PG;
typedef Creator_uniform_2< Point, Segment>           Creator;
typedef Join_input_iterator_2< PG, PG, Creator>      Segm_iterator;
typedef Counting_iterator<Segm_iterator,Segment>     Count_iterator;
typedef std::vector<Segment>             Vector;

int main() {
    // Create test segment set. Prepare a vector for 100 segments.
    Vector segs;
    segs.reserve(100);

    // A horizontal like fan.
    PG p1( Point(-250, -50), Point(-250, 50),50);   // Point generator.
    PG p2( Point( 250,-250), Point( 250,250),50);
    Segm_iterator  t1( p1, p2);                     // Segment generator.
    Count_iterator t1_begin( t1);                   // Finite range.
    Count_iterator t1_end( 50);
    std::copy( t1_begin, t1_end, std::back_inserter(segs));

    // A vertical like fan.
    PG p3( Point( -50,-250), Point(  50,-250),50);
    PG p4( Point(-250, 250), Point( 250, 250),50);
    Segm_iterator  t2( p3, p4);
    Count_iterator t2_begin( t2);
    Count_iterator t2_end( 50);
    std::copy( t2_begin, t2_end, std::back_inserter(segs));

    CGAL_assertion( segs.size() == 100);
    for ( Vector::iterator i = segs.begin(); i != segs.end(); i++){
	CGAL_assertion( i->source().x() <=  250);
	CGAL_assertion( i->source().x() >= -250);
	CGAL_assertion( i->source().y() <=  250);
	CGAL_assertion( i->source().y() >= -250);
	CGAL_assertion( i->target().x() <=  250);
	CGAL_assertion( i->target().x() >= -250);
	CGAL_assertion( i->target().y() <=  250);
	CGAL_assertion( i->target().y() >= -250);
    }
    return 0;
}


