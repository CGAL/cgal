#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/algorithm.h>

using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                       Point;
typedef Creator_uniform_2<double,Point>  Pt_creator;
typedef K::Segment_2                     Segment;
typedef std::vector<Segment>             Vector;


int main() {
    // Create test segment set. Prepare a vector for 200 segments.
    Vector segs;
    segs.reserve(200);

    // Prepare point generator for the horizontal segment, length 200.
    typedef  Random_points_on_segment_2<Point,Pt_creator>  P1;
    P1 p1( Point(-100,0), Point(100,0));

    // Prepare point generator for random points on circle, radius 250.
    typedef  Random_points_on_circle_2<Point,Pt_creator>  P2;
    P2 p2( 250);

    // Create 200 segments.
    typedef Creator_uniform_2< Point, Segment> Seg_creator;
    typedef Join_input_iterator_2< P1, P2, Seg_creator> Seg_iterator;
    Seg_iterator g( p1, p2);
    std::copy_n( g, 200, std::back_inserter(segs));

    assert( segs.size() == 200);
    for ( Vector::iterator i = segs.begin(); i != segs.end(); i++){
        assert( i->source().x() <=  100);
        assert( i->source().x() >= -100);
        assert( i->source().y() ==    0);
        assert( i->target().x() * i->target().x() +
                i->target().y() * i->target().y() <=  251*251);
        assert( i->target().x() * i->target().x() +
                i->target().y() * i->target().y() >=  249*249);
    }
    return 0;
}
