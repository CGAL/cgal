#include <CGAL/intersections.h>
#include <CGAL/iterator.h>
#include <boost/bind.hpp>

using namespace CGAL;

typedef Cartesian<double>    K;
typedef K::Point_2           Point;
typedef Creator_uniform_2<double,Point>  Pt_creator;
typedef K::Segment_2         Segment;

int main()
{
  std::vector<Segment> input;

  // Prepare point generator for the horizontal segment, length 200.
  typedef Random_points_on_segment_2<Point,Pt_creator> P1;
  P1 p1( Point(-100,0), Point(100,0));

  // Prepare point generator for random points on circle, radius 250.
  typedef Random_points_on_circle_2<Point,Pt_creator> P2;
  P2 p2( 250);

  // Create segments.
  typedef Creator_uniform_2< Point, Segment> Seg_creator;
  typedef Join_input_iterator_2< P1, P2, Seg_creator> Seg_iterator;
  Seg_iterator g( p1, p2);
  copy_n( g, seg_count, std::back_inserter(input));
  

  std::vector<Point> points;
  std::vector<Segment> segments;

  typedef Dispatch_output_iterator<
    cpp0x::tuple<Point,Segment>, cpp0x::tuple< std::back_insert_iterator<std::vector<Point> >,
                                               std::back_insert_iterator<std::vector<Segment> > > >
    Dispatcher;

  Dispatcher disp = dispatch_output<Point,Segment>( std::back_inserter(points),std::back_inserter(segments) );
  
  // intersections of many objects, directly dispatched
  std::transform(input.begin(), input.end(), disp,
                 // XXX fix binder 
                 boost::bind(intersection, input.front(), _1));
  std::cout << "Point intersections: " << points.size() << std::endl;
  std::cout << "Segment intersections: " << segments.size() << std::endl;

  // intersections of a single object
  IT2<K, Segment, Segment> v = intersection(input.back(), input.front());
  return 0;
}
