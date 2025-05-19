#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/iterator.h>
#include <CGAL/point_generators_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2                     Point;
typedef K::Segment_2                   Segment;

typedef CGAL::Creator_uniform_2<double,Point>              Pt_creator;
typedef CGAL::Random_points_on_segment_2<Point,Pt_creator> P1;
typedef CGAL::Random_points_on_circle_2<Point,Pt_creator>  P2;
typedef CGAL::Creator_uniform_2< Point, Segment>           Seg_creator;
typedef CGAL::Join_input_iterator_2< P1, P2, Seg_creator>  Seg_iterator;

struct Intersector{
  const Segment& s;
  K::Intersect_2 intersect;

  Intersector(const Segment& seg): s(seg) {}

  decltype(auto)
  operator() ( const Segment& other) const
  {
    return intersect(s, other);
  }
};

int main()
{
  std::vector<Segment> input;

  // Prepare point generator for the horizontal segment, length 200.
  P1 p1( Point(-100,0), Point(100,0));

  // Prepare point generator for random points on circle, radius 250.
  P2 p2( 250);

  // Create segments.
  Seg_iterator g( p1, p2);
  std::copy_n( g, 200, std::back_inserter(input));


  // splitting results with Dispatch_output_iterator
  std::vector<Point> points;
  std::vector<Segment> segments;

  typedef CGAL::Dispatch_output_iterator<
    std::tuple<Point,Segment>, std::tuple< std::back_insert_iterator<std::vector<Point> >,
                                               std::back_insert_iterator<std::vector<Segment> > > >
    Dispatcher;

  Dispatcher disp = CGAL::dispatch_output<Point,Segment>( std::back_inserter(points),
                                                          std::back_inserter(segments) );

  // intersects the first segment of input with all other segments
  // The resulting points or segments are written in the vectors with the same names
  std::transform( input.begin(), input.end(), disp,
                  Intersector(input.front()) );


  std::cout << "Point intersections: " << points.size() << std::endl;
  std::cout << "Segment intersections: " << segments.size() << std::endl;


  return 0;
}
