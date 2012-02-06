#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
#include <CGAL/iterator.h>
#include <CGAL/point_generators_2.h>


#include <boost/bind.hpp>

using namespace CGAL;

typedef CGAL::Simple_cartesian<double>    K;
typedef K::Point_2           Point;
typedef CGAL::Creator_uniform_2<double,Point>  Pt_creator;
typedef K::Segment_2         Segment;
typedef Random_points_on_segment_2<Point,Pt_creator> P1;
typedef Random_points_on_circle_2<Point,Pt_creator> P2;
typedef Creator_uniform_2< Point, Segment> Seg_creator;
typedef Join_input_iterator_2< P1, P2, Seg_creator> Seg_iterator;

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
  
  // intersection of two objects
#if !defined(CGAL_CFG_NO_CPP0X_AUTO)
  // with C++11
  auto v = intersection(input.back(), input.front());
#else
  // without C++11
  K::Intersect_2::Result<Segment, Segment>::Type v = intersection(input.back(), input.front());
#endif

  // splitting results with Dispatch_output_iterator
  std::vector<Point> points;
  std::vector<Segment> segments;

  typedef Dispatch_output_iterator<
    cpp0x::tuple<Point,Segment>, cpp0x::tuple< std::back_insert_iterator<std::vector<Point> >,
                                               std::back_insert_iterator<std::vector<Segment> > > >
    Dispatcher;
  
  Dispatcher disp = dispatch_output<Point,Segment>( std::back_inserter(points), 
                                                    std::back_inserter(segments) );
  
  // intersections of many objects, directly dispatched
#if !defined(CGAL_CFG_NO_CPP0X_LAMBDA)
  // with C++11 lambdas
  auto& s1 = input.front();
  std::transform(input.begin(), input.end(), disp,
                 [&s1] (const Segment& s) { return intersection(s1, s); });
#else
  // without
  std::transform(input.begin(), input.end(), disp,
                 boost::bind(static_cast<
                             K::Intersect_2::Result<Segment, Segment>::Type(*)(const Segment&, const Segment&)
                             >(&intersection),
                             input.front(), _1));
#endif

  std::cout << "Point intersections: " << points.size() << std::endl;
  std::cout << "Segment intersections: " << segments.size() << std::endl;


  return 0;
}
