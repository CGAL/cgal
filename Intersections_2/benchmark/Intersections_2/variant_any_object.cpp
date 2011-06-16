
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/intersections.h>
#include <CGAL/tuple.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Join_input_iterator.h>

#include <vector>
#include <functional>

#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/any.hpp>
#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>


template<typename Ret, typename K>
struct segment_point_2_visitor : public boost::static_visitor<Ret>
{
  segment_point_2_visitor(std::function<Ret(const typename K::Segment_2&)>&& t, 
			  std::function<Ret(const typename K::Point_2&)>&& u) 
    : t(t), u(u) {}

  const std::function<Ret(const typename K::Segment_2&)> t;
  const std::function<Ret(const typename K::Point_2&)> u;

  Ret operator()(const typename K::Segment_2& x) { return t(x); }
  Ret operator()(const typename K::Point_2& x) { return u(x); }
};

template<typename, typename, typename>
struct Intersection_traits;

template<typename K>
struct Intersection_traits<K, typename K::Segment_2, typename K::Segment_2> {
  typedef typename boost::optional< boost::variant<typename K::Segment_2, typename K::Point_2 > > result_type;
};


template <class K>
boost::optional< 
  boost::variant<typename K::Segment_2, typename K::Point_2>
  >
intersection_variant(const typename K::Segment_2 &seg1, 
	     const typename K::Segment_2 &seg2,
	     const K&)
{
  typedef CGAL::internal::Segment_2_Segment_2_pair<K> is_t;

  typedef boost::variant<typename K::Segment_2, typename K::Point_2> Variant;
  typedef boost::optional<Variant> OptVariant;

  is_t ispair(&seg1, &seg2);
  switch (ispair.intersection_type()) {
  case is_t::NO_INTERSECTION:
  default:
    return OptVariant();
  case is_t::POINT:
    return OptVariant(ispair.intersection_point());
  case is_t::SEGMENT:
    return OptVariant(ispair.intersection_segment());
  }
}

template <class K>
boost::any intersection_any(const typename K::Segment_2 &seg1, 
			    const typename K::Segment_2 &seg2,
			    const K&)
{
  typedef CGAL::internal::Segment_2_Segment_2_pair<K> is_t;

  is_t ispair(&seg1, &seg2);
  switch (ispair.intersection_type()) {
  case is_t::NO_INTERSECTION:
  default:
    return boost::any();
  case is_t::POINT:
    return boost::any(ispair.intersection_point());
  case is_t::SEGMENT:
    return boost::any(ispair.intersection_segment());
  }
}


using namespace CGAL;

typedef Cartesian<double>    K;
typedef K::Point_2           Point;
typedef Creator_uniform_2<double,Point>  Pt_creator;
typedef K::Segment_2         Segment;
typedef std::vector<Segment> Vector;


cpp0x::tuple<int, int, int> intersect_each_object(const Vector& segs) {
  // Calculate the intersections between each segment
  cpp0x::tuple<int, int, int> ret = cpp0x::make_tuple(0, 0, 0);
  for(Vector::const_iterator it = segs.begin(); it != segs.end(); ++it) {
    const Segment& seg_1 = *it;
    for(Vector::const_iterator it2 = segs.begin(); it2 != segs.end(); ++it2) {
      Object obj = intersection(seg_1, *it2);
      if (const Point * point = object_cast<Point>(&obj)) {
	(void)point;
	++(cpp0x::get<0>(ret));
      } else if (const Segment * segment = object_cast<Segment>(&obj)) {
	(void)segment;
	++(cpp0x::get<1>(ret));
      } else {
	++(cpp0x::get<2>(ret));
      }
    }
  }

  return ret;
}

struct raise_if : public boost::static_visitor<>
{
  raise_if(cpp0x::tuple<int, int, int>* t) : t(t) {}

  void operator()(const Point& p) const { (void)p; ++(cpp0x::get<0>(*t)); }  
  void operator()(const Segment& s) const { (void)s; ++(cpp0x::get<1>(*t)); }
  cpp0x::tuple<int, int, int>* t;
};

cpp0x::tuple<int, int, int> intersect_each_variant(const Vector& segs) {
  cpp0x::tuple<int, int, int> ret = cpp0x::make_tuple(0, 0, 0);
  typedef Intersection_traits<K, Segment, Segment>::result_type result_type;

  // Calculate the intersections between each segment
  for(Vector::const_iterator it = segs.begin(); it != segs.end(); ++it) {
    const Segment& seg_1 = *it;
    for(Vector::const_iterator it2 = segs.begin(); it2 != segs.end(); ++it2) {
      result_type obj = intersection_variant(seg_1, *it2, K());
      if(obj) {
 	// with c++0x
	// segment_point_2_visitor<void, K> visitor([&ret](const Segment& s) { ++(cpp0x::get<1>(ret)); },
	// 					 [&ret](const Point& s) { ++(cpp0x::get<0>(ret)); });
	// boost::apply_visitor(visitor, *obj);

	// without c++0x
	boost::apply_visitor(raise_if(&ret), *obj);
      } else {
	++(cpp0x::get<2>(ret));
      }
    }
  }
  
  return ret;
}

cpp0x::tuple<int, int, int> intersect_each_any(const Vector& segs) {
  cpp0x::tuple<int, int, int> ret = cpp0x::make_tuple(0, 0, 0);
  // Calculate the intersections between each segment
  for(Vector::const_iterator it = segs.begin(); it != segs.end(); ++it) {
    const Segment& seg_1 = *it;
    for(Vector::const_iterator it2 = segs.begin(); it2 != segs.end(); ++it2) {
      boost::any obj = intersection_any(seg_1, *it2, K());
      if (const Point * point = boost::any_cast<Point>(&obj)) {
	(void)point;
	++(cpp0x::get<0>(ret));
      } else if (const Segment * segment = boost::any_cast<Segment>(&obj)) {
	(void)segment;
	++(cpp0x::get<1>(ret));
      } else {
	++(cpp0x::get<2>(ret));
      }
      
    }
  }
  
  return ret;
}



int main(int argc, char* argv[]) {
  int repeats = 100;
  int segments = 200;

  if(argc > 1)
    repeats = boost::lexical_cast<int>(argv[1]);
  if(argc > 2)
    repeats = boost::lexical_cast<int>(argv[2]);

  // Create test segment set. Prepare a vector for 200 segments.
  Vector segs;
  segs.reserve(200);

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
  CGAL::copy_n( g, segments, std::back_inserter(segs));

  cpp0x::tuple<int, int, int> t(0, 0, 0);

  boost::timer timer;

  // variant vs object vs any

  timer.restart();
  for(int i = 0; i < repeats; ++i)
    t = intersect_each_object(segs);
  std::cout << "Time for object: " << timer.elapsed() << '\n';
  std::cout << cpp0x::get<0>(t) << " " << cpp0x::get<1>(t) << " " << cpp0x::get<2>(t) << '\n';

  t = cpp0x::make_tuple(0, 0, 0);

  timer.restart();
  for(int i = 0; i < repeats; ++i)
    t = intersect_each_variant(segs);
  std::cout << "Time for variant: " << timer.elapsed() << '\n';
  std::cout << cpp0x::get<0>(t) << " " << cpp0x::get<1>(t) << " " << cpp0x::get<2>(t) << '\n';

  t = cpp0x::make_tuple(0, 0, 0);

  timer.restart();
  for(int i = 0; i < repeats; ++i)
    t = intersect_each_any(segs);
  std::cout << "Time for any: " << timer.elapsed() << '\n';
  std::cout << cpp0x::get<0>(t) << " " << cpp0x::get<1>(t) << " " << cpp0x::get<2>(t) << '\n';


  std::cout << std::flush;

}
