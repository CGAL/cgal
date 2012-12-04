#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/intersections.h>
#include <CGAL/tuple.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/iterator.h>

#include <vector>
#include <functional>

#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/any.hpp>
#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>

#if !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE)
#include <tuple>
#include <functional>
#include <CGAL/Overload.h>
#endif

// Intersection_traits
template<typename, typename, typename>
struct Intersection_traits;

template<typename K>
struct Intersection_traits<K, typename K::Segment_2, typename K::Segment_2> {
  typedef typename boost::variant<typename K::Segment_2, typename K::Point_2 > variant_type;
  typedef typename boost::optional< variant_type > result_type;
};


template <class K, class OutputIterator>
OutputIterator intersect_do_iterator(const typename K::Segment_2 &seg1, 
                           const typename K::Segment_2 &seg2,
                           const K&, OutputIterator o) {
  typedef CGAL::internal::Segment_2_Segment_2_pair<K> is_t;

  is_t ispair(&seg1, &seg2);
  switch (ispair.intersection_type()) {
  case is_t::NO_INTERSECTION:
  default:
    return o;
  case is_t::POINT:
    *o++ = ispair.intersection_point();
    return o;
  case is_t::SEGMENT:
    *o++ = ispair.intersection_segment();
    return o;
  }
}


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

template<typename F>
void intersect_each(F f, const Vector& segs) {
  for(Vector::const_iterator it = segs.begin(); it != segs.end(); ++it) {
      const Segment& seg_1 = *it;
      for(Vector::const_iterator it2 = segs.begin(); it2 != segs.end(); ++it2) {
        f(seg_1, *it2);
      }
  }
}

struct Vec_holder {
  Vec_holder(std::vector<Point>* p, std::vector<Segment>* s) : p(p), s(s) { }
protected:
  std::vector<Point>* p;
  std::vector<Segment>* s;
};

struct Visitor : public boost::static_visitor<>, Vec_holder
{
  Visitor(std::vector<Point>* p, std::vector<Segment>* s) : 
    Vec_holder(p, s) { }

  void operator()(const Point& point) { p->push_back(point);  }  
  void operator()(const Segment& segment) { s->push_back(segment);  }
};

struct Variant_f {
  Variant_f(std::vector<Point>* p, std::vector<Segment>* s) : v(p, s)
    { }
  typedef Intersection_traits<K, Segment, Segment> Traits;
  typedef Traits::result_type result_type;

  Visitor v;
  void operator()(const Segment& s1, const Segment& s2) {
    result_type obj = intersection_variant(s1, s2, K());
    if(obj) {
      boost::apply_visitor(v, *obj);
    } 
  }
};

struct Object_f : Vec_holder {
  Object_f(std::vector<Point>* p, std::vector<Segment>* s) : 
    Vec_holder(p, s) { }
  
  void operator()(const Segment& s1, const Segment& s2) {
    Object obj = intersection(s1, s2);
      if (const Point * point = object_cast<Point>(&obj)) {
        p->push_back(*point);
      } else if (const Segment * segment = object_cast<Segment>(&obj)) {
        s->push_back(*segment);
      } 
  }
};

struct Any_f : Vec_holder {
  Any_f(std::vector<Point>* p, std::vector<Segment>* s) : 
    Vec_holder(p, s) { }
  
  void operator()(const Segment& s1, const Segment& s2) {
    boost::any obj = intersection_any(s1, s2, K());
    if (const Point * point = boost::any_cast<Point>(&obj)) {
       p->push_back(*point);
    } else if (const Segment * segment = boost::any_cast<Segment>(&obj)) {
      s->push_back(*segment);
    } 
  }
};  


struct Object_from_variant_f : Vec_holder {
  Object_from_variant_f(std::vector<Point>* p, std::vector<Segment>* s) : 
    Vec_holder(p, s) { }
  
  void operator()(const Segment& s1, const Segment& s2) {
    Object obj = intersection_variant(s1, s2, K());
      if (const Point * point = object_cast<Point>(&obj)) {
        p->push_back(*point);
      } else if (const Segment * segment = object_cast<Segment>(&obj)) {
        s->push_back(*segment);
      } 
  }
};

struct Do_f : Vec_holder {
  Do_f(std::vector<Point>* p, std::vector<Segment>* s) : 
    Vec_holder(p, s) { }

  typedef typename std::back_insert_iterator< std::vector<Point> >   Iter1;
  typedef typename std::back_insert_iterator< std::vector<Segment> > Iter2;

  void operator()(const Segment& s1, const Segment& s2) {

    CGAL::Dispatch_or_drop_output_iterator<CGAL::cpp0x::tuple<Point,Segment>,
                                           CGAL::cpp0x::tuple<Iter1,Iter2>
                                           > do_it(std::back_inserter(*p), std::back_inserter(*s));
    
    intersect_do_iterator(s1, s2, K(), do_it);
  }
};

#if !defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) && !defined(CGAL_CFG_NO_CPP0X_TUPLE)
cpp0x::tuple<int, int, int> intersect_each_variant_overload(const Vector& segs) {
  cpp0x::tuple<int, int, int> ret = cpp0x::make_tuple(0, 0, 0);
  typedef Intersection_traits<K, Segment, Segment> Traits;
  typedef Traits::result_type result_type;

  // Calculate the intersections between each segment
  for(Vector::const_iterator it = segs.begin(); it != segs.end(); ++it) {
    const Segment& seg_1 = *it;
    for(Vector::const_iterator it2 = segs.begin(); it2 != segs.end(); ++it2) {
      result_type obj = intersection_variant(seg_1, *it2, K());
      if(obj) {
 	// with c++0x
        auto v = make_overload(
          std::make_tuple(std::function<void(const Segment&)>(
                            [&ret](const Segment& s) { (void)s; ++(cpp0x::get<1>(ret)); }),
                          std::function<void(const Point&)>(
                            [&ret](const Point& p) { (void)p; ++(cpp0x::get<0>(ret)); })));

        boost::apply_visitor(v, *obj);
      } else {
	++(cpp0x::get<2>(ret));
      }
    }
  }
  
  return ret;
}
#endif

int main(int argc, char* argv[]) {
  int repeats = 100;
  int seg_count = 200;

  if(argc > 1)
    repeats = boost::lexical_cast<int>(argv[1]);
  if(argc > 2)
    seg_count = boost::lexical_cast<int>(argv[2]);

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
  CGAL::copy_n( g, seg_count, std::back_inserter(segs));

  std::vector<Point> points;
  std::vector<Segment> segments;
  points.clear(); segments.clear(); points.reserve(0); segments.reserve(0);
  //one run to get the size
  intersect_each(Variant_f(&points, &segments), segs);

  boost::timer timer;

  // variant vs object vs any

  points.clear(); segments.clear();
  timer.restart();

  for(int i = 0; i < repeats; ++i) {
    intersect_each(Object_f(&points, &segments), segs);
    points.clear(); segments.clear();
  }
  std::cout << "Time for object: " << timer.elapsed() << '\n';

  points.clear(); segments.clear();
  timer.restart();

  for(int i = 0; i < repeats; ++i) {
    intersect_each(Variant_f(&points, &segments), segs);
    points.clear(); segments.clear();
  }
  std::cout << "Time for variant: " << timer.elapsed() << '\n';

  points.clear(); segments.clear();
  timer.restart();

  for(int i = 0; i < repeats; ++i) {
    intersect_each(Any_f(&points, &segments), segs);
    points.clear(); segments.clear();
  }
  std::cout << "Time for any: " << timer.elapsed() << '\n';

  points.clear(); segments.clear();
  timer.restart();

  for(int i = 0; i < repeats; ++i) {
    intersect_each(Object_from_variant_f(&points, &segments), segs);
    points.clear(); segments.clear();
  }
  std::cout << "Time for object_from_variant: " << timer.elapsed() << '\n';

  points.clear(); segments.clear();
  timer.restart();

  for(int i = 0; i < repeats; ++i) {
    intersect_each(Do_f(&points, &segments), segs);
    points.clear(); segments.clear();
  }
  std::cout << "Time for dispatch_output: " << timer.elapsed() << '\n';

  std::cout << std::flush;
}
