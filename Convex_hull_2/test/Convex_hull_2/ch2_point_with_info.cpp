#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_traits_2.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;


template <class F, class Point>
struct Forward_bool_functor
  : public F
{
  const std::vector<Point>& points;

  Forward_bool_functor(const std::vector<Point>& points)
    : points(points)
  {}

  template <class Id>
  bool operator() (const Id& p, const Id& q) const
  {
    return static_cast<const F*>(this)->operator()(points[p], points[q]);
  }

  template <class Id>
  bool operator() (const Id& p, const Id& q, const Id& r) const
  {
    return static_cast<const F*>(this)->operator()(points[p], points[q], points[r]);
  }

  template <class Id>
  bool operator() (const Id& p, const Id& q, const Id& r, const Id& s) const
  {
    return static_cast<const F*>(this)->operator()(points[p], points[q], points[r], points[s]);
  }
};

template <class K>
struct CH_traits_for_point_ids
{
  const std::vector<typename K::Point_2>& points;

  CH_traits_for_point_ids(const std::vector<typename K::Point_2>& points)
    : points(points)
  {}

  typedef std::size_t Point_2;
  typedef CGAL::Convex_hull_traits_2<K> Base;
  typedef Forward_bool_functor<typename Base::Less_xy_2, typename K::Point_2> Less_xy_2;
  typedef Forward_bool_functor<typename Base::Less_yx_2, typename K::Point_2> Less_yx_2;
  typedef Forward_bool_functor<typename Base::Less_signed_distance_to_line_2, typename K::Point_2> Less_signed_distance_to_line_2;
  typedef Forward_bool_functor<typename Base::Less_rotate_ccw_2, typename K::Point_2> Less_rotate_ccw_2;
  typedef Forward_bool_functor<typename Base::Left_turn_2, typename K::Point_2> Left_turn_2;
  typedef Forward_bool_functor<typename Base::Equal_2, typename K::Point_2> Equal_2;

  struct Orientation_2
  {
    const std::vector<typename K::Point_2>& points;

    Orientation_2(const std::vector<typename K::Point_2>& points)
      : points(points)
    {}

    CGAL::Orientation
    operator()(Point_2 p, Point_2 q, Point_2 r) const
    {
      return typename Base::Orientation_2()(points[p], points[q], points[r]);
    }
  };

  Equal_2 equal_2_object () const
  {
    return Equal_2(points);
  }

  Less_xy_2 less_xy_2_object () const
  {
    return Less_xy_2(points);
  }

  Less_yx_2 less_yx_2_object () const
  {
    return Less_yx_2(points);
  }

  Less_signed_distance_to_line_2 less_signed_distance_to_line_2_object () const
  {
    return Less_signed_distance_to_line_2(points);
  }

  Less_rotate_ccw_2 less_rotate_ccw_2_object () const
  {
    return Less_rotate_ccw_2(points);
  }

  Left_turn_2 left_turn_2_object () const
  {
    return Left_turn_2(points);
  }

  Orientation_2 orientation_2_object () const
  {
    return Orientation_2(points);
  }
};

int main()
{
  std::vector<Point_2> input_points;
  std::vector<std::size_t> result;

  input_points.push_back( Point_2(0,0) );
  input_points.push_back( Point_2(0,1) );
  input_points.push_back( Point_2(1,0) );
  input_points.push_back( Point_2(0.25,0.25) );

  CGAL::convex_hull_2( boost::counting_iterator<std::size_t>(0),
                       boost::counting_iterator<std::size_t>(input_points.size()),
                       std::back_inserter(result), CH_traits_for_point_ids<K>(input_points) );

  for(std::size_t i : result)
  {
    std::cout << input_points[i] << " - " << i << "\n";
  }

  assert( result.size() == 3 );

  return 0;
}
