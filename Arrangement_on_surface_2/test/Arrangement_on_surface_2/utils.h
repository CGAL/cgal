
#include <iostream>
#include <CGAL/enum.h>

// TODO: Add constructors that acceps a pointer to the traits object

template <typename Traits>
class Point_compare
{
  typedef typename Traits::Point_2      Point_2;

public:
  bool operator()(const Point_2& p1, const Point_2& p2)
  {
    Traits tr;
    CGAL::Comparison_result res = tr.compare_xy_2_object()(p1, p2);
    return (res == CGAL::SMALLER) ? true : false;
  }
};


template <typename Traits>
class Curve_compare
{
  typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Traits::Point_2               Point_2;

public:
  bool operator()(const X_monotone_curve_2& c1, const X_monotone_curve_2& c2)
  {
    Traits tr;
    const Point_2& c1_left = tr.construct_min_vertex_2_object()(c1);
    const Point_2& c2_left = tr.construct_min_vertex_2_object()(c2);

    CGAL::Comparison_result res = tr.compare_xy_2_object()(c1_left, c2_left);
    
    if(res == CGAL::SMALLER)
      return true;

    if(res == CGAL::LARGER)
      return false;

    CGAL_assertion(res == CGAL::EQUAL);
    res = tr.compare_y_at_x_right_2_object()(c1, c2, c1_left);
    return (res == CGAL::SMALLER) ? true : false;
  }
};
