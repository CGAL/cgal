#include <iostream>

#include <CGAL/enum.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>

template <typename T_Geom_traits>
class Point_compare {
private:
  typedef T_Geom_traits                         Traits;
  typedef typename Traits::Left_side_category   Left_side_category;
  typedef typename Traits::Bottom_side_category Bottom_side_category;
  typedef typename Traits::Top_side_category    Top_side_category;
  typedef typename Traits::Right_side_category  Right_side_category;
  typedef typename
  CGAL::Arr_are_all_sides_oblivious_tag<Left_side_category,
                                        Bottom_side_category,
                                        Top_side_category,
                                        Right_side_category>::result
    Are_all_sides_oblivious_category;

  typedef typename CGAL::Arr_has_identified_sides<Left_side_category,
                                                  Bottom_side_category>::result
    Has_identified_sides_category;

  const Traits& m_traits;

public:
  typedef typename Traits::Point_2              Point_2;

  Point_compare(const Traits& traits) : m_traits(traits) {}

  bool operator()(const Point_2& p1, const Point_2& p2) const
  { return operator()(p1, p2, Are_all_sides_oblivious_category()); }

private:
  // The following set of operators is incomplete, but for now there are
  // no tests that require the missing ones. In particular, an operator
  // an operator for traits classes, where at least one boundary is either
  // identified or contracted and another boundary is open, is missing.

  // This function is invoked for traits classes where all boundaries
  // are oblivious.
  bool operator()(const Point_2& p1, const Point_2& p2,
                  CGAL::Arr_all_sides_oblivious_tag) const
  {
    return (m_traits.compare_xy_2_object()(p1, p2) == CGAL::SMALLER);
  }

  bool operator()(const Point_2& p1, const Point_2& p2,
                  CGAL::Arr_not_all_sides_oblivious_tag) const
  { return operator()(p1, p2, Has_identified_sides_category()); }

  // This function is invoked for traits classes where at least one
  // boundary is not oblivious and all boundaries are not identified.
  bool operator()(const Point_2& p1, const Point_2& p2,
                  boost::mpl::bool_<false>) const
  {
    return (m_traits.compare_xy_2_object()(p1, p2) == CGAL::SMALLER);
  }

  // This function should be invoked for traits classes where at least one
  // boundary is identified.
  bool operator()(const Point_2& p1, const Point_2& p2,
                  boost::mpl::bool_<true>) const
  {
    // Compare in y boundaries:
    CGAL::Arr_parameter_space ps_y1 =
      m_traits.parameter_space_in_y_2_object()(p1);
    CGAL::Arr_parameter_space ps_y2 =
      m_traits.parameter_space_in_y_2_object()(p2);
    if (ps_y1 == CGAL::ARR_BOTTOM_BOUNDARY)
      return (ps_y2 != CGAL::ARR_BOTTOM_BOUNDARY);
    if (ps_y1 == CGAL::ARR_TOP_BOUNDARY) return false;
    if (ps_y2 == CGAL::ARR_BOTTOM_BOUNDARY) return false;
    if (ps_y2 == CGAL::ARR_TOP_BOUNDARY) return true;

    // Compare in x boundaries:
    CGAL::Arr_parameter_space ps_x1 =
      m_traits.parameter_space_in_x_2_object()(p1);
    CGAL::Arr_parameter_space ps_x2 =
      m_traits.parameter_space_in_x_2_object()(p2);
    if (ps_x1 == CGAL::ARR_LEFT_BOUNDARY) {
      if (ps_x2 == CGAL::ARR_LEFT_BOUNDARY)
        return
          (m_traits.compare_y_on_boundary_2_object()(p1, p2) == CGAL::SMALLER);
      return true;
    }

    if (ps_x1 == CGAL::ARR_RIGHT_BOUNDARY) {
      if (ps_x2 == CGAL::ARR_RIGHT_BOUNDARY)
        return
          (m_traits.compare_y_on_boundary_2_object()(p1, p2) == CGAL::SMALLER);
      return false;
    }
    if (ps_x2 == CGAL::ARR_LEFT_BOUNDARY) return false;
    if (ps_x2 == CGAL::ARR_RIGHT_BOUNDARY) return true;

    // Compare interiors:
    return (m_traits.compare_xy_2_object()(p1, p2) == CGAL::SMALLER);
  }
};

template <typename T_Geom_traits>
class Curve_compare {
private:
  typedef T_Geom_traits      Traits;

  const Traits& m_traits;

public:
  typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Traits::Point_2               Point_2;

  Curve_compare(const Traits& traits) : m_traits(traits) {}

  bool operator()(const X_monotone_curve_2& c1, const X_monotone_curve_2& c2)
  {
    const Point_2& c1_left = m_traits.construct_min_vertex_2_object()(c1);
    const Point_2& c2_left = m_traits.construct_min_vertex_2_object()(c2);

    CGAL::Comparison_result res =
      m_traits.compare_xy_2_object()(c1_left, c2_left);

    if (res == CGAL::SMALLER) return true;
    if (res == CGAL::LARGER) return false;
    CGAL_assertion(res == CGAL::EQUAL);
    res = m_traits.compare_y_at_x_right_2_object()(c1, c2, c1_left);
    return (res == CGAL::SMALLER) ? true : false;
  }
};
