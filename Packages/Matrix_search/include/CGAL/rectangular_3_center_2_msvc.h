// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : rectangular_3_center_2_msvc.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : 3cover.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// MSVC Workarounds
// ============================================================================

#if ! (CGAL_RECTANGULAR_3_CENTER_2_MSVC_H)
#define CGAL_RECTANGULAR_3_CENTER_2_MSVC_H 1

template < class R >
struct Rectangular_3_center_2_type2_operations0 {
  typedef typename R::FT                         FT;
  typedef typename R::Point_2                    Point_2;
  typedef typename R::Iso_rectangle_2            Iso_rectangle_2;
  typedef typename R::Infinity_distance_2        Infinity_distance_2;
  typedef typename R::Less_x_2                   Less_x_2;
  typedef typename R::Less_y_2                   Less_y_2;
  typedef typename R::Greater_x_2                Greater_x_2;
  typedef typename R::Greater_y_2                Greater_y_2;
  typedef typename R::Construct_min_2            Construct_min_2;
  typedef typename R::Construct_max_2            Construct_max_2;
  typedef typename R::Construct_corner_2         Construct_corner_2;
  typedef typename R::Construct_iso_rectangle_2  C_rect;
  typedef typename R::Construct_projection_onto_horizontal_implicit_line_2
    Construct_projection_onto_horizontal_implicit_line_2;
  typedef typename R::Construct_iso_rectangle_2_below_left_point_2
    Construct_iso_rectangle_2_below_left_point_2;
  typedef typename R::Construct_iso_rectangle_2_above_left_point_2
    Construct_iso_rectangle_2_above_left_point_2;
  typedef typename R::Construct_iso_rectangle_2_below_right_point_2
    Construct_iso_rectangle_2_below_right_point_2;
  typedef typename R::Construct_iso_rectangle_2_above_right_point_2
    Construct_iso_rectangle_2_above_right_point_2;
  typedef std::binder1st< Infinity_distance_2 >   Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x() const { return r_.get_less_x_2(); }
  Less_y_2  less_y() const { return r_.get_less_y_2(); }
  Greater_x_2  greater_x() const { return r_.get_greater_x_2(); }
  Greater_y_2  greater_y() const { return r_.get_greater_y_2(); }
  Infinity_distance_2  distance() const
  { return r_.get_infinity_distance_2(); }
  Construct_min_2  rmin() const { return r_.get_construct_min_2(); }
  Construct_max_2  rmax() const { return r_.get_construct_max_2(); }
  Construct_corner_2  corner() const { return r_.get_construct_corner_2(); }
  C_rect rectangle() const { return r_.get_construct_iso_rectangle_2(); }
  
  Construct_iso_rectangle_2_below_left_point_2
  rect_b_l() const
  { return r_.get_construct_iso_rectangle_2_below_left_point_2(); }
  Construct_iso_rectangle_2_above_left_point_2
  rect_a_l() const
  { return r_.get_construct_iso_rectangle_2_above_left_point_2(); }
  Construct_iso_rectangle_2_below_right_point_2
  rect_b_r() const
  { return r_.get_construct_iso_rectangle_2_below_right_point_2(); }
  Construct_iso_rectangle_2_above_right_point_2
  rect_a_r() const
  { return r_.get_construct_iso_rectangle_2_above_right_point_2(); }
  
  typedef typename R::Min_x_2 Min_x_2;
  Min_x_2 minx() const { return r_.get_min_x_2(); }
  typedef typename R::Min_y_2 Min_y_2;
  Min_y_2 miny() const { return r_.get_min_y_2(); }
  typedef typename R::Max_x_2 Max_x_2;
  Max_x_2 maxx() const { return r_.get_max_x_2(); }
  typedef typename R::Max_y_2 Max_y_2;
  Max_y_2 maxy() const { return r_.get_max_y_2(); }
  
  Construct_projection_onto_horizontal_implicit_line_2
  build_point() const
  { return r_.get_construct_projection_onto_horizontal_implicit_line_2(); }
  
  private:
    R& r_;
    Delta delta_;
public:
  typedef Point_2              Point;
  typedef Iso_rectangle_2      Rectangle;
  typedef Less_x_2             X_compare;
  typedef Greater_y_2          Y_compare;
  typedef Infinity_distance_2  Distance;

  Rectangular_3_center_2_type2_operations0(R& r, const Point& p)
  : r_(r), delta_(std::bind1st(r.get_infinity_distance_2(), p))
  {}

  X_compare  compare_x() const { return less_x(); }
  Y_compare  compare_y() const { return greater_y(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    return build_point()(
      constraint_empty ? first_uncovered :
      minx()(first_uncovered, rmin()(constraint)),
        rmax()(bbox));
  }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    return build_point()(
      constraint_empty ? rmax()(bbox) : rmin()(constraint),
      rmax()(bbox));
  }
  
  Point place_x_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return build_point()(
      minx()(
        rmin()(rect_b_l()(rmax()(bbox), radius)), so_far),
        so_far);
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    return build_point()(
      rmax()(bbox),
      constraint_empty ? first_uncovered :
    miny()(first_uncovered, rmin()(constraint)));
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    return build_point()(
      rmax()(bbox),
      constraint_empty ? rmax()(bbox) : rmin()(constraint));
  }
  
  Point place_y_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return build_point()(
      so_far,
      miny()(rmin()(rect_b_l()(rmax()(bbox), radius)),
             so_far));
  }
  
  Point update_x_square(const Point& s, const Point& newp) const
  { return build_point()(minx()(s, newp), s); }
  
  Point update_y_square(const Point& s, const Point& newp) const
  { return build_point()(s, miny()(s, newp)); }
  
  FT compute_x_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 1)); }
  
  FT compute_y_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 3)); }
  
  Point construct_corner_square(const Rectangle& bbox, FT r) const
  { return rmax()(rect_a_r()(rmin()(bbox), r)); }
  
  Point construct_x_square(const Point& p, FT r) const
  { return corner()(rect_b_r()(p, r), 1); }
  
  Point construct_y_square(const Point& p, FT r) const
  { return corner()(rect_a_l()(p, r), 3); }
};
template < class R >
struct Rectangular_3_center_2_type2_operations1 {
  typedef typename R::FT                         FT;
  typedef typename R::Point_2                    Point_2;
  typedef typename R::Iso_rectangle_2            Iso_rectangle_2;
  typedef typename R::Infinity_distance_2        Infinity_distance_2;
  typedef typename R::Less_x_2                   Less_x_2;
  typedef typename R::Less_y_2                   Less_y_2;
  typedef typename R::Greater_x_2                Greater_x_2;
  typedef typename R::Greater_y_2                Greater_y_2;
  typedef typename R::Construct_min_2            Construct_min_2;
  typedef typename R::Construct_max_2            Construct_max_2;
  typedef typename R::Construct_corner_2         Construct_corner_2;
  typedef typename R::Construct_iso_rectangle_2  C_rect;
  typedef typename R::Construct_projection_onto_horizontal_implicit_line_2
    Construct_projection_onto_horizontal_implicit_line_2;
  typedef typename R::Construct_iso_rectangle_2_below_left_point_2
    Construct_iso_rectangle_2_below_left_point_2;
  typedef typename R::Construct_iso_rectangle_2_above_left_point_2
    Construct_iso_rectangle_2_above_left_point_2;
  typedef typename R::Construct_iso_rectangle_2_below_right_point_2
    Construct_iso_rectangle_2_below_right_point_2;
  typedef typename R::Construct_iso_rectangle_2_above_right_point_2
    Construct_iso_rectangle_2_above_right_point_2;
  typedef std::binder1st< Infinity_distance_2 >   Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x() const { return r_.get_less_x_2(); }
  Less_y_2  less_y() const { return r_.get_less_y_2(); }
  Greater_x_2  greater_x() const { return r_.get_greater_x_2(); }
  Greater_y_2  greater_y() const { return r_.get_greater_y_2(); }
  Infinity_distance_2  distance() const
  { return r_.get_infinity_distance_2(); }
  Construct_min_2  rmin() const { return r_.get_construct_min_2(); }
  Construct_max_2  rmax() const { return r_.get_construct_max_2(); }
  Construct_corner_2  corner() const { return r_.get_construct_corner_2(); }
  C_rect rectangle() const { return r_.get_construct_iso_rectangle_2(); }
  
  Construct_iso_rectangle_2_below_left_point_2
  rect_b_l() const
  { return r_.get_construct_iso_rectangle_2_below_left_point_2(); }
  Construct_iso_rectangle_2_above_left_point_2
  rect_a_l() const
  { return r_.get_construct_iso_rectangle_2_above_left_point_2(); }
  Construct_iso_rectangle_2_below_right_point_2
  rect_b_r() const
  { return r_.get_construct_iso_rectangle_2_below_right_point_2(); }
  Construct_iso_rectangle_2_above_right_point_2
  rect_a_r() const
  { return r_.get_construct_iso_rectangle_2_above_right_point_2(); }
  
  typedef typename R::Min_x_2 Min_x_2;
  Min_x_2 minx() const { return r_.get_min_x_2(); }
  typedef typename R::Min_y_2 Min_y_2;
  Min_y_2 miny() const { return r_.get_min_y_2(); }
  typedef typename R::Max_x_2 Max_x_2;
  Max_x_2 maxx() const { return r_.get_max_x_2(); }
  typedef typename R::Max_y_2 Max_y_2;
  Max_y_2 maxy() const { return r_.get_max_y_2(); }
  
  Construct_projection_onto_horizontal_implicit_line_2
  build_point() const
  { return r_.get_construct_projection_onto_horizontal_implicit_line_2(); }
  
  private:
    R& r_;
    Delta delta_;
public:
  typedef Point_2              Point;
  typedef Iso_rectangle_2      Rectangle;
  typedef Greater_x_2          X_compare;
  typedef Greater_y_2          Y_compare;
  typedef Infinity_distance_2  Distance;

  Rectangular_3_center_2_type2_operations1(R& r, const Point& p)
  : r_(r), delta_(std::bind1st(r.get_infinity_distance_2(), p))
  {}

  X_compare  compare_x() const { return greater_x(); }
  Y_compare  compare_y() const { return greater_y(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    return build_point()(
      constraint_empty ? first_uncovered :
      maxx()(first_uncovered, rmax()(constraint)),
        rmax()(bbox));
  }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    return build_point()(
      constraint_empty ? rmin()(bbox) : rmax()(constraint),
      rmax()(bbox));
  }
  
  Point place_x_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return build_point()(
      maxx()(so_far, rmax()(rect_a_r()(rmin()(bbox), radius))),
      so_far);
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    return build_point()(
      rmin()(bbox),
      constraint_empty ? first_uncovered :
    miny()(first_uncovered, rmin()(constraint)));
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    return build_point()(
      rmin()(bbox),
      constraint_empty ? rmax()(bbox) : rmin()(constraint));
  }
  
  Point place_y_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return build_point()(
      so_far,
      miny()(so_far, rmin()(rect_b_l()(rmax()(bbox), radius))));
  }
  
  Point update_x_square(const Point& s, const Point& newp) const
  { return build_point()(maxx()(s, newp), s); }
  
  Point update_y_square(const Point& s, const Point& newp) const
  { return build_point()(s, miny()(s, newp)); }
  
  FT compute_x_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 0)); }
  
  FT compute_y_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 2)); }
  
  Point construct_corner_square(const Rectangle& bbox, FT r) const
  { return corner()(rect_a_l()(corner()(bbox, 1), r), 3); }
  
  Point construct_x_square(const Point& p, FT r) const
  { return rmin()(rect_b_l()(p, r)); }
  
  Point construct_y_square(const Point& p, FT r) const
  { return rmax()(rect_a_r()(p, r)); }
};
template < class R >
struct Rectangular_3_center_2_type2_operations2 {
  typedef typename R::FT                         FT;
  typedef typename R::Point_2                    Point_2;
  typedef typename R::Iso_rectangle_2            Iso_rectangle_2;
  typedef typename R::Infinity_distance_2        Infinity_distance_2;
  typedef typename R::Less_x_2                   Less_x_2;
  typedef typename R::Less_y_2                   Less_y_2;
  typedef typename R::Greater_x_2                Greater_x_2;
  typedef typename R::Greater_y_2                Greater_y_2;
  typedef typename R::Construct_min_2            Construct_min_2;
  typedef typename R::Construct_max_2            Construct_max_2;
  typedef typename R::Construct_corner_2         Construct_corner_2;
  typedef typename R::Construct_iso_rectangle_2  C_rect;
  typedef typename R::Construct_projection_onto_horizontal_implicit_line_2
    Construct_projection_onto_horizontal_implicit_line_2;
  typedef typename R::Construct_iso_rectangle_2_below_left_point_2
    Construct_iso_rectangle_2_below_left_point_2;
  typedef typename R::Construct_iso_rectangle_2_above_left_point_2
    Construct_iso_rectangle_2_above_left_point_2;
  typedef typename R::Construct_iso_rectangle_2_below_right_point_2
    Construct_iso_rectangle_2_below_right_point_2;
  typedef typename R::Construct_iso_rectangle_2_above_right_point_2
    Construct_iso_rectangle_2_above_right_point_2;
  typedef std::binder1st< Infinity_distance_2 >   Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x() const { return r_.get_less_x_2(); }
  Less_y_2  less_y() const { return r_.get_less_y_2(); }
  Greater_x_2  greater_x() const { return r_.get_greater_x_2(); }
  Greater_y_2  greater_y() const { return r_.get_greater_y_2(); }
  Infinity_distance_2  distance() const
  { return r_.get_infinity_distance_2(); }
  Construct_min_2  rmin() const { return r_.get_construct_min_2(); }
  Construct_max_2  rmax() const { return r_.get_construct_max_2(); }
  Construct_corner_2  corner() const { return r_.get_construct_corner_2(); }
  C_rect rectangle() const { return r_.get_construct_iso_rectangle_2(); }
  
  Construct_iso_rectangle_2_below_left_point_2
  rect_b_l() const
  { return r_.get_construct_iso_rectangle_2_below_left_point_2(); }
  Construct_iso_rectangle_2_above_left_point_2
  rect_a_l() const
  { return r_.get_construct_iso_rectangle_2_above_left_point_2(); }
  Construct_iso_rectangle_2_below_right_point_2
  rect_b_r() const
  { return r_.get_construct_iso_rectangle_2_below_right_point_2(); }
  Construct_iso_rectangle_2_above_right_point_2
  rect_a_r() const
  { return r_.get_construct_iso_rectangle_2_above_right_point_2(); }
  
  typedef typename R::Min_x_2 Min_x_2;
  Min_x_2 minx() const { return r_.get_min_x_2(); }
  typedef typename R::Min_y_2 Min_y_2;
  Min_y_2 miny() const { return r_.get_min_y_2(); }
  typedef typename R::Max_x_2 Max_x_2;
  Max_x_2 maxx() const { return r_.get_max_x_2(); }
  typedef typename R::Max_y_2 Max_y_2;
  Max_y_2 maxy() const { return r_.get_max_y_2(); }
  
  Construct_projection_onto_horizontal_implicit_line_2
  build_point() const
  { return r_.get_construct_projection_onto_horizontal_implicit_line_2(); }
  
  private:
    R& r_;
    Delta delta_;
public:
  typedef Point_2              Point;
  typedef Iso_rectangle_2      Rectangle;
  typedef Greater_x_2          X_compare;
  typedef Less_y_2             Y_compare;
  typedef Infinity_distance_2  Distance;

  Rectangular_3_center_2_type2_operations2(R& r, const Point& p)
  : r_(r), delta_(std::bind1st(r.get_infinity_distance_2(), p))
  {}

  X_compare  compare_x() const { return greater_x(); }
  Y_compare  compare_y() const { return less_y(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    return build_point()(
      constraint_empty ? first_uncovered :
      maxx()(first_uncovered, rmax()(constraint)),
        rmin()(bbox));
  }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    return build_point()(
      constraint_empty ? rmin()(bbox) : rmax()(constraint),
      rmin()(bbox));
  }
  
  Point place_x_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return build_point()(
      maxx()(so_far, rmax()(rect_a_r()(rmin()(bbox), radius))),
      so_far);
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    return build_point()(
      rmin()(bbox),
      constraint_empty ? first_uncovered :
    maxy()(first_uncovered, rmax()(constraint)));
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    return build_point()(
      rmin()(bbox),
      constraint_empty ? rmin()(bbox) : rmax()(constraint));
  }
  
  Point place_y_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return build_point()(
      so_far,
      maxy()(rmax()(rect_a_r()(rmin()(bbox), radius)), so_far));
  }
  
  Point update_x_square(const Point& s, const Point& newp) const
  { return build_point()(maxx()(s, newp), s); }
  
  Point update_y_square(const Point& s, const Point& newp) const
  { return build_point()(s, maxy()(s, newp)); }
  
  FT compute_x_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 3)); }
  
  FT compute_y_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 1)); }
  
  Point construct_corner_square(const Rectangle& bbox, FT r) const
  { return rmin()(rect_b_l()(rmax()(bbox), r)); }
  
  Point construct_x_square(const Point& p, FT r) const
  { return corner()(rect_a_l()(p, r), 3); }
  
  Point construct_y_square(const Point& p, FT r) const
  { return corner()(rect_b_r()(p, r), 1); }
};
template < class R >
struct Rectangular_3_center_2_type2_operations3 {
  typedef typename R::FT                         FT;
  typedef typename R::Point_2                    Point_2;
  typedef typename R::Iso_rectangle_2            Iso_rectangle_2;
  typedef typename R::Infinity_distance_2        Infinity_distance_2;
  typedef typename R::Less_x_2                   Less_x_2;
  typedef typename R::Less_y_2                   Less_y_2;
  typedef typename R::Greater_x_2                Greater_x_2;
  typedef typename R::Greater_y_2                Greater_y_2;
  typedef typename R::Construct_min_2            Construct_min_2;
  typedef typename R::Construct_max_2            Construct_max_2;
  typedef typename R::Construct_corner_2         Construct_corner_2;
  typedef typename R::Construct_iso_rectangle_2  C_rect;
  typedef typename R::Construct_projection_onto_horizontal_implicit_line_2
    Construct_projection_onto_horizontal_implicit_line_2;
  typedef typename R::Construct_iso_rectangle_2_below_left_point_2
    Construct_iso_rectangle_2_below_left_point_2;
  typedef typename R::Construct_iso_rectangle_2_above_left_point_2
    Construct_iso_rectangle_2_above_left_point_2;
  typedef typename R::Construct_iso_rectangle_2_below_right_point_2
    Construct_iso_rectangle_2_below_right_point_2;
  typedef typename R::Construct_iso_rectangle_2_above_right_point_2
    Construct_iso_rectangle_2_above_right_point_2;
  typedef std::binder1st< Infinity_distance_2 >   Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x() const { return r_.get_less_x_2(); }
  Less_y_2  less_y() const { return r_.get_less_y_2(); }
  Greater_x_2  greater_x() const { return r_.get_greater_x_2(); }
  Greater_y_2  greater_y() const { return r_.get_greater_y_2(); }
  Infinity_distance_2  distance() const
  { return r_.get_infinity_distance_2(); }
  Construct_min_2  rmin() const { return r_.get_construct_min_2(); }
  Construct_max_2  rmax() const { return r_.get_construct_max_2(); }
  Construct_corner_2  corner() const { return r_.get_construct_corner_2(); }
  C_rect rectangle() const { return r_.get_construct_iso_rectangle_2(); }
  
  Construct_iso_rectangle_2_below_left_point_2
  rect_b_l() const
  { return r_.get_construct_iso_rectangle_2_below_left_point_2(); }
  Construct_iso_rectangle_2_above_left_point_2
  rect_a_l() const
  { return r_.get_construct_iso_rectangle_2_above_left_point_2(); }
  Construct_iso_rectangle_2_below_right_point_2
  rect_b_r() const
  { return r_.get_construct_iso_rectangle_2_below_right_point_2(); }
  Construct_iso_rectangle_2_above_right_point_2
  rect_a_r() const
  { return r_.get_construct_iso_rectangle_2_above_right_point_2(); }
  
  typedef typename R::Min_x_2 Min_x_2;
  Min_x_2 minx() const { return r_.get_min_x_2(); }
  typedef typename R::Min_y_2 Min_y_2;
  Min_y_2 miny() const { return r_.get_min_y_2(); }
  typedef typename R::Max_x_2 Max_x_2;
  Max_x_2 maxx() const { return r_.get_max_x_2(); }
  typedef typename R::Max_y_2 Max_y_2;
  Max_y_2 maxy() const { return r_.get_max_y_2(); }
  
  Construct_projection_onto_horizontal_implicit_line_2
  build_point() const
  { return r_.get_construct_projection_onto_horizontal_implicit_line_2(); }
  
  private:
    R& r_;
    Delta delta_;
public:
  typedef Point_2              Point;
  typedef Iso_rectangle_2      Rectangle;
  typedef Less_x_2             X_compare;
  typedef Less_y_2             Y_compare;
  typedef Infinity_distance_2  Distance;

  Rectangular_3_center_2_type2_operations3(R& r, const Point& p)
  : r_(r), delta_(std::bind1st(r.get_infinity_distance_2(), p))
  {}

  X_compare  compare_x() const { return less_x(); }
  Y_compare  compare_y() const { return less_y(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    return build_point()(
      constraint_empty ? first_uncovered :
      minx()(first_uncovered, rmin()(constraint)),
        rmin()(bbox));
  }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    return build_point()(
      constraint_empty ? rmax()(bbox) : rmin()(constraint),
      rmin()(bbox));
  }
  
  Point place_x_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return build_point()(
      minx()(rmin()(rect_b_l()(rmax()(bbox), radius)), so_far),
      so_far);
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    return build_point()(
      rmax()(bbox),
      constraint_empty ? first_uncovered :
    maxy()(first_uncovered, rmax()(constraint)));
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    return build_point()(
      rmax()(bbox),
      constraint_empty ? rmin()(bbox) : rmax()(constraint));
  }
  
  Point place_y_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return build_point()(
      so_far,
      maxy()(rmax()(rect_a_r()(rmin()(bbox), radius)), so_far));
  }
  
  Point update_x_square(const Point& s, const Point& newp) const
  { return build_point()(minx()(s, newp), s); }
  
  Point update_y_square(const Point& s, const Point& newp) const
  { return build_point()(s, maxy()(s, newp)); }
  
  FT compute_x_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 2)); }
  
  FT compute_y_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 0)); }
  
  Point construct_corner_square(const Rectangle& bbox, FT r) const
  { return corner()(rect_b_r()(corner()(bbox, 3), r), 1); }
  
  Point construct_x_square(const Point& p, FT r) const
  { return rmax()(rect_a_r()(p, r)); }
  
  Point construct_y_square(const Point& p, FT r) const
  { return rmin()(rect_b_l()(p, r)); }
};

#endif // ! (CGAL_RECTANGULAR_3_CENTER_2_MSVC_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

