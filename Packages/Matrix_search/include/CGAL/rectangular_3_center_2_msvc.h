// ============================================================================
//
// Copyright (c) 1998, 1999, 2000 The CGAL Consortium
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
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
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
  typedef typename R::Construct_min_point_2      Construct_min_point_2;
  typedef typename R::Construct_max_point_2      Construct_max_point_2;
  typedef typename R::Construct_corner_2         Construct_corner_2;
  typedef typename R::Construct_iso_rectangle_2  Construct_iso_rectangle_2;
  typedef typename R::Construct_projection_onto_horizontal_implicit_line_2
    Construct_projection_onto_horizontal_implicit_line_2;
  typedef typename R::Construct_point_2_above_right_implicit_point_2
    Construct_point_2_above_right_implicit_point_2;
  typedef typename R::Construct_point_2_above_left_implicit_point_2
    Construct_point_2_above_left_implicit_point_2;
  typedef typename R::Construct_point_2_below_right_implicit_point_2
    Construct_point_2_below_right_implicit_point_2;
  typedef typename R::Construct_point_2_below_left_implicit_point_2
    Construct_point_2_below_left_implicit_point_2;
  typedef std::binder1st< Infinity_distance_2 >   Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x_2_object() const { return r_.less_x_2_object(); }
  Less_y_2  less_y_2_object() const { return r_.less_y_2_object(); }
  Greater_x_2  greater_x_2_object() const { return r_.greater_x_2_object(); }
  Greater_y_2  greater_y_2_object() const { return r_.greater_y_2_object(); }
  Infinity_distance_2  distance() const
  { return r_.infinity_distance_2_object(); }
  Construct_min_point_2  construct_min_point_2_object() const
  { return r_.construct_min_point_2_object(); }
  Construct_max_point_2  construct_max_point_2_object() const
  { return r_.construct_max_point_2_object(); }
  Construct_corner_2  corner() const
  { return r_.construct_corner_2_object(); }
  Construct_iso_rectangle_2 construct_iso_rectangle_2_object() const
  { return r_.construct_iso_rectangle_2_object(); }
  
  Construct_point_2_below_left_implicit_point_2
  pt_b_l() const
  { return r_.construct_point_2_below_left_implicit_point_2_object(); }
  Construct_point_2_above_left_implicit_point_2
  pt_a_l() const
  { return r_.construct_point_2_above_left_implicit_point_2_object(); }
  Construct_point_2_below_right_implicit_point_2
  pt_b_r() const
  { return r_.construct_point_2_below_right_implicit_point_2_object(); }
  Construct_point_2_above_right_implicit_point_2
  pt_a_r() const
  { return r_.construct_point_2_above_right_implicit_point_2_object(); }
  
  typedef typename R::Min_x_2 Min_x_2;
  Min_x_2 minx() const { return r_.min_x_2_object(); }
  typedef typename R::Min_y_2 Min_y_2;
  Min_y_2 miny() const { return r_.min_y_2_object(); }
  typedef typename R::Max_x_2 Max_x_2;
  Max_x_2 maxx() const { return r_.max_x_2_object(); }
  typedef typename R::Max_y_2 Max_y_2;
  Max_y_2 maxy() const { return r_.max_y_2_object(); }
  
  Construct_projection_onto_horizontal_implicit_line_2
  construct_projection_onto_horizontal_implicit_line_2_object() const
  { return r_.construct_projection_onto_horizontal_implicit_line_2_object(); }
  
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
  : r_(r), delta_(std::bind1st(r.infinity_distance_2_object(), p))
  {}

  X_compare  compare_x() const { return less_x_2_object(); }
  Y_compare  compare_y() const { return greater_y_2_object(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
      Point bpt = constraint_empty ? first_uncovered :
      minx()(first_uncovered, construct_min_point_2_object()(constraint));
      return construct_projection_onto_horizontal_implicit_line_2_object()(
        bpt, construct_max_point_2_object()(bbox));
  #else
      return construct_projection_onto_horizontal_implicit_line_2_object()(
        constraint_empty ? first_uncovered :
        minx()(first_uncovered, construct_min_point_2_object()(constraint)),
          construct_max_point_2_object()(bbox));
  #endif
    }
  
    Point place_x_square(bool constraint_empty,
                         const Rectangle& constraint,
                         const Rectangle& bbox) const
    {
  #ifdef __BORLANDC__
      Point bpt = constraint_empty ? construct_max_point_2_object()(bbox) :
      construct_min_point_2_object()(constraint);
      return construct_projection_onto_horizontal_implicit_line_2_object()(
        bpt, construct_max_point_2_object()(bbox));
  #else
      return construct_projection_onto_horizontal_implicit_line_2_object()(
        constraint_empty ?
        construct_max_point_2_object()(bbox) :
        construct_min_point_2_object()(constraint),
        construct_max_point_2_object()(bbox));
  #endif
    }
  
    Point place_x_square(const Point& so_far,
                         const Rectangle& bbox,
                         FT radius) const
    {
      return construct_projection_onto_horizontal_implicit_line_2_object()(
        minx()(
          pt_b_l()(construct_max_point_2_object()(bbox),
                   construct_max_point_2_object()(bbox), radius), so_far),
          so_far);
    }
  
    Point place_y_square(bool constraint_empty,
                         const Rectangle& constraint,
                         const Point& first_uncovered,
                         const Rectangle& bbox) const
    {
  #ifdef __BORLANDC__
      Point bpt = constraint_empty ? first_uncovered :
      miny()(first_uncovered, construct_min_point_2_object()(constraint));
      return construct_projection_onto_horizontal_implicit_line_2_object()(
        construct_max_point_2_object()(bbox), bpt);
  #else
      return construct_projection_onto_horizontal_implicit_line_2_object()(
        construct_max_point_2_object()(bbox),
        constraint_empty ? first_uncovered :
      miny()(first_uncovered, construct_min_point_2_object()(constraint)));
  #endif
    }
  
    Point place_y_square(bool constraint_empty,
                         const Rectangle& constraint,
                         const Rectangle& bbox) const
    {
  #ifdef __BORLANDC__
      Point bpt = constraint_empty ? construct_max_point_2_object()(bbox) :
      construct_min_point_2_object()(constraint);
      return construct_projection_onto_horizontal_implicit_line_2_object()(
        construct_max_point_2_object()(bbox), bpt);
  #else
      return construct_projection_onto_horizontal_implicit_line_2_object()(
        construct_max_point_2_object()(bbox),
        constraint_empty ?
      construct_max_point_2_object()(bbox) :
      construct_min_point_2_object()(constraint));
  #endif
    }
  
    Point place_y_square(const Point& so_far,
                         const Rectangle& bbox,
                         FT radius) const
    {
      return construct_projection_onto_horizontal_implicit_line_2_object()(
        so_far,
        miny()(pt_b_l()(construct_max_point_2_object()(bbox),
                        construct_max_point_2_object()(bbox), radius),
               so_far));
    }
  
    Point update_x_square(const Point& s, const Point& newp) const
    { return construct_projection_onto_horizontal_implicit_line_2_object()(
      minx()(s, newp), s); }
  
    Point update_y_square(const Point& s, const Point& newp) const
    { return construct_projection_onto_horizontal_implicit_line_2_object()(
      s, miny()(s, newp)); }
  
    FT compute_x_distance(const Point& extreme,
                          const Rectangle& constraint) const
    { return distance()(extreme, corner()(constraint, 1)); }
  
    FT compute_y_distance(const Point& extreme,
                          const Rectangle& constraint) const
    { return distance()(extreme, corner()(constraint, 3)); }
  
    Point construct_corner_square(const Rectangle& bbox, FT r) const
    { return pt_a_r()(construct_min_point_2_object()(bbox),
                      construct_min_point_2_object()(bbox), r); }
  
    Point construct_x_square(const Point& p, FT r) const
    { return pt_b_r()(p, p, r); }
  
    Point construct_y_square(const Point& p, FT r) const
    { return pt_a_l()(p, p, r); }
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
  typedef typename R::Construct_min_point_2      Construct_min_point_2;
  typedef typename R::Construct_max_point_2      Construct_max_point_2;
  typedef typename R::Construct_corner_2         Construct_corner_2;
  typedef typename R::Construct_iso_rectangle_2  Construct_iso_rectangle_2;
  typedef typename R::Construct_projection_onto_horizontal_implicit_line_2
    Construct_projection_onto_horizontal_implicit_line_2;
  typedef typename R::Construct_point_2_above_right_implicit_point_2
    Construct_point_2_above_right_implicit_point_2;
  typedef typename R::Construct_point_2_above_left_implicit_point_2
    Construct_point_2_above_left_implicit_point_2;
  typedef typename R::Construct_point_2_below_right_implicit_point_2
    Construct_point_2_below_right_implicit_point_2;
  typedef typename R::Construct_point_2_below_left_implicit_point_2
    Construct_point_2_below_left_implicit_point_2;
  typedef std::binder1st< Infinity_distance_2 >   Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x_2_object() const { return r_.less_x_2_object(); }
  Less_y_2  less_y_2_object() const { return r_.less_y_2_object(); }
  Greater_x_2  greater_x_2_object() const { return r_.greater_x_2_object(); }
  Greater_y_2  greater_y_2_object() const { return r_.greater_y_2_object(); }
  Infinity_distance_2  distance() const
  { return r_.infinity_distance_2_object(); }
  Construct_min_point_2  construct_min_point_2_object() const
  { return r_.construct_min_point_2_object(); }
  Construct_max_point_2  construct_max_point_2_object() const
  { return r_.construct_max_point_2_object(); }
  Construct_corner_2  corner() const
  { return r_.construct_corner_2_object(); }
  Construct_iso_rectangle_2 construct_iso_rectangle_2_object() const
  { return r_.construct_iso_rectangle_2_object(); }
  
  Construct_point_2_below_left_implicit_point_2
  pt_b_l() const
  { return r_.construct_point_2_below_left_implicit_point_2_object(); }
  Construct_point_2_above_left_implicit_point_2
  pt_a_l() const
  { return r_.construct_point_2_above_left_implicit_point_2_object(); }
  Construct_point_2_below_right_implicit_point_2
  pt_b_r() const
  { return r_.construct_point_2_below_right_implicit_point_2_object(); }
  Construct_point_2_above_right_implicit_point_2
  pt_a_r() const
  { return r_.construct_point_2_above_right_implicit_point_2_object(); }
  
  typedef typename R::Min_x_2 Min_x_2;
  Min_x_2 minx() const { return r_.min_x_2_object(); }
  typedef typename R::Min_y_2 Min_y_2;
  Min_y_2 miny() const { return r_.min_y_2_object(); }
  typedef typename R::Max_x_2 Max_x_2;
  Max_x_2 maxx() const { return r_.max_x_2_object(); }
  typedef typename R::Max_y_2 Max_y_2;
  Max_y_2 maxy() const { return r_.max_y_2_object(); }
  
  Construct_projection_onto_horizontal_implicit_line_2
  construct_projection_onto_horizontal_implicit_line_2_object() const
  { return r_.construct_projection_onto_horizontal_implicit_line_2_object(); }
  
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
  : r_(r), delta_(std::bind1st(r.infinity_distance_2_object(), p))
  {}

  X_compare  compare_x() const { return greater_x_2_object(); }
  Y_compare  compare_y() const { return greater_y_2_object(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    maxx()(first_uncovered, construct_max_point_2_object()(constraint));
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      bpt, construct_max_point_2_object()(bbox));
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      constraint_empty ? first_uncovered :
      maxx()(first_uncovered, construct_max_point_2_object()(constraint)),
        construct_max_point_2_object()(bbox));
  #endif
  }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? construct_min_point_2_object()(bbox) :
    construct_max_point_2_object()(constraint);
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      bpt, construct_max_point_2_object()(bbox));
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      constraint_empty ?
      construct_min_point_2_object()(bbox) :
      construct_max_point_2_object()(constraint),
        construct_max_point_2_object()(bbox));
  #endif
  }
  
  Point place_x_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      maxx()(so_far, pt_a_r()(construct_min_point_2_object()(bbox),
                              construct_min_point_2_object()(bbox), radius)),
      so_far);
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    miny()(first_uncovered, construct_min_point_2_object()(constraint));
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_min_point_2_object()(bbox), bpt);
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_min_point_2_object()(bbox),
      constraint_empty ? first_uncovered :
    miny()(first_uncovered, construct_min_point_2_object()(constraint)));
  #endif
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? construct_max_point_2_object()(bbox) :
    construct_min_point_2_object()(constraint);
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_min_point_2_object()(bbox), bpt);
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_min_point_2_object()(bbox),
      constraint_empty ?
      construct_max_point_2_object()(bbox) :
      construct_min_point_2_object()(constraint));
  #endif
  }
  
  Point place_y_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      so_far,
      miny()(so_far, pt_b_l()(construct_max_point_2_object()(bbox),
                              construct_max_point_2_object()(bbox), radius)));
  }
  
  Point update_x_square(const Point& s, const Point& newp) const
  { return construct_projection_onto_horizontal_implicit_line_2_object()(
    maxx()(s, newp), s); }
  
  Point update_y_square(const Point& s, const Point& newp) const
  { return construct_projection_onto_horizontal_implicit_line_2_object()(
    s, miny()(s, newp)); }
  
  FT compute_x_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 0)); }
  
  FT compute_y_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 2)); }
  
  Point construct_corner_square(const Rectangle& bbox, FT r) const
  { return pt_a_l()(construct_max_point_2_object()(bbox),
                    construct_min_point_2_object()(bbox), r); }
  
  Point construct_x_square(const Point& p, FT r) const
  { return pt_b_l()(p, p, r); }
  
  Point construct_y_square(const Point& p, FT r) const
  { return pt_a_r()(p, p, r); }
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
  typedef typename R::Construct_min_point_2      Construct_min_point_2;
  typedef typename R::Construct_max_point_2      Construct_max_point_2;
  typedef typename R::Construct_corner_2         Construct_corner_2;
  typedef typename R::Construct_iso_rectangle_2  Construct_iso_rectangle_2;
  typedef typename R::Construct_projection_onto_horizontal_implicit_line_2
    Construct_projection_onto_horizontal_implicit_line_2;
  typedef typename R::Construct_point_2_above_right_implicit_point_2
    Construct_point_2_above_right_implicit_point_2;
  typedef typename R::Construct_point_2_above_left_implicit_point_2
    Construct_point_2_above_left_implicit_point_2;
  typedef typename R::Construct_point_2_below_right_implicit_point_2
    Construct_point_2_below_right_implicit_point_2;
  typedef typename R::Construct_point_2_below_left_implicit_point_2
    Construct_point_2_below_left_implicit_point_2;
  typedef std::binder1st< Infinity_distance_2 >   Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x_2_object() const { return r_.less_x_2_object(); }
  Less_y_2  less_y_2_object() const { return r_.less_y_2_object(); }
  Greater_x_2  greater_x_2_object() const { return r_.greater_x_2_object(); }
  Greater_y_2  greater_y_2_object() const { return r_.greater_y_2_object(); }
  Infinity_distance_2  distance() const
  { return r_.infinity_distance_2_object(); }
  Construct_min_point_2  construct_min_point_2_object() const
  { return r_.construct_min_point_2_object(); }
  Construct_max_point_2  construct_max_point_2_object() const
  { return r_.construct_max_point_2_object(); }
  Construct_corner_2  corner() const
  { return r_.construct_corner_2_object(); }
  Construct_iso_rectangle_2 construct_iso_rectangle_2_object() const
  { return r_.construct_iso_rectangle_2_object(); }
  
  Construct_point_2_below_left_implicit_point_2
  pt_b_l() const
  { return r_.construct_point_2_below_left_implicit_point_2_object(); }
  Construct_point_2_above_left_implicit_point_2
  pt_a_l() const
  { return r_.construct_point_2_above_left_implicit_point_2_object(); }
  Construct_point_2_below_right_implicit_point_2
  pt_b_r() const
  { return r_.construct_point_2_below_right_implicit_point_2_object(); }
  Construct_point_2_above_right_implicit_point_2
  pt_a_r() const
  { return r_.construct_point_2_above_right_implicit_point_2_object(); }
  
  typedef typename R::Min_x_2 Min_x_2;
  Min_x_2 minx() const { return r_.min_x_2_object(); }
  typedef typename R::Min_y_2 Min_y_2;
  Min_y_2 miny() const { return r_.min_y_2_object(); }
  typedef typename R::Max_x_2 Max_x_2;
  Max_x_2 maxx() const { return r_.max_x_2_object(); }
  typedef typename R::Max_y_2 Max_y_2;
  Max_y_2 maxy() const { return r_.max_y_2_object(); }
  
  Construct_projection_onto_horizontal_implicit_line_2
  construct_projection_onto_horizontal_implicit_line_2_object() const
  { return r_.construct_projection_onto_horizontal_implicit_line_2_object(); }
  
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
  : r_(r), delta_(std::bind1st(r.infinity_distance_2_object(), p))
  {}

  X_compare  compare_x() const { return greater_x_2_object(); }
  Y_compare  compare_y() const { return less_y_2_object(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    maxx()(first_uncovered, construct_max_point_2_object()(constraint));
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      bpt, construct_min_point_2_object()(bbox));
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      constraint_empty ? first_uncovered :
      maxx()(first_uncovered, construct_max_point_2_object()(constraint)),
        construct_min_point_2_object()(bbox));
  #endif
  }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? construct_min_point_2_object()(bbox) :
    construct_max_point_2_object()(constraint);
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      bpt, construct_min_point_2_object()(bbox));
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      constraint_empty ?
      construct_min_point_2_object()(bbox) :
      construct_max_point_2_object()(constraint),
      construct_min_point_2_object()(bbox));
  #endif
  }
  
  Point place_x_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      maxx()(so_far, pt_a_r()(construct_min_point_2_object()(bbox),
                              construct_min_point_2_object()(bbox), radius)),
      so_far);
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    maxy()(first_uncovered, construct_max_point_2_object()(constraint));
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_min_point_2_object()(bbox), bpt);
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_min_point_2_object()(bbox),
      constraint_empty ? first_uncovered :
    maxy()(first_uncovered, construct_max_point_2_object()(constraint)));
  #endif
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? construct_min_point_2_object()(bbox) :
    construct_max_point_2_object()(constraint);
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_min_point_2_object()(bbox), bpt);
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_min_point_2_object()(bbox),
      constraint_empty ?
      construct_min_point_2_object()(bbox) :
      construct_max_point_2_object()(constraint));
  #endif
  }
  
  Point place_y_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      so_far,
      maxy()(pt_a_r()(construct_min_point_2_object()(bbox),
                      construct_min_point_2_object()(bbox), radius), so_far));
  }
  
  Point update_x_square(const Point& s, const Point& newp) const
  { return construct_projection_onto_horizontal_implicit_line_2_object()(
    maxx()(s, newp), s); }
  
  Point update_y_square(const Point& s, const Point& newp) const
  { return construct_projection_onto_horizontal_implicit_line_2_object()(
    s, maxy()(s, newp)); }
  
  FT compute_x_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 3)); }
  
  FT compute_y_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 1)); }
  
  Point construct_corner_square(const Rectangle& bbox, FT r) const
  { return pt_b_l()(construct_max_point_2_object()(bbox),
                    construct_max_point_2_object()(bbox), r); }
  
  Point construct_x_square(const Point& p, FT r) const
  { return pt_a_l()(p, p, r); }
  
  Point construct_y_square(const Point& p, FT r) const
  { return pt_b_r()(p, p, r); }
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
  typedef typename R::Construct_min_point_2      Construct_min_point_2;
  typedef typename R::Construct_max_point_2      Construct_max_point_2;
  typedef typename R::Construct_corner_2         Construct_corner_2;
  typedef typename R::Construct_iso_rectangle_2  Construct_iso_rectangle_2;
  typedef typename R::Construct_projection_onto_horizontal_implicit_line_2
    Construct_projection_onto_horizontal_implicit_line_2;
  typedef typename R::Construct_point_2_above_right_implicit_point_2
    Construct_point_2_above_right_implicit_point_2;
  typedef typename R::Construct_point_2_above_left_implicit_point_2
    Construct_point_2_above_left_implicit_point_2;
  typedef typename R::Construct_point_2_below_right_implicit_point_2
    Construct_point_2_below_right_implicit_point_2;
  typedef typename R::Construct_point_2_below_left_implicit_point_2
    Construct_point_2_below_left_implicit_point_2;
  typedef std::binder1st< Infinity_distance_2 >   Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x_2_object() const { return r_.less_x_2_object(); }
  Less_y_2  less_y_2_object() const { return r_.less_y_2_object(); }
  Greater_x_2  greater_x_2_object() const { return r_.greater_x_2_object(); }
  Greater_y_2  greater_y_2_object() const { return r_.greater_y_2_object(); }
  Infinity_distance_2  distance() const
  { return r_.infinity_distance_2_object(); }
  Construct_min_point_2  construct_min_point_2_object() const
  { return r_.construct_min_point_2_object(); }
  Construct_max_point_2  construct_max_point_2_object() const
  { return r_.construct_max_point_2_object(); }
  Construct_corner_2  corner() const
  { return r_.construct_corner_2_object(); }
  Construct_iso_rectangle_2 construct_iso_rectangle_2_object() const
  { return r_.construct_iso_rectangle_2_object(); }
  
  Construct_point_2_below_left_implicit_point_2
  pt_b_l() const
  { return r_.construct_point_2_below_left_implicit_point_2_object(); }
  Construct_point_2_above_left_implicit_point_2
  pt_a_l() const
  { return r_.construct_point_2_above_left_implicit_point_2_object(); }
  Construct_point_2_below_right_implicit_point_2
  pt_b_r() const
  { return r_.construct_point_2_below_right_implicit_point_2_object(); }
  Construct_point_2_above_right_implicit_point_2
  pt_a_r() const
  { return r_.construct_point_2_above_right_implicit_point_2_object(); }
  
  typedef typename R::Min_x_2 Min_x_2;
  Min_x_2 minx() const { return r_.min_x_2_object(); }
  typedef typename R::Min_y_2 Min_y_2;
  Min_y_2 miny() const { return r_.min_y_2_object(); }
  typedef typename R::Max_x_2 Max_x_2;
  Max_x_2 maxx() const { return r_.max_x_2_object(); }
  typedef typename R::Max_y_2 Max_y_2;
  Max_y_2 maxy() const { return r_.max_y_2_object(); }
  
  Construct_projection_onto_horizontal_implicit_line_2
  construct_projection_onto_horizontal_implicit_line_2_object() const
  { return r_.construct_projection_onto_horizontal_implicit_line_2_object(); }
  
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
  : r_(r), delta_(std::bind1st(r.infinity_distance_2_object(), p))
  {}

  X_compare  compare_x() const { return less_x_2_object(); }
  Y_compare  compare_y() const { return less_y_2_object(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    minx()(first_uncovered, construct_min_point_2_object()(constraint));
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      bpt, construct_min_point_2_object()(bbox));
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      constraint_empty ? first_uncovered :
      minx()(first_uncovered, construct_min_point_2_object()(constraint)),
        construct_min_point_2_object()(bbox));
  #endif
  }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? construct_max_point_2_object()(bbox) :
    construct_min_point_2_object()(constraint);
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      bpt, construct_min_point_2_object()(bbox));
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      constraint_empty ?
      construct_max_point_2_object()(bbox) :
      construct_min_point_2_object()(constraint),
      construct_min_point_2_object()(bbox));
  #endif
  }
  
  Point place_x_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      minx()(pt_b_l()(construct_max_point_2_object()(bbox),
                      construct_max_point_2_object()(bbox), radius), so_far),
      so_far);
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    maxy()(first_uncovered, construct_max_point_2_object()(constraint));
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_max_point_2_object()(bbox), bpt);
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_max_point_2_object()(bbox),
      constraint_empty ? first_uncovered :
    maxy()(first_uncovered, construct_max_point_2_object()(constraint)));
  #endif
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? construct_min_point_2_object()(bbox) :
    construct_max_point_2_object()(constraint);
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_max_point_2_object()(bbox), bpt);
  #else
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      construct_max_point_2_object()(bbox),
      constraint_empty ?
      construct_min_point_2_object()(bbox) :
      construct_max_point_2_object()(constraint));
  #endif
  }
  
  Point place_y_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    return construct_projection_onto_horizontal_implicit_line_2_object()(
      so_far,
      maxy()(pt_a_r()(construct_min_point_2_object()(bbox),
                      construct_min_point_2_object()(bbox), radius), so_far));
  }
  
  Point update_x_square(const Point& s, const Point& newp) const
  { return construct_projection_onto_horizontal_implicit_line_2_object()(
    minx()(s, newp), s); }
  
  Point update_y_square(const Point& s, const Point& newp) const
  { return construct_projection_onto_horizontal_implicit_line_2_object()(
    s, maxy()(s, newp)); }
  
  FT compute_x_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 2)); }
  
  FT compute_y_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, corner()(constraint, 0)); }
  
  Point construct_corner_square(const Rectangle& bbox, FT r) const
  { return pt_b_r()(construct_min_point_2_object()(bbox),
                    construct_max_point_2_object()(bbox), r); }
  
  Point construct_x_square(const Point& p, FT r) const
  { return pt_a_r()(p, p, r); }
  
  Point construct_y_square(const Point& p, FT r) const
  { return pt_b_l()(p, p, r); }
};

#endif // ! (CGAL_RECTANGULAR_3_CENTER_2_MSVC_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

