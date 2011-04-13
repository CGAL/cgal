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
  typedef typename Swap< Less_x_2 >::Type        Greater_x_2;
  typedef typename Swap< Less_y_2 >::Type        Greater_y_2;
  typedef Min< Point_2, Less_x_2 >               Min_x_2;
  typedef Max< Point_2, Less_x_2 >               Max_x_2;
  typedef Min< Point_2, Less_y_2 >               Min_y_2;
  typedef Max< Point_2, Less_y_2 >               Max_y_2;
  typedef typename R::Construct_vertex_2         Construct_vertex_2;
  typedef typename R::Construct_iso_rectangle_2  Construct_iso_rectangle_2;
  typedef typename R::Construct_point_2_above_right_implicit_point_2
    Construct_point_2_above_right_implicit_point_2;
  typedef typename R::Construct_point_2_above_left_implicit_point_2
    Construct_point_2_above_left_implicit_point_2;
  typedef typename R::Construct_point_2_below_right_implicit_point_2
    Construct_point_2_below_right_implicit_point_2;
  typedef typename R::Construct_point_2_below_left_implicit_point_2
    Construct_point_2_below_left_implicit_point_2;
  typedef typename Bind< Infinity_distance_2, Point_2, 1 >::Type  Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x_2_object() const { return r_.less_x_2_object(); }
  Less_y_2  less_y_2_object() const { return r_.less_y_2_object(); }
  Greater_x_2  greater_x_2_object() const
  { return swap_1(less_x_2_object()); }
  Greater_y_2  greater_y_2_object() const
  { return swap_1(less_y_2_object()); }
  Infinity_distance_2  distance() const
  { return r_.infinity_distance_2_object(); }
  Construct_vertex_2  construct_vertex_2_object() const
  { return r_.construct_vertex_2_object(); }
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
  
  Min_x_2 minx() const { return Min_x_2(less_x_2_object()); }
  Min_y_2 miny() const { return Min_y_2(less_y_2_object()); }
  Max_x_2 maxx() const { return Max_x_2(less_x_2_object()); }
  Max_y_2 maxy() const { return Max_y_2(less_y_2_object()); }
  
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
  : r_(r), delta_(bind_1(r.infinity_distance_2_object(), p))
  {}

  X_compare  compare_x() const { return less_x_2_object(); }
  Y_compare  compare_y() const { return greater_y_2_object(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
      Point bpt = constraint_empty ? first_uncovered :
      minx()(first_uncovered, v(constraint, 0));
      return v(rect(bpt, v(bbox, 2)), 3);
  #else
      return v(rect(constraint_empty ? first_uncovered :
                      minx()(first_uncovered, v(constraint, 0)),
                    v(bbox, 2)),
               3);
  #endif
    }
  
    Point place_x_square(bool constraint_empty,
                         const Rectangle& constraint,
                         const Rectangle& bbox) const
    {
      Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
      Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
      Point bpt = constraint_empty ? v(bbox, 2) : v(constraint, 0);
      return v(rect(bpt, v(bbox, 2)), 3);
  #else
      return v(rect(constraint_empty ? v(bbox, 2) : v(constraint, 0),
                    v(bbox, 2)),
               3);
  #endif
    }
  
    Point place_x_square(const Point& so_far,
                         const Rectangle& bbox,
                         FT radius) const
    {
      Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
      Construct_vertex_2        v    = construct_vertex_2_object();
      return v(rect(minx()(pt_b_l()(v(bbox, 2), v(bbox, 2), radius), so_far),
                    so_far),
             3);
    }
  
    Point place_y_square(bool constraint_empty,
                         const Rectangle& constraint,
                         const Point& first_uncovered,
                         const Rectangle& bbox) const
    {
      Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
      Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
      Point bpt = constraint_empty ? first_uncovered :
      miny()(first_uncovered, v(constraint, 0));
      return v(rect(v(bbox, 2), bpt), 1);
  #else
      return v(rect(v(bbox, 2),
                    constraint_empty ? first_uncovered :
                      miny()(first_uncovered, v(constraint, 0))),
               1);
  #endif
    }
  
    Point place_y_square(bool constraint_empty,
                         const Rectangle& constraint,
                         const Rectangle& bbox) const
    {
      Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
      Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
      Point bpt = constraint_empty ? v(bbox, 2) : v(constraint, 0);
      return v(rect(v(bbox, 2), bpt), 1);
  #else
      return v(rect(v(bbox, 2),
                    constraint_empty ? v(bbox, 2) : v(constraint, 0)),
               1);
  #endif
    }
  
    Point place_y_square(const Point& so_far,
                         const Rectangle& bbox,
                         FT radius) const
    {
      Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
      Construct_vertex_2        v    = construct_vertex_2_object();
      return v(rect(so_far,
                    miny()(pt_b_l()(v(bbox, 2), v(bbox, 2), radius),
                           so_far)),
               1);
    }
  
    Point update_x_square(const Point& s, const Point& newp) const
    {
      Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
      Construct_vertex_2        v    = construct_vertex_2_object();
  
      return v(rect(minx()(s, newp), s), 3);
    }
  
    Point update_y_square(const Point& s, const Point& newp) const {
      Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
      Construct_vertex_2        v    = construct_vertex_2_object();
  
      return v(rect(s, miny()(s, newp)), 1);
    }
  
    FT compute_x_distance(const Point& extreme,
                          const Rectangle& constraint) const
    { return distance()(extreme, construct_vertex_2_object()(constraint, 1)); }
  
    FT compute_y_distance(const Point& extreme,
                          const Rectangle& constraint) const
    { return distance()(extreme, construct_vertex_2_object()(constraint, 3)); }
  
    Point construct_corner_square(const Rectangle& bbox, FT r) const
    { return pt_a_r()(construct_vertex_2_object()(bbox, 0),
                      construct_vertex_2_object()(bbox, 0), r); }
  
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
  typedef typename Swap< Less_x_2 >::Type        Greater_x_2;
  typedef typename Swap< Less_y_2 >::Type        Greater_y_2;
  typedef Min< Point_2, Less_x_2 >               Min_x_2;
  typedef Max< Point_2, Less_x_2 >               Max_x_2;
  typedef Min< Point_2, Less_y_2 >               Min_y_2;
  typedef Max< Point_2, Less_y_2 >               Max_y_2;
  typedef typename R::Construct_vertex_2         Construct_vertex_2;
  typedef typename R::Construct_iso_rectangle_2  Construct_iso_rectangle_2;
  typedef typename R::Construct_point_2_above_right_implicit_point_2
    Construct_point_2_above_right_implicit_point_2;
  typedef typename R::Construct_point_2_above_left_implicit_point_2
    Construct_point_2_above_left_implicit_point_2;
  typedef typename R::Construct_point_2_below_right_implicit_point_2
    Construct_point_2_below_right_implicit_point_2;
  typedef typename R::Construct_point_2_below_left_implicit_point_2
    Construct_point_2_below_left_implicit_point_2;
  typedef typename Bind< Infinity_distance_2, Point_2, 1 >::Type  Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x_2_object() const { return r_.less_x_2_object(); }
  Less_y_2  less_y_2_object() const { return r_.less_y_2_object(); }
  Greater_x_2  greater_x_2_object() const
  { return swap_1(less_x_2_object()); }
  Greater_y_2  greater_y_2_object() const
  { return swap_1(less_y_2_object()); }
  Infinity_distance_2  distance() const
  { return r_.infinity_distance_2_object(); }
  Construct_vertex_2  construct_vertex_2_object() const
  { return r_.construct_vertex_2_object(); }
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
  
  Min_x_2 minx() const { return Min_x_2(less_x_2_object()); }
  Min_y_2 miny() const { return Min_y_2(less_y_2_object()); }
  Max_x_2 maxx() const { return Max_x_2(less_x_2_object()); }
  Max_y_2 maxy() const { return Max_y_2(less_y_2_object()); }
  
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
  : r_(r), delta_(bind_1(r.infinity_distance_2_object(), p))
  {}

  X_compare  compare_x() const { return greater_x_2_object(); }
  Y_compare  compare_y() const { return greater_y_2_object(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
      Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
      Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    maxx()(first_uncovered, v(constraint, 2));
    return v(rect(bpt, v(bbox, 3)), 2);
  #else
    return v(rect(constraint_empty ? first_uncovered :
                    maxx()(first_uncovered, v(constraint, 2)),
                  v(bbox, 3)),
             2);
  #endif
  }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? v(bbox, 0) : v(constraint, 2);
    return v(rect(bpt, v(bbox, 3)), 2);
  #else
    return v(rect(constraint_empty ? v(bbox, 0) : v(constraint, 2),
                  v(bbox, 3)),
             2);
  #endif
  }
  
  Point place_x_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(maxx()(so_far, pt_a_r()(v(bbox, 0),
                                          v(bbox, 0), radius)),
                  so_far),
             2);
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    miny()(first_uncovered, v(constraint, 0));
    return v(rect(v(bbox, 3), bpt), 0);
  #else
    return v(rect(v(bbox, 3),
                  constraint_empty ? first_uncovered :
                    miny()(first_uncovered, v(constraint, 0))),
             0);
  #endif
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? v(bbox, 2) : v(constraint, 0);
    return v(rect(v(bbox, 3), bpt), 0);
  #else
    return v(rect(v(bbox, 3),
                  constraint_empty ? v(bbox, 2) : v(constraint, 0)),
             0);
  #endif
  }
  
  Point place_y_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(so_far,
                  miny()(so_far, pt_b_l()(v(bbox, 2),
                                          v(bbox, 2), radius))),
             0);
  }
  
  Point update_x_square(const Point& s, const Point& newp) const {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(maxx()(s, newp), s), 2);
  }
  
  Point update_y_square(const Point& s, const Point& newp) const {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(s, miny()(s, newp)), 0);
  }
  
  FT compute_x_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, construct_vertex_2_object()(constraint, 0)); }
  
  FT compute_y_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, construct_vertex_2_object()(constraint, 2)); }
  
  Point construct_corner_square(const Rectangle& bbox, FT r) const
  { return pt_a_l()(construct_vertex_2_object()(bbox, 2),
                    construct_vertex_2_object()(bbox, 0), r); }
  
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
  typedef typename Swap< Less_x_2 >::Type        Greater_x_2;
  typedef typename Swap< Less_y_2 >::Type        Greater_y_2;
  typedef Min< Point_2, Less_x_2 >               Min_x_2;
  typedef Max< Point_2, Less_x_2 >               Max_x_2;
  typedef Min< Point_2, Less_y_2 >               Min_y_2;
  typedef Max< Point_2, Less_y_2 >               Max_y_2;
  typedef typename R::Construct_vertex_2         Construct_vertex_2;
  typedef typename R::Construct_iso_rectangle_2  Construct_iso_rectangle_2;
  typedef typename R::Construct_point_2_above_right_implicit_point_2
    Construct_point_2_above_right_implicit_point_2;
  typedef typename R::Construct_point_2_above_left_implicit_point_2
    Construct_point_2_above_left_implicit_point_2;
  typedef typename R::Construct_point_2_below_right_implicit_point_2
    Construct_point_2_below_right_implicit_point_2;
  typedef typename R::Construct_point_2_below_left_implicit_point_2
    Construct_point_2_below_left_implicit_point_2;
  typedef typename Bind< Infinity_distance_2, Point_2, 1 >::Type  Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x_2_object() const { return r_.less_x_2_object(); }
  Less_y_2  less_y_2_object() const { return r_.less_y_2_object(); }
  Greater_x_2  greater_x_2_object() const
  { return swap_1(less_x_2_object()); }
  Greater_y_2  greater_y_2_object() const
  { return swap_1(less_y_2_object()); }
  Infinity_distance_2  distance() const
  { return r_.infinity_distance_2_object(); }
  Construct_vertex_2  construct_vertex_2_object() const
  { return r_.construct_vertex_2_object(); }
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
  
  Min_x_2 minx() const { return Min_x_2(less_x_2_object()); }
  Min_y_2 miny() const { return Min_y_2(less_y_2_object()); }
  Max_x_2 maxx() const { return Max_x_2(less_x_2_object()); }
  Max_y_2 maxy() const { return Max_y_2(less_y_2_object()); }
  
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
  : r_(r), delta_(bind_1(r.infinity_distance_2_object(), p))
  {}

  X_compare  compare_x() const { return greater_x_2_object(); }
  Y_compare  compare_y() const { return less_y_2_object(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    maxx()(first_uncovered, v(constraint, 2));
    return v(rect(bpt, v(bbox, 0)), 1);
  #else
    return v(rect(constraint_empty ? first_uncovered :
                    maxx()(first_uncovered, v(constraint, 2)),
                  v(bbox, 0)),
             1);
  #endif
  }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? v(bbox, 0) : v(constraint, 2);
    return v(rect(bpt, v(bbox, 0)), 1);
  #else
    return v(rect(constraint_empty ? v(bbox, 0) : v(constraint, 2),
                  v(bbox, 0)),
             1);
  #endif
  }
  
  Point place_x_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(maxx()(so_far, pt_a_r()(v(bbox, 0),
                                          v(bbox, 0), radius)),
                  so_far),
             1);
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    maxy()(first_uncovered, v(constraint, 2));
    return v(rect(v(bbox, 0), bpt), 3);
  #else
    return v(rect(v(bbox, 0),
                  constraint_empty ? first_uncovered :
                    maxy()(first_uncovered, v(constraint, 2))),
             3);
  #endif
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? v(bbox, 0) : v(constraint, 2);
    return v(rect(v(bbox, 0), bpt), 3);
  #else
    return v(rect(v(bbox, 0),
                  constraint_empty ? v(bbox, 0) : v(constraint, 2)),
             3);
  #endif
  }
  
  Point place_y_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(so_far,
                  maxy()(pt_a_r()(v(bbox, 0),
                                  v(bbox, 0), radius), so_far)),
             3);
  }
  
  Point update_x_square(const Point& s, const Point& newp) const {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(maxx()(s, newp), s), 1);
  }
  
  Point update_y_square(const Point& s, const Point& newp) const {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(s, maxy()(s, newp)), 3);
  }
  
  FT compute_x_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, construct_vertex_2_object()(constraint, 3)); }
  
  FT compute_y_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, construct_vertex_2_object()(constraint, 1)); }
  
  Point construct_corner_square(const Rectangle& bbox, FT r) const
  { return pt_b_l()(construct_vertex_2_object()(bbox, 2),
                    construct_vertex_2_object()(bbox, 2), r); }
  
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
  typedef typename Swap< Less_x_2 >::Type        Greater_x_2;
  typedef typename Swap< Less_y_2 >::Type        Greater_y_2;
  typedef Min< Point_2, Less_x_2 >               Min_x_2;
  typedef Max< Point_2, Less_x_2 >               Max_x_2;
  typedef Min< Point_2, Less_y_2 >               Min_y_2;
  typedef Max< Point_2, Less_y_2 >               Max_y_2;
  typedef typename R::Construct_vertex_2         Construct_vertex_2;
  typedef typename R::Construct_iso_rectangle_2  Construct_iso_rectangle_2;
  typedef typename R::Construct_point_2_above_right_implicit_point_2
    Construct_point_2_above_right_implicit_point_2;
  typedef typename R::Construct_point_2_above_left_implicit_point_2
    Construct_point_2_above_left_implicit_point_2;
  typedef typename R::Construct_point_2_below_right_implicit_point_2
    Construct_point_2_below_right_implicit_point_2;
  typedef typename R::Construct_point_2_below_left_implicit_point_2
    Construct_point_2_below_left_implicit_point_2;
  typedef typename Bind< Infinity_distance_2, Point_2, 1 >::Type  Delta;
  
  Delta  delta() const { return delta_; }
  Less_x_2  less_x_2_object() const { return r_.less_x_2_object(); }
  Less_y_2  less_y_2_object() const { return r_.less_y_2_object(); }
  Greater_x_2  greater_x_2_object() const
  { return swap_1(less_x_2_object()); }
  Greater_y_2  greater_y_2_object() const
  { return swap_1(less_y_2_object()); }
  Infinity_distance_2  distance() const
  { return r_.infinity_distance_2_object(); }
  Construct_vertex_2  construct_vertex_2_object() const
  { return r_.construct_vertex_2_object(); }
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
  
  Min_x_2 minx() const { return Min_x_2(less_x_2_object()); }
  Min_y_2 miny() const { return Min_y_2(less_y_2_object()); }
  Max_x_2 maxx() const { return Max_x_2(less_x_2_object()); }
  Max_y_2 maxy() const { return Max_y_2(less_y_2_object()); }
  
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
  : r_(r), delta_(bind_1(r.infinity_distance_2_object(), p))
  {}

  X_compare  compare_x() const { return less_x_2_object(); }
  Y_compare  compare_y() const { return less_y_2_object(); }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    minx()(first_uncovered, v(constraint, 0));
    return v(rect(bpt, v(bbox, 1)), 0);
  #else
    return v(rect(constraint_empty ? first_uncovered :
                    minx()(first_uncovered, v(constraint, 0)),
                  v(bbox, 1)),
             0);
  #endif
  }
  
  Point place_x_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? v(bbox, 2) : v(constraint, 0);
    return v(rect(bpt, v(bbox, 1)), 0);
  #else
    return v(rect(constraint_empty ? v(bbox, 2) : v(constraint, 0),
                  v(bbox, 1)),
             0);
  #endif
  }
  
  Point place_x_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(minx()(pt_b_l()(v(bbox, 2),
                                  v(bbox, 2), radius), so_far),
                  so_far),
             0);
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Point& first_uncovered,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? first_uncovered :
    maxy()(first_uncovered, v(constraint, 2));
    return v(rect(v(bbox, 1), bpt), 2);
  #else
    return v(rect(v(bbox, 1),
                  constraint_empty ? first_uncovered :
                    maxy()(first_uncovered, v(constraint, 2))),
             2);
  #endif
  }
  
  Point place_y_square(bool constraint_empty,
                       const Rectangle& constraint,
                       const Rectangle& bbox) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
  #ifdef __BORLANDC__
    Point bpt = constraint_empty ? v(bbox, 0) : v(constraint, 2);
    return v(rect(v(bbox, 1), bpt), 2);
  #else
    return v(rect(v(bbox, 1),
                  constraint_empty ? v(bbox, 0) : v(constraint, 2)),
             2);
  #endif
  }
  
  Point place_y_square(const Point& so_far,
                       const Rectangle& bbox,
                       FT radius) const
  {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(so_far,
                  maxy()(pt_a_r()(v(bbox, 0),
                                  v(bbox, 0), radius), so_far)),
             2);
  }
  
  Point update_x_square(const Point& s, const Point& newp) const {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(minx()(s, newp), s), 0);
  }
  
  Point update_y_square(const Point& s, const Point& newp) const {
    Construct_iso_rectangle_2 rect = construct_iso_rectangle_2_object();
    Construct_vertex_2        v    = construct_vertex_2_object();
    return v(rect(s, maxy()(s, newp)), 2);
  }
  
  FT compute_x_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, construct_vertex_2_object()(constraint, 2)); }
  
  FT compute_y_distance(const Point& extreme,
                        const Rectangle& constraint) const
  { return distance()(extreme, construct_vertex_2_object()(constraint, 0)); }
  
  Point construct_corner_square(const Rectangle& bbox, FT r) const
  { return pt_b_r()(construct_vertex_2_object()(bbox, 0),
                    construct_vertex_2_object()(bbox, 2), r); }
  
  Point construct_x_square(const Point& p, FT r) const
  { return pt_a_r()(p, p, r); }
  
  Point construct_y_square(const Point& p, FT r) const
  { return pt_b_l()(p, p, r); }
};

#endif // ! (CGAL_RECTANGULAR_3_CENTER_2_MSVC_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

