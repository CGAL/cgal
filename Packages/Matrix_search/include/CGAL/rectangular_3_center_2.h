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
// file          : rectangular_3_center_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : 3cover.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2,3-Center Computation for Axis-Parallel 2D-Rectangles
// ============================================================================

#if ! (CGAL_RECTANGULAR_3_CENTER_2_H)
#define CGAL_RECTANGULAR_3_CENTER_2_H 1

#include <CGAL/Cartesian.h>
#include <CGAL/function_objects.h>
#include <CGAL/algorithm.h>
#ifdef CGAL_REP_CLASS_DEFINED
#include <CGAL/Rectangular_p_center_traits_2.h>
#endif // CGAL_REP_CLASS_DEFINED
#include <algorithm>
#include <vector>
#include <functional>

#ifdef _MSC_VER
// that compiler cannot even distinguish between global
// and class scope, so ...
#define Base B_B_Base
#endif // _MSC_VER

CGAL_BEGIN_NAMESPACE
/*
struct Wastebasket
: public CGAL_STD::iterator< std::output_iterator_tag, void >
{
  typedef Wastebasket< T > iterator;

  template < class T >
  iterator operator=( const T&) { return *this; }

  iterator operator*()     { return *this; }
  iterator operator++()    { return *this; }
  iterator operator++(int) { return *this; }
};
*/

template < class ForwardIterator, class OutputIterator,
           class FT, class Traits >
OutputIterator
rectangular_2_center_2(
  ForwardIterator f,
  ForwardIterator l,
  OutputIterator o,
  FT& r,
  Traits& t)
{
#ifndef CGAL_CFG_NO_NAMESPACE
  using std::pair;
  using std::greater;
  using std::less;
  using std::bind1st;
  using std::binder1st;
  using std::binary_compose;
#endif
  typedef typename Traits::Iso_rectangle_2     Rectangle;
  typedef typename Traits::Point_2             Point;
  typedef typename Traits::Infinity_distance_2 Dist;
  typedef typename Traits::Construct_min_2     Min_pt;
  typedef typename Traits::Construct_max_2     Max_pt;
  typedef typename Traits::Construct_corner_2  Corner;
  typedef typename Traits::Construct_iso_rectangle_2_above_right_point_2
    C_rect_above_right;
  typedef typename Traits::Construct_iso_rectangle_2_above_left_point_2
    C_rect_above_left;
  typedef typename Traits::Construct_iso_rectangle_2_below_right_point_2
    C_rect_below_right;
  typedef typename Traits::Construct_iso_rectangle_2_below_left_point_2
    C_rect_below_left;
  typedef binary_compose< Min< FT >, binder1st< Dist >, binder1st< Dist > >
    Gamma;

  // fetch function objects from traits class
  Min_pt              minpt    = t.get_construct_min_2();
  Max_pt              maxpt    = t.get_construct_max_2();
  Corner              corner   = t.get_construct_corner_2();
  Dist                dist     = t.get_infinity_distance_2();
  C_rect_above_left   rect_a_l =
    t.get_construct_iso_rectangle_2_above_left_point_2();
  C_rect_above_right  rect_a_r =
    t.get_construct_iso_rectangle_2_above_right_point_2();
  C_rect_below_left   rect_b_l =
    t.get_construct_iso_rectangle_2_below_left_point_2();
  C_rect_below_right  rect_b_r =
    t.get_construct_iso_rectangle_2_below_right_point_2();

  // compute bounding box
  Rectangle bb = bounding_box_2(f, l);

  // two cases: top-left & bottom-right or top-right & bottom-left
  Min< FT > minft;
  Gamma gamma1(minft, bind1st(dist, corner(bb, 0)),
               bind1st(dist, corner(bb, 2)));
  Gamma gamma2(minft, bind1st(dist, corner(bb, 1)),
               bind1st(dist, corner(bb, 3)));

  pair< ForwardIterator, ForwardIterator > cand =
    min_max_element(f, l,
                    compose2_2(greater< FT >(), gamma1, gamma1),
                    compose2_2(less< FT >(), gamma2, gamma2));

  // return the result
  r = gamma1(*cand.first);
  Rectangle r1 = rect_a_r(corner(bb, 0), r);
  Rectangle r2 = rect_b_l(corner(bb, 2), r);
  if (gamma1(*cand.first) >= gamma2(*cand.second)) {
    r = gamma2(*cand.second);
    r1 = rect_a_l(corner(bb, 1), r);
    r2 = rect_b_r(corner(bb, 3), r);
  }
  *o++ = midpoint(minpt(r1), maxpt(r1));
  *o++ = midpoint(minpt(r2), maxpt(r2));
  return o;
}
template < class RandomAccessIterator,
           class OutputIterator,
           class Traits >
OutputIterator
rectangular_3_center_2_type1(
  RandomAccessIterator f,
  RandomAccessIterator l,
  const typename Traits::Iso_rectangle_2& r,
  OutputIterator o,
  typename Traits::FT& rad,
  Traits& t)
{
#ifndef CGAL_CFG_NO_NAMESPACE
  using std::max;
  using std::less;
  using std::nth_element;
  using std::bind1st;
  using std::binder1st;
  using std::binary_compose;
#endif
  typedef typename Traits::FT                         FT;
  typedef typename Traits::Iso_rectangle_2            Rectangle;
  typedef typename Traits::Point_2                    Point;
  typedef typename Traits::Infinity_distance_2        Dist;
  typedef typename Traits::Signed_infinity_distance_2 Sdist;
  typedef typename Traits::Construct_min_2            Min_pt;
  typedef typename Traits::Construct_max_2            Max_pt;
  typedef typename Traits::Construct_corner_2         Corner;
  typedef typename Traits::Construct_iso_rectangle_2_above_right_point_2
    C_rect_above_right;
  typedef typename Traits::Construct_iso_rectangle_2_above_left_point_2
    C_rect_above_left;
  typedef typename Traits::Construct_iso_rectangle_2_below_right_point_2
    C_rect_below_right;
  typedef typename Traits::Construct_iso_rectangle_2_below_left_point_2
    C_rect_below_left;
  typedef binary_compose< Min< FT >, binder1st< Dist >, binder1st< Dist > >
    Gamma;

  // fetch function objects from traits class
  Min_pt              minpt    = t.get_construct_min_2();
  Max_pt              maxpt    = t.get_construct_max_2();
  Corner              corner   = t.get_construct_corner_2();
  Dist                dist     = t.get_infinity_distance_2();
  Sdist               sdist    = t.get_signed_infinity_distance_2();
  C_rect_above_left   rect_a_l =
    t.get_construct_iso_rectangle_2_above_left_point_2();
  C_rect_above_right  rect_a_r =
    t.get_construct_iso_rectangle_2_above_right_point_2();
  C_rect_below_left   rect_b_l =
    t.get_construct_iso_rectangle_2_below_left_point_2();
  C_rect_below_right  rect_b_r =
    t.get_construct_iso_rectangle_2_below_right_point_2();

  // initialize best radius so far
  rad = sdist(maxpt(r), minpt(r));
  // init to prevent default constructor requirement
  Point bestpoint = *f;
  unsigned int bestrun;

  // two cases: top-left & bottom-right or top-right & bottom-left

  // init to prevent default constructor requirement
  Rectangle b = bounding_box_2(f, f + 1);
  for (unsigned int i = 0; i < 2; ++i) {

    // the range [s, e) defines the point set Pt
    RandomAccessIterator s = f;
    RandomAccessIterator e = l;
    bool b_empty = true;
    Min< FT > minft;
    Gamma gamma(minft,
                bind1st(dist, corner(r, i)),
                bind1st(dist, corner(r, 2 + i)));

    while (e - s > 1) {
      // step (a)
      RandomAccessIterator m = s + (e - s - 1) / 2;
      nth_element(s, m, e, compose2_2(less< FT >(), gamma, gamma));

      // step (b)
      Rectangle b_prime = bounding_box_2(m + 1, e);
      if (!b_empty)
        b_prime = construct_bounding_box_union(b, b_prime);

      // step (c) / (d)
      if (sdist(maxpt(b_prime), minpt(b_prime)) > gamma(*m))
        s = m + 1;
      else {
        e = m + 1;
        b_empty = false;
        b = b_prime;
      }
    }

    // check whether s or (s-1) is the solution
    Rectangle b_prime =
      construct_bounding_box_union(b, Rectangle(*s, *s));
    FT b_prime_size = sdist(maxpt(b_prime), minpt(b_prime));
    if (b_prime_size < gamma(*s)) {
      if (b_prime_size < rad) {
        rad = b_prime_size;
        bestpoint = midpoint(minpt(b_prime), maxpt(b_prime));
        bestrun = i;
      }
    } else if (gamma(*s) < rad) {
      rad = gamma(*s);
      bestpoint = midpoint(minpt(b), maxpt(b));
      bestrun = i;
    }
  }

  // return the result
  *o++ = bestpoint;
  Rectangle r1 = rect_a_r(corner(r, 0), rad);
  Rectangle r2 = rect_b_l(corner(r, 2), rad);
  if (bestrun == 1) {
    r1 = rect_a_l(corner(r, 1), rad);
    r2 = rect_b_r(corner(r, 3), rad);
  }
  *o++ = midpoint(minpt(r1), maxpt(r1));
  *o++ = midpoint(minpt(r2), maxpt(r2));
  return o;
}
#ifndef _MSC_VER

template < class R >
struct Rectangular_3_center_2_type2_operations_base {

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

  Rectangular_3_center_2_type2_operations_base(R& r, const Point_2& p)
  : r_(r), delta_(std::bind1st(r.get_infinity_distance_2(), p))
  {}

};
template < class R >
struct Rectangular_3_center_2_type2_operations0
: public Rectangular_3_center_2_type2_operations_base< R >
{
  typedef Rectangular_3_center_2_type2_operations0< R >      This;
  typedef Rectangular_3_center_2_type2_operations_base< R >  Base;
  typedef typename Base::FT                   FT;
  typedef typename Base::Point_2              Point;
  typedef typename Base::Iso_rectangle_2      Rectangle;
  typedef typename Base::Less_x_2             X_compare;
  typedef typename Base::Greater_y_2          Y_compare;
  typedef typename Base::Infinity_distance_2  Distance;

  Rectangular_3_center_2_type2_operations0(R& r, const Point& p)
  : Rectangular_3_center_2_type2_operations_base< R >(r, p)
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
struct Rectangular_3_center_2_type2_operations1
: public Rectangular_3_center_2_type2_operations_base< R >
{
  typedef Rectangular_3_center_2_type2_operations_base< R >  Base;
  typedef typename Base::FT                   FT;
  typedef typename Base::Point_2              Point;
  typedef typename Base::Iso_rectangle_2      Rectangle;
  typedef typename Base::Infinity_distance_2  Distance;
  typedef typename Base::Greater_x_2          X_compare;
  typedef typename Base::Greater_y_2          Y_compare;

  Rectangular_3_center_2_type2_operations1(R& r, const Point& p)
  : Rectangular_3_center_2_type2_operations_base< R >(r, p)
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
struct Rectangular_3_center_2_type2_operations2
: public Rectangular_3_center_2_type2_operations_base< R >
{
  typedef Rectangular_3_center_2_type2_operations_base< R >  Base;
  typedef typename Base::FT                   FT;
  typedef typename Base::Point_2              Point;
  typedef typename Base::Iso_rectangle_2      Rectangle;
  typedef typename Base::Infinity_distance_2  Distance;
  typedef typename Base::Greater_x_2          X_compare;
  typedef typename Base::Less_y_2             Y_compare;

  Rectangular_3_center_2_type2_operations2(R& r, const Point& p)
  : Rectangular_3_center_2_type2_operations_base< R >(r, p)
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
struct Rectangular_3_center_2_type2_operations3
: public Rectangular_3_center_2_type2_operations_base< R >
{
  typedef Rectangular_3_center_2_type2_operations_base< R >  Base;
  typedef typename Base::FT                   FT;
  typedef typename Base::Point_2              Point;
  typedef typename Base::Iso_rectangle_2      Rectangle;
  typedef typename Base::Infinity_distance_2  Distance;
  typedef typename Base::Less_x_2             X_compare;
  typedef typename Base::Less_y_2             Y_compare;

  Rectangular_3_center_2_type2_operations3(R& r, const Point& p)
  : Rectangular_3_center_2_type2_operations_base< R >(r, p)
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

#else

#include <CGAL/rectangular_3_center_2_msvc.h>

#endif // _MSC_VER
template < class RandomAccessIterator,
           class Rectangle,
           class OutputIterator,
           class FT,
           class Operations >
OutputIterator
rectangular_3_center_2_type2(
  RandomAccessIterator f,
  RandomAccessIterator l,
  const Rectangle& r,
  OutputIterator o,
  FT& rad,
  Operations op)
{
#ifndef CGAL_CFG_NO_NAMESPACE
  using std::max;
  using std::less;
  using std::greater_equal;
  using std::not_equal_to;
  using std::logical_and;
  using std::max_element;
  using std::nth_element;
  using std::find_if;
  using std::sort;
  using std::partition;
  using std::pair;
  using std::binder1st;
  using std::bind1st;
  using std::compose1;
  using std::compose2;
  using std::binary_compose;
#endif

  typedef typename Operations::Point                       Point;
  typedef typename Operations::Distance                 Distance;
  typedef pair< RandomAccessIterator, RandomAccessIterator >  IP;

  // constant fraction to be excluded on every iteration (1/.)
  const unsigned int fraction = 7;

  // the range [s, e) defines the point set P
  RandomAccessIterator s = f;
  RandomAccessIterator e = l;
  // a flag to indicate whether we need to search the radius
  // in the final brute-force step or not
  bool radius_is_known = false;

  // the bounding boxes of assigned points
  Rectangle Q_t, Q_r;
  bool Q_t_empty = true, Q_r_empty = true;

  // lower bound for the diameter (2 * radius)
  // also store the corresponding positions of q_t and q_r
  FT rho_max = 0, rho_min = 0, q_t_q_r_cover_at_rho_min = 0;
  Point q_t_at_rho_max, q_r_at_rho_max, q_t_at_rho_min, q_r_at_rho_min;
  RandomAccessIterator s_at_rho_min = s, e_at_rho_min = s;

#ifndef CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
  {
    // First try whether the best radius so far can be reached at all
    RandomAccessIterator m =
      partition(f, l, compose1(bind1st(greater_equal< FT >(), rad),
                               op.delta()));
    IP pos = min_max_element(m, l, op.compare_x(), op.compare_y());
    // extreme points of the two other squares
    Point q_t =
      op.place_x_square(op.place_x_square(Q_t_empty, Q_t, *pos.first, r),
                        r,
                        rad);
    Point q_r =
      op.place_y_square(op.place_y_square(Q_r_empty, Q_r, *pos.second, r),
                        r,
                        rad);
    binder1st< greater_equal< FT > >
      le_rad = bind1st(greater_equal< FT >(), rad);
    RandomAccessIterator b1 =
      partition(m, l, compose1(le_rad, bind1st(op.distance(), q_t)));
    RandomAccessIterator b2 =
      partition(b1, l, compose1(le_rad, bind1st(op.distance(), q_r)));

    if (b2 != l)
      return o;
  }
#endif // ! CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST

#ifndef CGAL_3COVER_NO_PREFILTER
  // Prefiltering heuristic
  while (e - s > 6) {
    int cutoff = (e - s) / 2;
    RandomAccessIterator m = s + cutoff - 1;
    nth_element(s, m, e, compose2_2(less< FT >(), op.delta(), op.delta()));

    // step (b)
    IP pos = min_max_element(m + 1, e, op.compare_x(), op.compare_y());
    // extreme points of the two other squares
    // (consider Q_i and pos [move as far as possible])
    Point q_t_afap = op.place_x_square(Q_t_empty, Q_t, *pos.first, r);
    Point q_r_afap = op.place_y_square(Q_r_empty, Q_r, *pos.second, r);
    // now consider also that we must not leave the bbox r
    Point q_t = op.place_x_square(q_t_afap, r, op.delta()(*m));
    Point q_r = op.place_y_square(q_r_afap, r, op.delta()(*m));

    // check for covering
    binder1st< greater_equal< FT > >
      le_delta_m(bind1st(greater_equal< FT >(), op.delta()(*m)));
    RandomAccessIterator b1 =
      partition(m + 1, e, compose1(le_delta_m, bind1st(op.distance(), q_t)));
    RandomAccessIterator b2 =
      partition(b1, e, compose1(le_delta_m, bind1st(op.distance(), q_r)));

    if (b2 != e)
      s = m;
    else
      break;
  }
#endif // ! CGAL_3COVER_NO_PREFILTER

#ifdef CGAL_PCENTER_DEBUG
  leda_window W(730, 690);
  cgalize(W);
  W.init(-1.5, 1.5, -1.2);
  W.set_node_width(3);
  W.display();
  typedef Ostream_iterator< leda_window > Window_stream_iterator;
  Window_stream_iterator wout(W);
#endif // CGAL_PCENTER_DEBUG

  while (e - s > 6) {
    /*
    cerr << e - s << " points (" << e - s - (e - s) / fraction
         << ")" << endl;
    */
    // step (a)
    int cutoff = (e - s) / fraction;
    RandomAccessIterator m = s + cutoff - 1;
    nth_element(s, m, e, compose2_2(less< FT >(), op.delta(), op.delta()));

    // step (b)
    IP pos = min_max_element(m + 1, e, op.compare_x(), op.compare_y());
    // extreme points of the two other squares
    // (consider Q_i and pos [move as far as possible])
    Point q_t_afap = op.place_x_square(Q_t_empty, Q_t, *pos.first, r);
    Point q_r_afap = op.place_y_square(Q_r_empty, Q_r, *pos.second, r);
    // now consider also that we must not leave the bbox r
    Point q_t = op.place_x_square(q_t_afap, r, op.delta()(*m));
    Point q_r = op.place_y_square(q_r_afap, r, op.delta()(*m));

#ifdef CGAL_PCENTER_DEBUG
    W << GREEN << r;
    if (!Q_t_empty)
      W << ORANGE << Q_t;
    if (!Q_r_empty)
      W << Q_r;
    W << BLACK << Rectangle(r[0], Point(r[0].x() + op.delta()(*m),
                                        r[0].y() + op.delta()(*m)));
    W.set_node_width(5);
    W << BLUE;
    std::copy(s, e, wout);
    W.set_node_width(2);
    W << WHITE << *pos.first << *pos.second;
    { Point p; W >> p; }
#endif // CGAL_PCENTER_DEBUG

#ifdef CGAL_3COVER_CHECK
    // check whether the points in [e,l) which have been assigned
    // to Q_t and Q_r are covered by q_t and q_r
    if ((Q_t_empty || op.compute_x_distance(q_t, Q_t) <= op.delta()(*m)) &&
        (Q_r_empty || op.compute_y_distance(q_r, Q_r) <= op.delta()(*m))) {
      binder1st< less< FT > >
        greater_delta_m(bind1st(less< FT >(), op.delta()(*m)));
      RandomAccessIterator iii =
        find_if(e,
                l,
                compose2(logical_and< bool >(),
                         compose1(greater_delta_m,
                                  bind1st(op.distance(), q_t)),
                         compose1(greater_delta_m,
                                  bind1st(op.distance(), q_r))));
      CGAL_assertion(iii == l);
    }
    // check whether the points in [f,s) are covered
    {
      binder1st< greater_equal< FT > >
      le_delta_m(bind1st(greater_equal< FT >(), op.delta()(*m)));
      RandomAccessIterator iii =
        partition(f, s, compose1(le_delta_m, op.delta()));
      iii = partition(iii, s,
                      compose1(le_delta_m, bind1st(op.distance(), q_t)));
      iii = partition(iii, s,
                      compose1(le_delta_m, bind1st(op.distance(), q_r)));
      CGAL_assertion(iii == s);
    }
#endif // CGAL_3COVER_CHECK

    // partition the range [m+1, e) into ranges
    // [m+1, b1), [b1, b2),   [b2, b3) and [b3, e)
    //     R      G cap q_t  G cap q_r      none
    binder1st< greater_equal< FT > >
    le_delta_m(bind1st(greater_equal< FT >(), op.delta()(*m)));
    RandomAccessIterator b2 =
      partition(m + 1, e, compose1(le_delta_m, bind1st(op.distance(), q_t)));
    RandomAccessIterator b1 =
      partition(m + 1, b2, compose1(le_delta_m, bind1st(op.distance(), q_r)));
    RandomAccessIterator b3 =
      partition(b2, e, compose1(le_delta_m, bind1st(op.distance(), q_r)));

#ifdef CGAL_PCENTER_DEBUG
    W << RED; std::copy(m+1, b1, wout);
    W << GREEN; std::copy(b1, b2, wout);
    W << ORANGE; std::copy(b2, b3, wout);
    W << VIOLET; std::copy(b3, e, wout);
    W << BLACK << q_r;
    { Point p; W << WHITE >> p; }
#endif // CGAL_PCENTER_DEBUG

    // step (c)
    if (b3 != e ||
        !Q_t_empty && op.compute_x_distance(q_t, Q_t) > op.delta()(*m) ||
        !Q_r_empty && op.compute_y_distance(q_r, Q_r) > op.delta()(*m))
      {
        // not covered
        s = b1;
        rho_min = op.delta()(*m);
        q_t_q_r_cover_at_rho_min = 0;
        if (!Q_t_empty)
          q_t_q_r_cover_at_rho_min =
            max(q_t_q_r_cover_at_rho_min, op.compute_x_distance(q_t, Q_t));
        if (!Q_r_empty)
          q_t_q_r_cover_at_rho_min =
            max(q_t_q_r_cover_at_rho_min, op.compute_y_distance(q_r, Q_r));
        q_t_at_rho_min = q_t, q_r_at_rho_min = q_r;
        s_at_rho_min = s, e_at_rho_min = e;
        continue;
      }

    // step (d) [covered]
    if (b3 - b1 >= cutoff) { // enough points in G
      e = b1;
      // adjust Q_t
      if (b1 != b2)
        if (Q_t_empty) {
          Q_t = bounding_box_2(b1, b2);
          Q_t_empty = false;
        } else
          Q_t = construct_bounding_box_union(Q_t, bounding_box_2(b1, b2));
      // adjust Q_r
      if (b2 != b3)
        if (Q_r_empty) {
          Q_r = bounding_box_2(b2, b3);
          Q_r_empty = false;
        } else
          Q_r = construct_bounding_box_union(Q_r, bounding_box_2(b2, b3));
      continue;
    }

    // step (e) [not enough points in G]
    CGAL_assertion(b1 - (m + 1) >= 5 * cutoff);

    // compute the four cutting lines for R
    nth_element(m + 1, m + 1 + cutoff, b1, op.less_x());
    Point x_min_cutoff = *(m + 1 + cutoff);
    nth_element(m + 1, m + 1 + cutoff, b1, op.greater_x());
    Point x_max_cutoff = *(m + 1 + cutoff);
    nth_element(m + 1, m + 1 + cutoff, b1, op.less_y());
    Point y_min_cutoff = *(m + 1 + cutoff);
    nth_element(m + 1, m + 1 + cutoff, b1, op.greater_y());
    Point y_max_cutoff = *(m + 1 + cutoff);
    Rectangle B = op.rectangle()(x_min_cutoff, y_min_cutoff,
                                 x_max_cutoff, y_max_cutoff);

    // Algorithm search_E:

    // the range [s_b, s_e) defines the point set S
    RandomAccessIterator s_b = s;
    RandomAccessIterator s_e = m + 1;

    while (s_e - s_b > 1) {
      // step 1
      RandomAccessIterator s_m = s_b + (s_e - s_b - 1) / 2;
      nth_element(s_b, s_m, s_e,
                  compose2_2(less< FT >(), op.delta(), op.delta()));

      // step 2 (as above)
      Point q_t_m = q_t_afap;
      Point q_r_m = q_r_afap;
      if (s_m + 1 != s_e) {
        pos = min_max_element(s_m + 1, s_e, op.compare_x(), op.compare_y());
        q_t_m = op.update_x_square(q_t_m, *pos.first);
        q_r_m = op.update_y_square(q_r_m, *pos.second);
      }

      // step 3/4
      if (op.compute_x_distance(
          op.place_x_square(q_t_m, r, op.delta()(*s_m)), B) <=
          op.delta()(*s_m) &&
          op.compute_y_distance(
            op.place_y_square(q_r_m, r, op.delta()(*s_m)), B) <=
          op.delta()(*s_m)) {
        s_e = s_m + 1;
        q_t_afap = q_t_m;
        q_r_afap = q_r_m;
      } else
        s_b = s_m + 1;
    }

    // now s_b corresponds to the first moment in [s, m+1)
    // where q_t and q_r cover B

    // place q_t and q_r
    q_t = op.place_x_square(q_t_afap, r, op.delta()(*s_b));
    q_r = op.place_y_square(q_r_afap, r, op.delta()(*s_b));

    // partition the range [s_b+1, e) into ranges
    // [s_b+1, b1), [b1, b2),   [b2, b3) and [b3, e)
    //     R      G cap q_t  G cap q_r      none
    le_delta_m = bind1st(greater_equal< FT >(), op.delta()(*s_b));
    b2 = partition(s_b + 1, e, compose1(le_delta_m,
                                        bind1st(op.distance(), q_t)));
    b1 = partition(s_b + 1, b2, compose1(le_delta_m,
                                         bind1st(op.distance(), q_r)));
    b3 = partition(b2, e, compose1(le_delta_m, bind1st(op.distance(), q_r)));

    if (b3 != e ||
        !Q_t_empty && op.compute_x_distance(q_t, Q_t) > op.delta()(*s_b) ||
        !Q_r_empty && op.compute_y_distance(q_r, Q_r) > op.delta()(*s_b)) {
      // no covering
      CGAL_assertion(b1 - s >= cutoff);
      s = b1;
      rho_min = op.delta()(*s_b);
      q_t_at_rho_min = q_t, q_r_at_rho_min = q_r;
      q_t_q_r_cover_at_rho_min = 0;
      if (!Q_t_empty)
        q_t_q_r_cover_at_rho_min =
          max(q_t_q_r_cover_at_rho_min, op.compute_x_distance(q_t, Q_t));
      if (!Q_r_empty)
        q_t_q_r_cover_at_rho_min =
          max(q_t_q_r_cover_at_rho_min, op.compute_y_distance(q_r, Q_r));
      s_at_rho_min = s, e_at_rho_min = e;
      continue;
    }

    // we still have a covering

    if (s_b == s) {
      CGAL_expensive_assertion_code(
        RandomAccessIterator ii =
        partition(f, l, compose1(le_delta_m, op.delta()));
        pos = min_max_element(ii, l, op.compare_x(), op.compare_y());
        )
      CGAL_expensive_assertion(!op.compare_x()(*pos.first, q_t));
      CGAL_expensive_assertion(!op.compare_y()(*pos.second, q_r));

      // we are done
      rho_max = op.delta()(*s);
      q_t_at_rho_max = q_t, q_r_at_rho_max = q_r;
      radius_is_known = true;
      break;
    }

    // if there are enough points in G ...
    if (b3 - b1 >= cutoff) {
      e = b1;
      // adjust Q_t
      if (b1 != b2)
        if (Q_t_empty) {
          Q_t = bounding_box_2(b1, b2);
          Q_t_empty = false;
        } else
          Q_t = construct_bounding_box_union(Q_t, bounding_box_2(b1, b2));
      // adjust Q_r
      if (b2 != b3)
        if (Q_r_empty) {
          Q_r = bounding_box_2(b2, b3);
          Q_r_empty = false;
        } else
          Q_r = construct_bounding_box_union(Q_r, bounding_box_2(b2, b3));
      continue;
    }

    // we have to take the next smaller radius
    RandomAccessIterator next =
      max_element_if(s, s_b,
                     compose2_2(less< FT >(), op.delta(), op.delta()),
                     compose1(bind1st(not_equal_to< FT >(),
                                      op.delta()(*s_b)),
                     op.delta()));
    rho_max = op.delta()(*s_b);
    q_t_at_rho_max = q_t, q_r_at_rho_max = q_r;
    CGAL_assertion(op.delta()(*next) < op.delta()(*s_b));
    q_t_afap = op.update_x_square(q_t_afap, *s_b);
    q_r_afap = op.update_y_square(q_r_afap, *s_b);
    q_t = op.place_x_square(q_t_afap, r, op.delta()(*next));
    q_r = op.place_y_square(q_r_afap, r, op.delta()(*next));

    // again check for covering
    le_delta_m = bind1st(greater_equal< FT >(), op.delta()(*next));
    b2 = partition(s_b, e,
                   compose1(le_delta_m, bind1st(op.distance(), q_t)));
    b1 = partition(s_b, b2,
                   compose1(le_delta_m, bind1st(op.distance(), q_r)));
    b3 = partition(b2, e,
                   compose1(le_delta_m, bind1st(op.distance(), q_r)));

    if (b3 != e ||
        !Q_t_empty && op.compute_x_distance(q_t, Q_t) > op.delta()(*next) ||
        !Q_r_empty && op.compute_y_distance(q_r, Q_r) > op.delta()(*next)) {
      // no covering
      rho_min = op.delta()(*next);
      q_t_q_r_cover_at_rho_min = 0;
      if (!Q_t_empty)
        q_t_q_r_cover_at_rho_min =
          max(q_t_q_r_cover_at_rho_min, op.compute_x_distance(q_t, Q_t));
      if (!Q_r_empty)
        q_t_q_r_cover_at_rho_min =
          max(q_t_q_r_cover_at_rho_min, op.compute_y_distance(q_r, Q_r));
      q_t_at_rho_min = q_t, q_r_at_rho_min = q_r;
      s_at_rho_min = b3, e_at_rho_min = e;
      radius_is_known = true;
      break;
    }

    // still a covering, but now there must be enough points in G
    CGAL_assertion(b3 - b1 >= cutoff);
    e = b1;
    // adjust Q_t
    if (b1 != b2)
      if (Q_t_empty) {
        Q_t = bounding_box_2(b1, b2);
        Q_t_empty = false;
      } else
        Q_t = construct_bounding_box_union(Q_t, bounding_box_2(b1, b2));
    // adjust Q_r
    if (b2 != b3)
      if (Q_r_empty) {
        Q_r = bounding_box_2(b2, b3);
        Q_r_empty = false;
      } else
        Q_r = construct_bounding_box_union(Q_r, bounding_box_2(b2, b3));

  } // while (e - s > 6)

  // compute the solution brute-force
  if (!radius_is_known) {
    RandomAccessIterator t = e;
    Point q_t_afap = op.place_x_square(Q_t_empty, Q_t, r);
    Point q_r_afap = op.place_y_square(Q_r_empty, Q_r, r);
    if (s != e) {
      sort(s, e, compose2_2(less< FT >(), op.delta(), op.delta()));
      rho_max = op.delta()(*--t);
    } else
      rho_max = rho_min;
    Point q_t = op.place_x_square(q_t_afap, r, rho_max);
    Point q_r = op.place_y_square(q_r_afap, r, rho_max);

#ifdef CGAL_PCENTER_DEBUG
    W << BLACK << Rectangle(r[0], Point(r[0].x() + rho_max,
                                        r[0].y() + rho_max));
    W << GREEN << r << ORANGE;
    if (!Q_t_empty)
      W << Q_t;
    if (!Q_r_empty)
      W << Q_r;
    W.set_node_width(5);
    W << VIOLET << q_t_afap << q_r_afap;
    W.set_node_width(4);
    W << RED;
    std::copy(f, l, wout);
    W.set_node_width(3);
    W << BLACK;
    std::copy(s, e, wout);
    W << *s << BLUE << q_t << q_r;
    { Point p; W << WHITE >> p; }
#endif // CGAL_PCENTER_DEBUG

    if (!Q_t_empty && op.compute_x_distance(q_t, Q_t) > rho_max ||
        !Q_r_empty && op.compute_y_distance(q_r, Q_r) > rho_max) {
      rho_max = max(op.compute_x_distance(q_t, Q_t),
                    op.compute_y_distance(q_r, Q_r));
#ifndef CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
      CGAL_assertion(rho_max <= rad);
#endif // ! CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
      rad = rho_max;
      *o++ = op.construct_corner_square(r, rad / FT(2));
      *o++ = op.construct_x_square(q_t, rad / FT(2));
      *o++ = op.construct_y_square(q_r, rad / FT(2));
      return o;
    }
    CGAL_assertion(s != e);

    // find the first diameter where covering is possible
    for (;;) {
      q_t_at_rho_max = q_t, q_r_at_rho_max = q_r;

      // these get uncovered now
      do {
        q_t_afap = op.update_x_square(q_t_afap, *t);
        q_r_afap = op.update_y_square(q_r_afap, *t);
      } while (t != s && op.delta()(*--t) == rho_max);

      if (op.delta()(*t) == rho_max)
        break;

      // try the next possible diameter value
      FT try_rho = op.delta()(*t);
      CGAL_assertion(try_rho < rho_max);
      q_t = op.place_x_square(q_t_afap, r, try_rho);
      q_r = op.place_y_square(q_r_afap, r, try_rho);

      // check for covering
      binder1st< less< FT > >
        greater_rho_max(bind1st(less< FT >(), try_rho));
      if (!Q_t_empty && op.compute_x_distance(q_t, Q_t) > try_rho ||
          !Q_r_empty && op.compute_y_distance(q_r, Q_r) > try_rho ||
          e != find_if(t + 1,
                       e,
                       compose2(logical_and< bool >(),
                                compose1(greater_rho_max,
                                         bind1st(op.distance(), q_t)),
                                compose1(greater_rho_max,
                                         bind1st(op.distance(), q_r))))) {
        rho_min = try_rho;
        q_t_q_r_cover_at_rho_min = 0;
        if (!Q_t_empty)
          q_t_q_r_cover_at_rho_min =
            max(q_t_q_r_cover_at_rho_min, op.compute_x_distance(q_t, Q_t));
        if (!Q_r_empty)
          q_t_q_r_cover_at_rho_min =
            max(q_t_q_r_cover_at_rho_min, op.compute_y_distance(q_r, Q_r));
        q_t_at_rho_min = q_t, q_r_at_rho_min = q_r;
        s_at_rho_min = t + 1, e_at_rho_min = e;
        break;
      }
      rho_max = try_rho;
    } // for (;;)
  } // if (!radius_is_known)

  // now we have the following:
  // rho_min corresponds to a non-covering with
  //   - q_t_at_rho_min is the corr. position of q_t,
  //   - q_r_at_rho_min is the corr. position of q_r and
  //   - q_t_q_r_cover_at_rho_min is the radius needed to
  //     cover Q_t and Q_r
  //   - the range [s_at_rho_min, e_at_rho_min) contains the points
  //     still to be covered
  // rho_min corresponds to a covering with
  //   - q_t_at_rho_min is the corr. position of q_t and
  //   - q_r_at_rho_min is the corr. position of q_r.

  // try rho_min
  CGAL_assertion(rho_min <= rho_max);
  FT rad2 = q_t_q_r_cover_at_rho_min;
  typedef binder1st< Distance > Dist_bind;
  if (s_at_rho_min != e_at_rho_min) {
    binary_compose< Min< FT >, Dist_bind, Dist_bind >
      mydist(compose2(Min< FT >(),
                      bind1st(op.distance(), q_t_at_rho_min),
                      bind1st(op.distance(), q_r_at_rho_min)));
    rad2 =
      max(rad2,
          mydist(*max_element(s_at_rho_min,
                              e_at_rho_min,
                              compose2_2(less< FT >(), mydist, mydist))));
  }
  CGAL_assertion(rad2 == 0 || rad2 > rho_min);

  // if a covering with rho == 0 is possible,
  // it will be catched in the type1 functions
  Point q_t, q_r;
  if (rad2 > rho_max || rho_min == 0) {
    // it is rho_max ...
    q_t = q_t_at_rho_max, q_r = q_r_at_rho_max;
    rad2 = rho_max;
  } else
    q_t = q_t_at_rho_min, q_r = q_r_at_rho_min;

#ifndef CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
  CGAL_assertion(rad2 <= rad);
#endif // ! CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
  rad = rad2;
  *o++ = op.construct_corner_square(r, rad / FT(2));
  *o++ = op.construct_x_square(q_t, rad / FT(2));
  *o++ = op.construct_y_square(q_r, rad / FT(2));
  return o;
} // rectangular_3_center_2_type2( ... )
template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
rectangular_3_center_2(
  ForwardIterator f,
  ForwardIterator l,
  OutputIterator o,
  typename Traits::FT& r,
  Traits& t)
{
  CGAL_precondition(f != l);
  typedef typename Traits::FT                                    FT;
  typedef typename Traits::Point_2                            Point;
  typedef typename Traits::Iso_rectangle_2                Rectangle;
  typedef Rectangular_3_center_2_type2_operations0< Traits >    Op0;
  typedef Rectangular_3_center_2_type2_operations1< Traits >    Op1;
  typedef Rectangular_3_center_2_type2_operations2< Traits >    Op2;
  typedef Rectangular_3_center_2_type2_operations3< Traits >    Op3;

  std::vector< Point > points(f, l);
  Rectangle bb = bounding_box_2(points.begin(), points.end());

  // try to place two squares at opposite corners of bb
  Point ptst[3];
  rectangular_3_center_2_type1(
    points.begin(), points.end(), bb, ptst, r, t);

  // try to place one square at a corner and the others
  // at the two remaining sides of bb
  Point pts0[3];
  Point pts1[3];
  Point pts2[3];
  Point pts3[3];
  Point* pts = ptst;
  FT rmin = r;

#ifdef CGAL_PCENTER_DEBUG
  std::cerr << "1: rmin = " << rmin << std::endl;
#endif // CGAL_PCENTER_DEBUG

  rectangular_3_center_2_type2(
    points.begin(), points.end(), bb, pts0, r, Op0(t, bb[0]));
  if (r < rmin)
    pts = pts0, rmin = r;
#ifdef CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
  else
    r = rmin;
#endif // CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST

#ifdef CGAL_PCENTER_DEBUG
  std::cerr << "2: rmin = " << rmin << std::endl;
#endif // CGAL_PCENTER_DEBUG

  rectangular_3_center_2_type2(
    points.begin(), points.end(), bb, pts1, r, Op1(t, bb[1]));
  if (r < rmin)
    pts = pts1, rmin = r;
#ifdef CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
  else
    r = rmin;
#endif // CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST

#ifdef CGAL_PCENTER_DEBUG
  std::cerr << "3: rmin = " << rmin << std::endl;
#endif // CGAL_PCENTER_DEBUG

  rectangular_3_center_2_type2(
    points.begin(), points.end(), bb, pts2, r, Op2(t, bb[2]));
  if (r < rmin)
    pts = pts2, rmin = r;
#ifdef CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
  else
    r = rmin;
#endif // CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST

#ifdef CGAL_PCENTER_DEBUG
  std::cerr << "4: rmin = " << rmin << std::endl;
#endif // CGAL_PCENTER_DEBUG

  rectangular_3_center_2_type2(
    points.begin(), points.end(), bb, pts3, r, Op3(t, bb[3]));
  if (r < rmin)
    pts = pts3;
#ifdef CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST
  else
    r = rmin;
#endif // CGAL_3COVER_NO_CHECK_OPTIMUM_FIRST

  *o++ = pts[0];
  *o++ = pts[1];
  *o++ = pts[2];
  return o;

} // rectangular_3_center_2( ... )

CGAL_END_NAMESPACE

#ifdef _MSC_VER
#undef Base
#endif // _MSC_VER

#endif // ! (CGAL_RECTANGULAR_3_CENTER_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

