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
// file          : pierce_rectangles_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2-4-Piercing Axis-Parallel 2D-Rectangles
// ============================================================================

#if ! (CGAL_PIERCE_RECTANGLES_2_H)
#define CGAL_PIERCE_RECTANGLES_2_H 1

#include <CGAL/optimisation_assertions.h>
#include <CGAL/circulator.h>
#include <CGAL/function_objects.h>
#include <CGAL/algorithm.h>
#include <algorithm>
#include <vector>

CGAL_BEGIN_NAMESPACE

//!!! STL-extensions
template < class T >
struct Wastebasket : public CGAL_STD::output_iterator
{
  typedef Wastebasket< T > iterator;

  iterator
  operator=( const T&)
  { return *this; }

  iterator
  operator*()
  { return *this; }

  iterator
  operator++()
  { return *this; }

  iterator
  operator++( int)
  { return *this; }
};

template < class _Traits >
struct Loc_domain {
  // ---------------------------------------------
  // types:

  typedef _Traits                               Traits;
  typedef typename Traits::FT                   FT;
  typedef typename Traits::Point_2              Point_2;

  typedef std::vector< Point_2 >                Container;
  typedef typename Container::iterator          Iterator;
  typedef typename Container::const_iterator    Citerator;
  typedef typename Container::reverse_iterator  Riterator;
  typedef typename Container::const_reference   Creference;
  typedef typename Container::size_type         size_type;

  // ---------------------------------------------
  // creation:

  void
  update(int j, Citerator i)
  {
    CGAL_optimisation_precondition(j >= 0 && j < 4);
    if (j < 2)
      if (j == 0) {
        if (traits.get_less_x_2()(*i, minx)) minx = *i;
        if (traits.get_less_y_2()(*i, miny)) miny = *i;
      }
      else {
        if (traits.get_less_y_2()(*i, miny)) miny = *i;
        if (traits.get_less_x_2()(maxx, *i)) maxx = *i;
      }
    else
      if (j == 2) {
        if (traits.get_less_x_2()(maxx, *i)) maxx = *i;
        if (traits.get_less_y_2()(maxy, *i)) maxy = *i;
      }
      else {
        if (traits.get_less_y_2()(maxy, *i)) maxy = *i;
        if (traits.get_less_x_2()(*i, minx)) minx = *i;
      }
  }

  template < class InputIC >
  Loc_domain(InputIC b, InputIC e, Traits t)
  : pts(b, e),
    end_(pts.end()),
    minx(pts.front()),
    miny(pts.front()),
    maxx(pts.front()),
    maxy(pts.front()),
    traits(t)
  {
    CGAL_optimisation_precondition(b != e);
    Iterator i = pts.begin();
    CGAL_optimisation_assertion(i != pts.end());
    while (++i != pts.end()) {
      if (traits.get_less_x_2()(*i, minx)) minx = *i;
      if (traits.get_less_x_2()(maxx, *i)) maxx = *i;
      if (traits.get_less_y_2()(*i, miny)) miny = *i;
      if (traits.get_less_y_2()(maxy, *i)) maxy = *i;
    }
  }

  // ---------------------------------------------
  // access operations:

  Point_2
  operator[](int i) const
  // return corner points (0 <-> bottom-left, 1 <-> bottom-right)
  {
    CGAL_optimisation_precondition(i >= 0 && i < 4);
    if (i == 0)
      return traits.get_construct_point_2_above_right_implicit_point_2()(
        minx, miny, r);
    else if (i == 1)
      return traits.get_construct_point_2_above_left_implicit_point_2()(
        maxx, miny, r);
    else if (i == 2)
      return traits.get_construct_point_2_below_left_implicit_point_2()(
        maxx, maxy, r);
    return traits.get_construct_point_2_below_right_implicit_point_2()(
      minx, maxy, r);
  }

  Point_2
  extreme(int i) const
  // return extreme points (0 <-> left, 1 <-> bottom)
  {
    CGAL_optimisation_precondition(i >= 0 && i < 4);
    if (i > 1) return i == 2 ? maxx : maxy;
    return i == 0 ? minx : miny;
  }

  Point_2&
  extreme(int i)
  // return extreme points (0 <-> left, 1 <-> bottom)
  {
    CGAL_optimisation_precondition(i >= 0 && i < 4);
    if (i > 1) return i == 2 ? maxx : maxy;
    return i == 0 ? minx : miny;
  }

  Citerator begin() const { return pts.begin(); }
  Iterator begin()        { return pts.begin(); }
#ifndef _MSC_VER
  Citerator end() const   { return end_; }
#else
  // Yet another really great MSVC feature ...
  // without that static_cast to itself it does not compile :-)
  Citerator end() const   { return static_cast<Iterator>(end_); }
#endif
  Iterator& end()         { return end_; }
  Iterator real_end()     { return pts.end(); }

  bool empty() const        { return begin() == end(); }
  size_type size() const    { return end() - begin(); }
  Creference front() const  { return pts.front(); }

  // ---------------------------------------------
  // check operation:

  void
  check() const {
    CGAL_optimisation_expensive_assertion_code(
      Iterator i = pts.begin();
      do {
        CGAL_optimisation_assertion(!traits.get_less_x_2()(*i, minx));
        CGAL_optimisation_assertion(!traits.get_less_x_2()(maxx, *i));
        CGAL_optimisation_assertion(!traits.get_less_y_2()(*i, miny));
        CGAL_optimisation_assertion(!traits.get_less_y_2()(maxy, *i));
      } while (++i != end);
      )
  }

protected:
  // points container
  Container pts;
  // const iterator to past-the-end
  Iterator end_;
public:
  // actual center radius
  FT r;
  // (copies of) elements with minimal/maximal x/y coordinate
  Point_2 minx, miny, maxx, maxy;
  // Traits class
  Traits traits;

}; // class Loc_domain
template < class _Traits >
struct Staircases : public Loc_domain< _Traits > {
  typedef _Traits                           Traits;
  typedef Loc_domain< Traits >              Base;
  typedef typename Base::Container          Container;
  typedef typename Base::Iterator           Iterator;
  typedef typename Base::Citerator          Citerator;
  typedef typename Base::Riterator          Riterator;
  typedef typename Traits::Point_2          Point_2;
  typedef typename Traits::FT               FT;
  typedef std::pair< Point_2, Point_2 >     Intervall;

  template < class InputIC >
  Staircases(InputIC b, InputIC e, Traits t)
  : Base(b, e, t)
  {
#ifndef CGAL_CFG_NO_NAMESPACE
    using std::sort;
    using std::find_if;
    using std::bind1st;
    using std::bind2nd;
#endif // ! CGAL_CFG_NO_NAMESPACE

    // build top-left and bottom-right staircases
    sort(pts.begin(), pts.end(), traits.get_less_y_2());
    // bottom-right
    Iterator i = pts.begin();
    do {
      brstc.push_back(*i++);
      i = find_if(i, pts.end(),
                  bind1st(traits.get_less_x_2(), brstc.back()));
    } while (i != pts.end());
    // top-left
    Riterator j = pts.rbegin();
    do {
      tlstc.push_back(*j++);
      j = find_if(j, pts.rend(),
                  bind2nd(traits.get_less_x_2(), tlstc.back()));
    } while (j != pts.rend());

    // build left-bottom and right-top staircases
    sort(pts.begin(), pts.end(), traits.get_less_x_2());
    // left-bottom
    i = pts.begin();
    do {
      lbstc.push_back(*i++);
      i = find_if(i, pts.end(),
                  bind2nd(traits.get_less_y_2(), lbstc.back()));
    } while (i != pts.end());
    // right-top
    j = pts.rbegin();
    do {
      rtstc.push_back(*j++);
      j = find_if(j, pts.rend(),
                  bind1st(traits.get_less_y_2(), rtstc.back()));
    } while (j != pts.rend());
  } // Staircases(b, e, t)

  bool is_middle_empty() const {
    //!!! the "middle" point could be precomputed in advance
    Citerator i = pts.begin();
    do
      if (traits.get_signed_x_distance_2()(maxx, *i) > FT(2) * r &&
          traits.get_signed_x_distance_2()(*i, minx) > FT(2) * r &&
          traits.get_signed_y_distance_2()(*i, miny) > FT(2) * r &&
          traits.get_signed_y_distance_2()(maxy, *i) > FT(2) * r)
        return false;
    while (++i != pts.end());
    return true;
  } // is_middle()

  Citerator tlstc_begin() const { return tlstc.begin(); }
  Citerator tlstc_end() const   { return tlstc.end(); }
  Citerator lbstc_begin() const { return lbstc.begin(); }
  Citerator lbstc_end() const   { return lbstc.end(); }
  Citerator brstc_begin() const { return brstc.begin(); }
  Citerator brstc_end() const   { return brstc.end(); }
  Citerator rtstc_begin() const { return rtstc.begin(); }
  Citerator rtstc_end() const   { return rtstc.end(); }

  Intervall top_intervall() const {
    Point_2 p =
      traits.get_construct_point_2_above_right_implicit_point_2()(
        minx, miny, FT(2) * r);
    Point_2 q =
      traits.get_construct_point_2_above_left_implicit_point_2()(
        maxx, miny, FT(2) * r);

    Citerator i =
      min_element_if(
        pts.begin(), pts.end(),
        traits.get_less_x_2(),
        std::compose2(std::logical_and< bool >(),
                      std::bind1st(traits.get_less_x_2(), p),
                      std::bind1st(traits.get_less_y_2(), p)));
    Citerator j =
      max_element_if(
        pts.begin(), pts.end(),
        traits.get_less_x_2(),
        std::compose2(std::logical_and< bool >(),
                      std::bind2nd(traits.get_less_x_2(), q),
                      std::bind1st(traits.get_less_y_2(), q)));
    return Intervall(i == pts.end() ? maxx : *i,
                     j == pts.end() ? minx : *j);
  } // top_intervall()

  Intervall bottom_intervall() const {
    Point_2 p =
      traits.get_construct_point_2_below_right_implicit_point_2()(
        minx, maxy, FT(2) * r);
    Point_2 q =
      traits.get_construct_point_2_below_left_implicit_point_2()(
        maxx, maxy, FT(2) * r);

    Citerator i =
      min_element_if(
        pts.begin(), pts.end(),
        traits.get_less_x_2(),
        std::compose2(std::logical_and< bool >(),
                      std::bind1st(traits.get_less_x_2(), p),
                      std::bind2nd(traits.get_less_y_2(), p)));
    Citerator j =
      max_element_if(
        pts.begin(), pts.end(),
        traits.get_less_x_2(),
        std::compose2(std::logical_and< bool >(),
                      std::bind2nd(traits.get_less_x_2(), q),
                      std::bind2nd(traits.get_less_y_2(), q)));
    return Intervall(i == pts.end() ? maxx : *i,
                     j == pts.end() ? minx : *j);
  } // bottom_intervall()

  Intervall left_intervall() const {
    Point_2 p =
      traits.get_construct_point_2_above_left_implicit_point_2()(
        maxx, miny, FT(2) * r);
    Point_2 q =
      traits.get_construct_point_2_below_left_implicit_point_2()(
        maxx, maxy, FT(2) * r);

    Citerator i =
      min_element_if(
        pts.begin(), pts.end(),
        traits.get_less_y_2(),
        std::compose2(std::logical_and< bool >(),
                      std::bind2nd(traits.get_less_x_2(), p),
                      std::bind1st(traits.get_less_y_2(), p)));
    Citerator j =
      max_element_if(
        pts.begin(), pts.end(),
        traits.get_less_y_2(),
        std::compose2(std::logical_and< bool >(),
                      std::bind2nd(traits.get_less_x_2(), q),
                      std::bind2nd(traits.get_less_y_2(), q)));
    return Intervall(i == pts.end() ? maxy : *i,
                     j == pts.end() ? miny : *j);
  } // left_intervall()

  Intervall right_intervall() const {
    Point_2 p =
      traits.get_construct_point_2_above_right_implicit_point_2()(
        minx, miny, FT(2) * r);
    Point_2 q =
      traits.get_construct_point_2_below_right_implicit_point_2()(
        minx, maxy, FT(2) * r);

    Citerator i =
      min_element_if(
        pts.begin(), pts.end(),
        traits.get_less_y_2(),
        std::compose2(std::logical_and< bool >(),
                      std::bind1st(traits.get_less_x_2(), p),
                      std::bind1st(traits.get_less_y_2(), p)));
    Citerator j =
      max_element_if(
        pts.begin(), pts.end(),
        traits.get_less_y_2(),
        std::compose2(std::logical_and< bool >(),
                      std::bind1st(traits.get_less_x_2(), q),
                      std::bind2nd(traits.get_less_y_2(), q)));
    return Intervall(i == pts.end() ? maxy : *i,
                     j == pts.end() ? miny : *j);
  } // right_intervall()

private:
  Container tlstc, lbstc, brstc, rtstc;
};

CGAL_END_NAMESPACE
//#ifdef CGAL_REP_CLASS_DEFINED
//#include <CGAL/Pierce_rectangles_2_traits.h>
//#endif // CGAL_REP_CLASS_DEFINED
CGAL_BEGIN_NAMESPACE

template < class InputIC, class OutputIterator, class Traits >
inline OutputIterator
two_cover_points(
  InputIC f, InputIC l, OutputIterator o, bool& ok, const Traits& t)
{
  CGAL_optimisation_precondition(f != l);

  // compute location domain:
  Loc_domain< Traits > d(f, l, t);

  return two_cover_points(d, o, ok, t);
} // two_cover_points(f, l, o, ok, t)
template < class InputIC, class OutputIterator, class Traits >
inline OutputIterator
three_cover_points(
  InputIC f, InputIC l, OutputIterator o, bool& ok, const Traits& t)
{
  CGAL_optimisation_precondition(f != l);

  // compute location domain:
  Loc_domain< Traits > d(f, l, t);

  return three_cover_points(d, o, ok, t);
} // three_cover_points(f, l, o, ok, t)
template < class InputIC, class OutputIterator, class Traits >
inline OutputIterator
four_cover_points(
  InputIC f, InputIC l, OutputIterator o, bool& ok, const Traits& t)
{
  CGAL_optimisation_precondition(f != l);

  // compute location domain:
  Loc_domain< Traits > d(f, l, t);

  return four_cover_points(d, o, ok, t);
} // four_cover_points(f, l, o, ok, t)
template < class OutputIterator, class Traits >
OutputIterator
two_cover_points(
  const Loc_domain< Traits >& d,
  OutputIterator o,
  bool& ok)
{
#ifndef CGAL_CFG_NO_NAMESPACE
  using std::compose1;
  using std::compose2;
  using std::bind1st;
  using std::find_if;
  using std::less;
#endif

  typedef typename Traits::FT           FT;
  typedef typename Traits::Point_2      Point_2;
  typename Traits::Infinity_distance_2 dist =
    d.traits.get_infinity_distance_2();
  typename Traits::Signed_infinity_distance_2 sdist =
    d.traits.get_signed_infinity_distance_2();

  Min< FT > minft;
  less< FT > lessft;

  if (sdist(d[2], d[0]) <= FT(0)) {
    // the location domain is degenerate and [f,l) is one-pierceable
    *o++ = d[0];
    ok = true;
    return o;
  }
  // check whether {d[0], d[2]} forms a piercing set
  if (d.end() ==
      find_if(d.begin(),
              d.end(),
              compose1(
                bind1st(lessft, d.r),
                compose2(
                  minft, bind1st(dist, d[0]), bind1st(dist, d[2])))))
  {
    *o++ = d[0];
    *o++ = d[2];
    ok = true;
    return o;
  }
  // check whether {d[1], d[3]} forms a piercing set
  if (d.end() ==
      find_if(d.begin(),
              d.end(),
              compose1(
                bind1st(lessft, d.r),
                compose2(
                  minft, bind1st(dist, d[1]), bind1st(dist, d[3])))))
  {
    *o++ = d[1];
    *o++ = d[3];
    ok = true;
    return o;
  }

  // no piercing set exists:
  ok = false;
  return o;
} // two_cover_points(d, o, ok, t)
template < class OutputIterator, class Traits >
OutputIterator
three_cover_points(
  Loc_domain< Traits >& d,
  OutputIterator o,
  bool& ok)
{
#ifndef CGAL_CFG_NO_NAMESPACE
  using std::find_if;
  using std::bind1st;
  using std::compose1;
  using std::less;
  using std::iter_swap;
#endif
  CGAL_optimisation_precondition(!d.empty());

  // typedefs:
  typedef typename Traits::FT                      FT;
  typedef typename Traits::Point_2                 Point_2;
  typedef typename Loc_domain< Traits >::Iterator  Iterator;
  typename Traits::Infinity_distance_2 dist =
    d.traits.get_infinity_distance_2();
  less< FT > lessft;

  // test the four corners:
  for (int k = 0; k < 4; ++k) {
    
    // extract all points which are close enough to d[k]
    Point_2 corner = d[k];
    
    // find first point not covered by the rectangle at d[k]
    Iterator i = find_if(d.begin(), d.end(),
                         compose1(bind1st(lessft, d.r),
                                  bind1st(dist, corner)));
    
    // are all points already covered?
    if (i == d.end()) {
      CGAL_optimisation_assertion(k == 0);
      *o++ = d[0];
      ok = true;
      return o;
    } // if (i == d.end())
    
    // save changing sides of d:
    Point_2 save_side1 = d.extreme(k);
    Point_2 save_side2 = d.extreme((k+1) % 4);
    Iterator save_end = d.end();
    
    // now run through it:
    // initialize the two (possibly) changing sides of d
    d.extreme(k) = d.extreme((k+1) % 4) = *i;
    
    // is there any point covered?
    if (i == d.begin()) {
      while (++i != d.end() && dist(corner, *i) > d.r)
        d.update(k, i);
      d.end() = i;
    } else {
      // move *i to the begin of pts
      iter_swap(i, d.begin());
      d.end() = d.begin() + 1;
    }
    
    // [d.begin(), d.end()) shall be the range of uncovered points
    if (i != save_end)
      while (++i != save_end)
        if (dist(corner, *i) > d.r) {
          d.update(k, i);
          iter_swap(i, d.end());
          ++d.end();
        } // if (dist(corner, *i) > d.r)
    
    // check disjoint for two-pierceability:
    
    CGAL_optimisation_assertion(
      save_end == find_if(d.end(), save_end,
                          compose1(bind1st(lessft, d.r),
                                   bind1st(dist, corner))));
    CGAL_optimisation_assertion(
      d.end() == find_if(d.begin(), d.end(),
                         compose1(bind1st(std::greater_equal<FT>(), d.r),
                                  bind1st(dist, corner))));
    
    #if 0
    if (false && d.r <= 0.27710729236714543023 &&
        d.r >= 0.27710729236714543021)
      {
        typedef Ostream_iterator< Point_2, leda_window > OIP;
        leda_window W(730, 690);
        cgalize(W);
        W.set_node_width(2);
        W.init(-1.5, 1.5, -1.2);
        W.display();
        OIP oip(W);
        typedef typename Traits::Iso_rectangle_2 Iso_rectangle_2;
    
        FT rsav = d.r;
        d.r = 0;
        Point_2 ss1 = d.extreme(k);
        Point_2 ss2 = d.extreme((k+1) % 4);
        d.extreme(k) = save_side1;
        d.extreme((k+1) % 4) = save_side2;
        W << GREEN << Iso_rectangle_2(d[0], d[2]);
        Point_2 psav = d[k];
        d.r = FT(2) * rsav;
        W << ORANGE << Iso_rectangle_2(psav, d[k]);
        d.extreme(k) = ss1;
        d.extreme((k+1) % 4) = ss2;
        d.r = 0;
        W << VIOLET << Iso_rectangle_2(d[0], d[2]);
        d.r = rsav;
    
        W << RED; std::copy(d.begin(), d.end(), oip);
        W << BLUE; std::copy(d.end(), save_end, oip);
        W << BLACK; std::copy(save_end, d.real_end(), oip);
    
        { Point_2 du; W >> du; }
      }
    #endif // 0
    
    two_cover_points(d, o, ok);
    
    // restore saved sides of d:
    d.extreme(k) = save_side1;
    d.extreme((k+1) % 4) = save_side2;
    
    if (ok) {
      // does any rectangle contain the corner?
      if (d.end() != save_end) {
        *o++ = d[k];
        d.end() = save_end;
      }
      return o;
    } // if (ok)
    
    d.end() = save_end;
  } // for (int k = 0; k < 4; ++k)

  // no piercing set exists:
  ok = false;
  return o;

} // three_cover_points(d, o, ok)
CGAL_END_NAMESPACE
CGAL_BEGIN_NAMESPACE
template < class OutputIterator, class Traits >
OutputIterator
four_cover_points(Staircases< Traits >& d, OutputIterator o, bool& ok)
{

  #ifndef CGAL_CFG_NO_NAMESPACE
  using std::less;
  using std::iter_swap;
  using std::find_if;
  using std::bind1st;
  using std::bind2nd;
  using std::compose1;
  #endif
  
  typedef typename Traits::Point_2                  Point_2;
  typedef typename Traits::FT                       FT;
  typedef typename Traits::Less_x_2                 Less_x_2;
  typedef typename Traits::Less_y_2                 Less_y_2;
  typedef typename Traits::Signed_x_distance_2      Signed_x_distance_2;
  typedef typename Traits::Signed_y_distance_2      Signed_y_distance_2;
  typedef typename Staircases< Traits >::Iterator   Iterator;
  typedef typename Staircases< Traits >::Citerator  Citerator;
  typedef typename Staircases< Traits >::Intervall  Intervall;
  
  less< FT > lessft;
  typename Traits::Infinity_distance_2 dist =
    d.traits.get_infinity_distance_2();
  
  

  // test the four corners:
  for (int k = 0; k < 4; ++k) {
    
    // extract all points which are close enough to d[k]
    Point_2 corner = d[k];
    
    // find first point not covered by the rectangle at d[k]
    Iterator i = find_if(d.begin(), d.end(),
                         compose1(bind1st(lessft, d.r),
                                  bind1st(dist, corner)));
    
    // are all points already covered?
    if (i == d.end()) {
      CGAL_optimisation_assertion(k == 0);
      *o++ = d[0];
      ok = true;
      return o;
    } // if (i == d.end())
    
    // save changing sides of d:
    Point_2 save_side1 = d.extreme(k);
    Point_2 save_side2 = d.extreme((k+1) % 4);
    Iterator save_end = d.end();
    
    // now run through it:
    // initialize the two (possibly) changing sides of d
    d.extreme(k) = d.extreme((k+1) % 4) = *i;
    
    // is there any point covered?
    if (i == d.begin()) {
      while (++i != d.end() && dist(corner, *i) > d.r)
        d.update(k, i);
      d.end() = i;
    } else {
      // move *i to the begin of pts
      iter_swap(i, d.begin());
      d.end() = d.begin() + 1;
    }
    
    // [d.begin(), d.end()) shall be the range of uncovered points
    if (i != save_end)
      while (++i != save_end)
        if (dist(corner, *i) > d.r) {
          d.update(k, i);
          iter_swap(i, d.end());
          ++d.end();
        } // if (dist(corner, *i) > d.r)
    
    // check disjoint for two-pierceability:
    
    CGAL_optimisation_assertion(
      save_end == find_if(d.end(), save_end,
                          compose1(bind1st(lessft, d.r),
                                   bind1st(dist, corner))));
    CGAL_optimisation_assertion(
      d.end() == find_if(d.begin(), d.end(),
                         compose1(bind1st(std::greater_equal<FT>(), d.r),
                                  bind1st(dist, corner))));
    
    #if 0
    if (false && d.r <= 0.27710729236714543023 &&
        d.r >= 0.27710729236714543021)
      {
        typedef Ostream_iterator< Point_2, leda_window > OIP;
        leda_window W(730, 690);
        cgalize(W);
        W.set_node_width(2);
        W.init(-1.5, 1.5, -1.2);
        W.display();
        OIP oip(W);
        typedef typename Traits::Iso_rectangle_2 Iso_rectangle_2;
    
        FT rsav = d.r;
        d.r = 0;
        Point_2 ss1 = d.extreme(k);
        Point_2 ss2 = d.extreme((k+1) % 4);
        d.extreme(k) = save_side1;
        d.extreme((k+1) % 4) = save_side2;
        W << GREEN << Iso_rectangle_2(d[0], d[2]);
        Point_2 psav = d[k];
        d.r = FT(2) * rsav;
        W << ORANGE << Iso_rectangle_2(psav, d[k]);
        d.extreme(k) = ss1;
        d.extreme((k+1) % 4) = ss2;
        d.r = 0;
        W << VIOLET << Iso_rectangle_2(d[0], d[2]);
        d.r = rsav;
    
        W << RED; std::copy(d.begin(), d.end(), oip);
        W << BLUE; std::copy(d.end(), save_end, oip);
        W << BLACK; std::copy(save_end, d.real_end(), oip);
    
        { Point_2 du; W >> du; }
      }
    #endif // 0
    
    three_cover_points(d, o, ok);
    
    // restore saved sides of d:
    d.extreme(k) = save_side1;
    d.extreme((k+1) % 4) = save_side2;
    
    if (ok) {
      // does any rectangle contain the corner?
      if (d.end() != save_end) {
        *o++ = d[k];
        d.end() = save_end;
      }
      return o;
    } // if (ok)
    
    d.end() = save_end;
  } // for (int k = 0; k < 4; ++k)
  
  // test if four covering rectangles can be placed
  // on the boundary of d, one on each side
  
  // if there is any point that cannot be covered in this way, stop here
  if (d.is_middle_empty()) {
  
    // now try to position the bottom piercing point in each
    // of the intervalls formed by S_bt and S_br
    // (no need to consider S_bl, since we move from left
    // to right and leaving a rectangle won't make piercing easier)
  
    
    Intervall top_i    = d.top_intervall();
    Intervall left_i   = d.left_intervall();
    Intervall bottom_i = d.bottom_intervall();
    Intervall right_i  = d.right_intervall();
    
    Citerator tl = d.tlstc_begin();
    Citerator lb = d.lbstc_begin();
    Citerator br = d.brstc_begin();
    Citerator rt = d.rtstc_begin();
    
    Less_x_2 lessx = d.traits.get_less_x_2();
    Less_y_2 lessy = d.traits.get_less_y_2();
    Signed_x_distance_2 sdistx = d.traits.get_signed_x_distance_2();
    Signed_y_distance_2 sdisty = d.traits.get_signed_y_distance_2();
    
    typename Traits::Construct_projection_onto_horizontal_implicit_line_2
    cpohil =
      d.traits.get_construct_projection_onto_horizontal_implicit_line_2();
    typename Traits::Construct_point_2_above_right_implicit_point_2
    cparip =
      d.traits.get_construct_point_2_above_right_implicit_point_2();
    typename Traits::Construct_point_2_above_left_implicit_point_2
    cpalip =
      d.traits.get_construct_point_2_above_left_implicit_point_2();
    typename Traits::Construct_point_2_below_right_implicit_point_2
    cpbrip =
      d.traits.get_construct_point_2_below_right_implicit_point_2();
    
    // make sure the top intervall is covered (left endpoint)
    // (it might be that top_i.first determines the placement of
    //  the top square)
    Point_2 top = top_i.first;
    if (lessx(top, *tl))
      for (;;) {
        if (++tl != d.tlstc_end() || lessx(*tl, top))
          break;
      }
    else
      top = *tl++;
    
    
    if (tl != d.tlstc_end()) {
      for (;;) {
    
        // make sure the top intervall is covered (right endpoint)
        if (sdistx(top_i.second, top) > FT(2) * d.r)
          break;
    
        // compute position of left square
        Point_2 left = lessy(left_i.second, *tl) ? *tl : left_i.second;
    
        // make sure the left intervall is covered
        if (sdisty(left, left_i.first) <= FT(2) * d.r) {
    
          // compute position of bottom square
          while (lb != d.lbstc_end() && sdisty(left, *lb) <= FT(2) * d.r)
            ++lb;
          // has the left square reached the bottom-left corner?
          if (lb == d.lbstc_end())
            break;
          Point_2 bottom = lessx(bottom_i.first, *lb) ? bottom_i.first : *lb;
    
          // make sure the bottom intervall is covered
          if (sdistx(bottom_i.second, bottom) <= FT(2) * d.r) {
    
            // compute position of right square
            while (br != d.brstc_end() && sdistx(*br, bottom) <= FT(2) * d.r)
              ++br;
            // has the bottom square reached the bottom-right corner?
            if (br == d.brstc_end())
              break;
            Point_2 right = lessy(right_i.first, *br) ? right_i.first : *br;
    
            // make sure the right intervall is covered
            if (sdisty(right_i.second, right) <= FT(2) * d.r) {
    
              // compute right bound for top square
              while (rt != d.rtstc_end() && sdisty(*rt, right) <= FT(2) * d.r)
                ++rt;
              // has the right square reached the top-right corner?
              if (rt == d.rtstc_end())
                break;
    
              // Finally: Do we have a covering?
              if (sdistx(*rt, top) <= FT(2) * d.r) {
                *o++ = cpbrip(d.minx, left, d.r);
                *o++ = cparip(bottom, d.miny, d.r);
                *o++ = cpalip(d.maxx, right, d.r);
                *o++ = cpbrip(top, d.maxy, d.r);
                ok = true;
    
    
                return o;
              } // if (covering)
    
            } // if (sdisty(right_i.second, right) <= FT(2) * d.r)
    
          } // if (sdistx(bottom_i.second, bottom) <= FT(2) * d.r)
    
        } // if (sdisty(left, left_i.first) <= FT(2) * d.r)
    
        top = *tl;
        if (lessy(left_i.second, *tl) || ++tl == d.tlstc_end() ||
            lessx(bottom_i.first, *lb) || lessy(right_i.first, *br))
          break;
    
      } // for (;;)
    } // if (tl != d.tlstc_end())
    
  
  } // if (!d.is_middle_empty())

  ok = false;
  return o;

} // four_cover_points(d, o, ok)

struct Two_covering_algorithm {
  template < class Traits, class OutputIterator >
  OutputIterator
  operator()(Staircases< Traits >& d,
             OutputIterator o,
             bool& ok) const
  { return two_cover_points(d, o, ok); }
}; // class Two_covering_algorithm
struct Three_covering_algorithm {
  template < class Traits, class OutputIterator >
  OutputIterator
  operator()(Staircases< Traits >& d,
             OutputIterator o,
             bool& ok) const
  { return three_cover_points(d, o, ok); }
}; // class Three_covering_algorithm
struct Four_covering_algorithm {
  template < class Traits, class OutputIterator >
  OutputIterator
  operator()(Staircases< Traits >& d,
             OutputIterator o,
             bool& ok) const
  { return four_cover_points(d, o, ok); }
}; // class Four_covering_algorithm
CGAL_END_NAMESPACE

#endif // ! (CGAL_PIERCE_RECTANGLES_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

