// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : min_quadrilateral_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Min_quadrilaterals $
// source        : oops.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch> and
//                 Emo Welzl <emo@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Computing minimum enclosing quadrilaterals of a convex point set
// ============================================================================

#if ! (CGAL_MIN_QUADRILATERAL_2_H)
#define CGAL_MIN_QUADRILATERAL_2_H 1

#include <CGAL/basic.h>
#include <CGAL/Optimisation/assertions.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE

template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
convex_bounding_box_2(
  ForwardIterator f, ForwardIterator l, OutputIterator o, Traits& t)
// PRE:
//   * f != l
//   * value type of ForwardIterator is Traits::Point_2
//   * [f,l) form a the vertices of a convex polygon
//     oriented counterclockwise
//   * OutputIterator accepts ForwardIterator as value type
// POST:
//   writes to o iterators from [f,l) referring to the last points with
//    - smallest x coordinate
//    - smallest y coordinate
//    - largest x coordinate
//    - largest y coordinate
//   in that order.
{
  CGAL_precondition(f != l);

  // make sure that we have two distinct points, such that it
  // can be determined in which quadrant of the polygon we are
  ForwardIterator first;
  do {
    first = f;
    // catch the one-element case:
    if (++f == l) {
      f = first;
      break;
    }
  } while (t.equal_2_object()(*first, *f));

  // Four extremes
  ForwardIterator minx = first;
  ForwardIterator maxx = first;
  ForwardIterator miny = first;
  ForwardIterator maxy = first;

  if (t.less_x_2_object()(*minx, *f) ||
      t.less_y_2_object()(*minx, *f) && !t.less_x_2_object()(*f, *minx))
    if (t.less_y_2_object()(*minx, *f))
      // first quadrant
      for (;;) {
        maxx = f;
        if (++f == l) {
          maxy = minx = miny = maxx;
          break;
        }
        if (t.less_x_2_object()(*f, *maxx)) {
          f = maxx;
          for (;;) {
            maxy = f;
            if (++f == l) {
              minx = miny = maxy;
              break;
            }
            if (t.less_y_2_object()(*f, *maxy)) {
              f = maxy;
              for (;;) {
                minx = f;
                if (++f == l) {
                  miny = minx;
                  break;
                }
                if (t.greater_x_2_object()(*f, *minx)) {
                  f = minx;
                  do
                    miny = f;
                  while (++f != l && !t.greater_y_2_object()(*f, *miny));
                  break;
                }
              } // for (;;)
              break;
            } // if (t.less_y_2_object()(*f, *maxy))
          } // for (;;)
          break;
        } // if (t.less_x_2_object()(*f, *maxx))
      } // for (;;)
    else
      // fourth quadrant
      for (;;) {
        miny = f;
        if (++f == l) {
          maxx = maxy = minx = miny;
          break;
        }
        if (t.greater_y_2_object()(*f, *miny)) {
          f = miny;
          for (;;) {
            maxx = f;
            if (++f == l) {
              maxy = minx = maxx;
              break;
            }
            if (t.less_x_2_object()(*f, *maxx)) {
              f = maxx;
              for (;;) {
                maxy = f;
                if (++f == l) {
                  minx = maxy;
                  break;
                }
                if (t.less_y_2_object()(*f, *maxy)) {
                  f = maxy;
                  do
                    minx = f;
                  while (++f != l && !t.greater_x_2_object()(*f, *minx));
                  break;
                }
              } // for (;;)
              break;
            } // if (t.less_x_2_object()(*f, *maxx))
          } // for (;;)
          break;
        } // if (t.greater_y_2_object()(*f, *miny))
      } // for (;;)
  else
    if (t.less_y_2_object()(*f, *minx))
      // third quadrant
      for (;;) {
        minx = f;
        if (++f == l) {
          miny = maxx = maxy = minx;
          break;
        }
        if (t.greater_x_2_object()(*f, *minx)) {
          f = minx;
          for (;;) {
            miny = f;
            if (++f == l) {
              maxx = maxy = miny;
              break;
            }
            if (t.greater_y_2_object()(*f, *miny)) {
              f = miny;
              for (;;) {
                maxx = f;
                if (++f == l) {
                  maxy = maxx;
                  break;
                }
                if (t.less_x_2_object()(*f, *maxx)) {
                  f = maxx;
                  do
                    maxy = f;
                  while (++f != l && !t.less_y_2_object()(*f, *maxy));
                  break;
                }
              } // for (;;)
              break;
            } // if (t.greater_y_2_object()(*f, *miny))
          } // for (;;)
          break;
        } // if (t.greater_x_2_object()(*f, *minx))
      } // for (;;)
    else
      // second quadrant
      for (;;) {
        maxy = f;
        if (++f == l) {
          minx = miny = maxx = maxy;
          break;
        }
        if (t.less_y_2_object()(*f, *maxy)) {
          f = maxy;
          for (;;) {
            minx = f;
            if (++f == l) {
              miny = maxx = minx;
              break;
            }
            if (t.greater_x_2_object()(*f, *minx)) {
              f = minx;
              for (;;) {
                miny = f;
                if (++f == l) {
                  maxx = miny;
                  break;
                }
                if (t.greater_y_2_object()(*f, *miny)) {
                  f = miny;
                  do
                    maxx = f;
                  while (++f != l && !t.less_x_2_object()(*f, *maxx));
                  break;
                }
              } // for (;;)
              break;
            } // if (t.greater_x_2_object()(*f, *minx))
          } // for (;;)
          break;
        } // if (t.less_y_2_object()(*f, *maxy))
      } // for (;;)

  // Output
  *o++ = t.less_x_2_object()(*first, *minx) ? first : minx;
  *o++ = t.less_y_2_object()(*first, *miny) ? first : miny;
  *o++ = t.less_x_2_object()(*maxx, *first) ? first : maxx;
  *o++ = t.less_y_2_object()(*maxy, *first) ? first : maxy;
  return o;
} // convex_bounding_box_2(f, l, o, t)

template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
min_rectangle_2(
  ForwardIterator f,
  ForwardIterator l,
  OutputIterator o,
  Traits& t)
{
  // check for trivial cases
  if (f == l) return o;
  ForwardIterator tst = f;
  if (++tst == l) {
    *o++ = *f;
    return o;
  }

  // types from the traits class
  typedef typename Traits::Rectangle_2  Rectangle_2;
  typedef typename Traits::Direction_2  Direction_2;

  // quadruple of points defining the current rectangle
  ForwardIterator curr[4];
  // initialised to the points defining the bounding box
  convex_bounding_box_2(f, l, curr, t);

  // curr[i] can be advanced (cyclically) until it reaches limit[i]
  ForwardIterator limit[4];
  limit[0] = curr[1], limit[1] = curr[2],
    limit[2] = curr[3], limit[3] = curr[0];

  // quadruple of direction candidates defining the current rectangle
  Direction_2  dir[4];
  for (int i = 0; i < 4; i++) {
    ForwardIterator cp = curr[i];
    if (++cp == l)
      cp = f;
    dir[i] = t.construct_direction_2_object()(*(curr[i]), *cp);
    dir[i] = t.rotate_direction_by_multiple_of_pi_2_object()(dir[i], i);
  }

  int  yet_to_finish = 0;
  for (int i1 = 0; i1 < 4; ++i1) {
    CGAL_optimisation_assertion(limit[i1] != l);
    if (curr[i1] != limit[i1])
      ++yet_to_finish;
  }

  int low = t.less_rotate_ccw_2_object()(dir[0], dir[1]) ? 0 : 1;
  int upp = t.less_rotate_ccw_2_object()(dir[2], dir[3]) ? 2 : 3;

  int event =
    t.less_rotate_ccw_2_object()(dir[low], dir[upp]) ? low : upp;

  Rectangle_2 rect_so_far =
    t.construct_rectangle_2_object()(
      *(curr[0]), dir[event], *(curr[1]), *(curr[2]), *(curr[3]));

  for (;;) {
    if (++curr[event] == l)
      curr[event] = f;
    ForwardIterator cp = curr[event];
    if (++cp == l)
      cp = f;

    dir[event] = t.construct_direction_2_object()(*(curr[event]), *cp);
    dir[event] = t.rotate_direction_by_multiple_of_pi_2_object()(
      dir[event], event);

    if (curr[event] == limit[event])
      if (--yet_to_finish == 0)
        break;

    if (event < 2)
      low = t.less_rotate_ccw_2_object()(dir[0], dir[1]) ? 0 : 1;
    else
      upp = t.less_rotate_ccw_2_object()(dir[2], dir[3]) ? 2 : 3;

    event = t.less_rotate_ccw_2_object()(dir[low], dir[upp]) ? low : upp;

    Rectangle_2 test_rect = t.construct_rectangle_2_object()(
      *(curr[0]), dir[event], *(curr[1]), *(curr[2]), *(curr[3]));
    if (t.area_less_rectangle_2_object()(test_rect, rect_so_far))
      rect_so_far = test_rect;

  } // for (;;)

  return t.copy_rectangle_vertices_2(rect_so_far, o);

} // min_rectangle_2( f, l, o , t)
CGAL_END_NAMESPACE
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/leda_window.h>
#include <CGAL/leda_real.h>
CGAL_BEGIN_NAMESPACE

template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
min_parallelogram_2(ForwardIterator f,
                    ForwardIterator l,
                    OutputIterator o,
                    Traits& t)
{
  // check for trivial cases
  if (f == l) return o;
  ForwardIterator tst = f;
  if (++tst == l) {
    *o++ = *f;
    return o;
  }

  // types from the traits class
  typedef typename Traits::Parallelogram_2  Parallelogram_2;
  typedef typename Traits::Direction_2      Direction_2;

  // quadruple of points defining the bounding box
  ForwardIterator curr[4];
  // initialised to the points defining the bounding box
  convex_bounding_box_2(f, l, curr, t);

#ifdef CGAL_TRACE
  /*
  ForwardIterator mmix = std::min_element(f, l, t.less_x_2_object());
  ForwardIterator mmax = std::max_element(f, l, t.less_x_2_object());
  ForwardIterator mmiy = std::min_element(f, l, t.less_y_2_object());
  ForwardIterator mmay = std::max_element(f, l, t.less_y_2_object());
  CGAL_assertion(!t.less_x_2_object()(*mmix, *(curr[0])));
  CGAL_assertion(!t.less_x_2_object()(*(curr[0]), *mmix));
  CGAL_assertion(!t.less_x_2_object()(*mmax, *(curr[2])));
  CGAL_assertion(!t.less_x_2_object()(*(curr[2]), *mmax));
  CGAL_assertion(!t.less_y_2_object()(*mmiy, *(curr[1])));
  CGAL_assertion(!t.less_y_2_object()(*(curr[1]), *mmiy));
  CGAL_assertion(!t.less_y_2_object()(*mmay, *(curr[3])));
  CGAL_assertion(!t.less_y_2_object()(*(curr[3]), *mmay));
  */
#endif

  ForwardIterator low   = curr[1];
  ForwardIterator upp   = curr[3];
  ForwardIterator right = low;
  ForwardIterator left  = upp;

  int yet_to_finish = 2;

  // initialize parallelogram
  ForwardIterator ln = low;
  do
    if (++ln == l)
      ln = f;
  while (t.equal_2_object()(*ln, *low));
  Direction_2 d_low = t.construct_direction_2_object()(*low, *ln);
  ForwardIterator un = upp;
  do
    if (++un == l)
      un = f;
  while (t.equal_2_object()(*un, *upp));
  Direction_2 d_upp = t.construct_direction_2_object()(*un, *upp);

  bool low_goes_next = t.less_rotate_ccw_2_object()(d_low, d_upp);
  Direction_2 next_dir = low_goes_next ? d_low : d_upp;

  Direction_2 d_leftright = next_dir;
  for (;;) {
    // compute the next left/right candidate and store it to d_leftright
    ForwardIterator rig = right;
    do
      if (++rig == l)
        rig = f;
    while (t.equal_2_object()(*rig, *right));
    Direction_2 d_right = t.construct_direction_2_object()(*right, *rig);
  
    ForwardIterator len = left;
    do
      if (++len == l)
        len = f;
    while (t.equal_2_object()(*len, *left));
    Direction_2 d_left = t.construct_direction_2_object()(*len, *left);
  
    if (t.less_rotate_ccw_2_object()(d_right, d_left))
      if (t.right_of_implicit_line_2_object()(*rig, *left, next_dir))
        right = rig;
      else {
        d_leftright = d_right;
        break;
      }
    else
      if (t.right_of_implicit_line_2_object()(*right, *len, next_dir))
        left = len;
      else {
        d_leftright = d_left;
        break;
      }
  } // for (;;)

  Parallelogram_2 para_so_far =
    t.construct_parallelogram_2_object()(
      *low, next_dir, *right, d_leftright, *upp, *left);

  for (;;) {
    if (low_goes_next) {
      low = ln;
      if (low == curr[3])
        if (--yet_to_finish == 0)
          break;
    } else {
      upp = un;
      if (upp == curr[1])
        if (--yet_to_finish == 0)
          break;
    }

    // compute the next lower/upper candidate
    ln = low;
    do
      if (++ln == l)
        ln = f;
    while (t.equal_2_object()(*ln, *low));
    d_low = t.construct_direction_2_object()(*low, *ln);
    un = upp;
    do
      if (++un == l)
        un = f;
    while (t.equal_2_object()(*un, *upp));
    d_upp = t.construct_direction_2_object()(*un, *upp);

    low_goes_next = t.less_rotate_ccw_2_object()(d_low, d_upp);
    next_dir = low_goes_next ? d_low : d_upp;

    for (;;) {
      // compute the next left/right candidate and store it to d_leftright
      ForwardIterator rig = right;
      do
        if (++rig == l)
          rig = f;
      while (t.equal_2_object()(*rig, *right));
      Direction_2 d_right = t.construct_direction_2_object()(*right, *rig);
    
      ForwardIterator len = left;
      do
        if (++len == l)
          len = f;
      while (t.equal_2_object()(*len, *left));
      Direction_2 d_left = t.construct_direction_2_object()(*len, *left);
    
      if (t.less_rotate_ccw_2_object()(d_right, d_left))
        if (t.right_of_implicit_line_2_object()(*rig, *left, next_dir))
          right = rig;
        else {
          d_leftright = d_right;
          break;
        }
      else
        if (t.right_of_implicit_line_2_object()(*right, *len, next_dir))
          left = len;
        else {
          d_leftright = d_left;
          break;
        }
    } // for (;;)

    // check whether we found a smaller parallelogram
    Parallelogram_2 test_para =
      t.construct_parallelogram_2_object()(
        *low, next_dir, *right, d_leftright, *upp, *left);

#ifdef CGAL_TRACE
    {
      typedef typename
        std::iterator_traits< ForwardIterator >::value_type Point;
      typedef Polygon_traits_2< typename Traits::R > P_traits;
      typedef std::vector< Point >                   Cont;
      typedef CGAL::Polygon_2< P_traits, Cont >      Polygon_2;
      Polygon_2 p;
      t.copy_parallelogram_vertices_2(test_para, std::back_inserter(p));
      CGAL_assertion(p.is_simple());
      CGAL_assertion(p.is_convex());
      cout << "p_area = " << p.area() << endl;
      for (ForwardIterator ii = f; ii != l; ++ii)
        CGAL_assertion(!p.has_on_unbounded_side(*ii));
    }
#endif // CGAL_TRACE

    if (t.area_less_parallelogram_2_object()(test_para, para_so_far))
      para_so_far = test_para;

  } // for (;;)

#ifdef CGAL_TRACE
   typedef typename
     std::iterator_traits< ForwardIterator >::value_type Point;
   Point p[4];
   t.copy_parallelogram_vertices_2(para_so_far, p);
   leda_window w;
   w.init(-50, 450, -35);
   w.display();
   Ostream_iterator< Point, leda_window > oip(w);
   //std::ostream_iterator< Point > oipc(std::cerr, "\n");
   std::copy(f, l, oip);
   w << YELLOW;
   w.set_node_width(7);
   {
     ForwardIterator ii = curr[0];
     while (ii != curr[2]) {
       *oip++ = *ii;
       if (++ii == l) ii = f;
     }
     *oip++ = *ii;
   }
   w.set_node_width(5);
   w << GREEN << para_so_far.p1 << para_so_far.p2
     << para_so_far.p3 << para_so_far.p4;
   {
     typedef typename Traits::Line_2 Line_2;
     Line_2 l1(para_so_far.p1, para_so_far.d1);
     Line_2 l2(para_so_far.p2, para_so_far.d2);
     Line_2 l3(para_so_far.p3, para_so_far.d1);
     Line_2 l4(para_so_far.p4, para_so_far.d2);
     if (l1 == l2) cout << "l1 == l2" << endl;
     if (l1 == l3) cout << "l1 == l3" << endl;
     if (l1 == l4) cout << "l1 == l4" << endl;
     w << BLUE << l1 << l2 << l3 << l4;
   }
   w.set_node_width(3);
   w << RED << p[0] << p[1] << p[2] << p[3];
   w.read_mouse();
   std::cerr << "ZAP" << std::endl;
   *o++ = p[0];
   *o++ = p[1];
   *o++ = p[2];
   *o++ = p[3];
   return o;
 #else
   return t.copy_parallelogram_vertices_2(para_so_far, o);
 #endif

 } // min_parallelogram_2(f, l, o , t)
template < class ForwardIterator, class OutputIterator, class Traits >
OutputIterator
min_strip_2(ForwardIterator f,
            ForwardIterator l,
            OutputIterator o,
            Traits& t)
{
  // check for trivial cases
  if (f == l) return o;
  ForwardIterator tst = f;
  if (++tst == l) return o;

  // types from the traits class
  typedef typename Traits::Strip_2        Strip_2;
  typedef typename Traits::Direction_2    Direction_2;

  // quadruple of points defining the bounding box
  ForwardIterator curr[4];
  // initialised to the points defining the bounding box
  convex_bounding_box_2(f, l, curr, t);

  ForwardIterator low = curr[0];
  ForwardIterator upp = curr[2];

  int yet_to_finish = 2;

  ForwardIterator nlow = low;
  if (++nlow == l)
    nlow = f;
  Direction_2 low_dir = t.construct_direction_2_object()(*low, *nlow);
  ForwardIterator nupp = upp;
  if (++nupp == l)
    nupp = f;
  Direction_2 upp_dir = t.construct_direction_2_object()(*nupp, *upp);

  bool low_goes_next = t.less_rotate_ccw_2_object()(low_dir, upp_dir);
  Strip_2 strip_so_far =
    low_goes_next ?
      t.construct_strip_2_object()(*low, low_dir, *upp) :
      t.construct_strip_2_object()(*low, upp_dir, *upp);

  for (;;) {
    // compute next direction
    if (low_goes_next) {
      low = nlow;
      if (low == curr[2])
        if (--yet_to_finish == 0)
          break;
      if (++nlow == l)
        nlow = f;
      low_dir = t.construct_direction_2_object()(*low, *nlow);
    } else {
      upp = nupp;
      if (upp == curr[0])
        if (--yet_to_finish == 0)
          break;
      if (++nupp == l)
        nupp = f;
      upp_dir = t.construct_direction_2_object()(*nupp, *upp);
    }

    low_goes_next = t.less_rotate_ccw_2_object()(low_dir, upp_dir);
    Strip_2 test_strip =
    low_goes_next ?
      t.construct_strip_2_object()(*low, low_dir, *upp) :
      t.construct_strip_2_object()(*low, upp_dir, *upp);
    if (t.width_less_strip_2_object()(test_strip, strip_so_far))
      strip_so_far = test_strip;

  } // for (;;)

  // return the result
  return t.copy_strip_lines_2(strip_so_far, o);

} // min_strip_2(f, l, o, t)


#ifdef CGAL_REP_CLASS_DEFINED

CGAL_END_NAMESPACE
#include <CGAL/Min_quadrilateral_traits_2.h>
CGAL_BEGIN_NAMESPACE

template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
min_rectangle_2(ForwardIterator f,
         ForwardIterator l,
         OutputIterator o)
{
  typedef std::iterator_traits< ForwardIterator >::value_type VT;
  typedef typename VT::R R;
  Min_quadrilateral_default_traits_2< R > t;
  return min_rectangle_2(f, l, o, t);
} // min_rectangle_2(f, l, o)

// backwards compatibility
template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
minimum_enclosing_rectangle_2(ForwardIterator f,
                       ForwardIterator l,
                       OutputIterator o)
{ return min_rectangle_2(f, l, o); }
template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
min_parallelogram_2(ForwardIterator f,
         ForwardIterator l,
         OutputIterator o)
{
  typedef std::iterator_traits< ForwardIterator >::value_type VT;
  typedef typename VT::R R;
  Min_quadrilateral_default_traits_2< R > t;
  return min_parallelogram_2(f, l, o, t);
} // min_parallelogram_2(f, l, o)

// backwards compatibility
template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
minimum_enclosing_parallelogram_2(ForwardIterator f,
                       ForwardIterator l,
                       OutputIterator o)
{ return min_parallelogram_2(f, l, o); }
template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
min_strip_2(ForwardIterator f,
         ForwardIterator l,
         OutputIterator o)
{
  typedef std::iterator_traits< ForwardIterator >::value_type VT;
  typedef typename VT::R R;
  Min_quadrilateral_default_traits_2< R > t;
  return min_strip_2(f, l, o, t);
} // min_strip_2(f, l, o)

// backwards compatibility
template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
minimum_enclosing_strip_2(ForwardIterator f,
                       ForwardIterator l,
                       OutputIterator o)
{ return min_strip_2(f, l, o); }

#endif // CGAL_REP_CLASS_DEFINED

CGAL_END_NAMESPACE

#endif // ! (CGAL_MIN_QUADRILATERAL_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

