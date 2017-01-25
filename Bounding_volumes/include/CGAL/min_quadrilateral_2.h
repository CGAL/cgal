// Copyright (c) 1999-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch> and
//                 Emo Welzl <emo@inf.ethz.ch>

#ifndef CGAL_MIN_QUADRILATERAL_2_H
#define CGAL_MIN_QUADRILATERAL_2_H 1

#include <CGAL/license/Bounding_volumes.h>


#include <CGAL/basic.h>
#include <CGAL/Optimisation/assertions.h>
#include <iterator>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#ifdef CGAL_OPTIMISATION_EXPENSIVE_PRECONDITION_TAG
#include <CGAL/Polygon_2_algorithms.h>
#endif

namespace CGAL {

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
//    - smallest y coordinate
//    - largest x coordinate
//    - largest y coordinate
//    - smallest x coordinate
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
  ForwardIterator maxx;
  ForwardIterator miny;
  ForwardIterator maxy;

  typedef typename Traits::Point_2                Point_2;
  typedef typename Traits::Less_xy_2              Less_xy_2;
  typedef typename Traits::Less_yx_2              Less_yx_2;
  typedef boost::function2<bool,Point_2,Point_2>  Greater_xy_2;
  typedef boost::function2<bool,Point_2,Point_2>  Greater_yx_2;

  Less_xy_2    less_xy_2    = t.less_xy_2_object();
  Less_yx_2    less_yx_2    = t.less_yx_2_object();
  Greater_xy_2 greater_xy_2 = boost::bind(less_xy_2, _2, _1);
  Greater_yx_2 greater_yx_2 = boost::bind(less_yx_2, _2, _1);

  if (less_xy_2(*minx, *f) ||
      (less_yx_2(*minx, *f) && !less_xy_2(*f, *minx)))
    if (less_yx_2(*minx, *f))
      // first quadrant
      for (;;) {
        maxx = f;
        if (++f == l) {
          maxy = minx = miny = maxx;
          break;
        }
        if (less_xy_2(*f, *maxx)) {
          f = maxx;
          for (;;) {
            maxy = f;
            if (++f == l) {
              minx = miny = maxy;
              break;
            }
            if (less_yx_2(*f, *maxy)) {
              f = maxy;
              for (;;) {
                minx = f;
                if (++f == l) {
                  miny = minx;
                  break;
                }
                if (greater_xy_2(*f, *minx)) {
                  f = minx;
                  do
                    miny = f;
                  while (++f != l && !greater_yx_2(*f, *miny));
                  break;
                }
              } // for (;;)
              break;
            } // if (less_yx_2(*f, *maxy))
          } // for (;;)
          break;
        } // if (less_xy_2(*f, *maxx))
      } // for (;;)
    else
      // fourth quadrant
      for (;;) {
        miny = f;
        if (++f == l) {
          maxx = maxy = minx = miny;
          break;
        }
        if (greater_yx_2(*f, *miny)) {
          f = miny;
          for (;;) {
            maxx = f;
            if (++f == l) {
              maxy = minx = maxx;
              break;
            }
            if (less_xy_2(*f, *maxx)) {
              f = maxx;
              for (;;) {
                maxy = f;
                if (++f == l) {
                  minx = maxy;
                  break;
                }
                if (less_yx_2(*f, *maxy)) {
                  f = maxy;
                  do
                    minx = f;
                  while (++f != l && !greater_xy_2(*f, *minx));
                  break;
                }
              } // for (;;)
              break;
            } // if (less_xy_2(*f, *maxx))
          } // for (;;)
          break;
        } // if (greater_yx_2(*f, *miny))
      } // for (;;)
  else
    if (less_yx_2(*f, *minx))
      // third quadrant
      for (;;) {
        minx = f;
        if (++f == l) {
          miny = maxx = maxy = minx;
          break;
        }
        if (greater_xy_2(*f, *minx)) {
          f = minx;
          for (;;) {
            miny = f;
            if (++f == l) {
              maxx = maxy = miny;
              break;
            }
            if (greater_yx_2(*f, *miny)) {
              f = miny;
              for (;;) {
                maxx = f;
                if (++f == l) {
                  maxy = maxx;
                  break;
                }
                if (less_xy_2(*f, *maxx)) {
                  f = maxx;
                  do
                    maxy = f;
                  while (++f != l && !less_yx_2(*f, *maxy));
                  break;
                }
              } // for (;;)
              break;
            } // if (greater_yx_2(*f, *miny))
          } // for (;;)
          break;
        } // if (greater_xy_2(*f, *minx))
      } // for (;;)
    else
      // second quadrant
      for (;;) {
        maxy = f;
        if (++f == l) {
          minx = miny = maxx = maxy;
          break;
        }
        if (less_yx_2(*f, *maxy)) {
          f = maxy;
          for (;;) {
            minx = f;
            if (++f == l) {
              miny = maxx = minx;
              break;
            }
            if (greater_xy_2(*f, *minx)) {
              f = minx;
              for (;;) {
                miny = f;
                if (++f == l) {
                  maxx = miny;
                  break;
                }
                if (greater_yx_2(*f, *miny)) {
                  f = miny;
                  do
                    maxx = f;
                  while (++f != l && !less_xy_2(*f, *maxx));
                  break;
                }
              } // for (;;)
              break;
            } // if (greater_xy_2(*f, *minx))
          } // for (;;)
          break;
        } // if (less_yx_2(*f, *maxy))
      } // for (;;)

  // Output
  *o++ = less_yx_2(*first, *miny) ? first : miny;
  *o++ = less_xy_2(*maxx, *first) ? first : maxx;
  *o++ = less_yx_2(*maxy, *first) ? first : maxy;
  *o++ = less_xy_2(*first, *minx) ? first : minx;
  return o;
} // convex_bounding_box_2(f, l, o, t)

namespace Optimisation {
  // Adds certain redundant functionality for convenience
  template < typename Traits >
  struct Min_quadrilateral_traits_wrapper : public Traits
  {
    typedef Traits                                      Base;
    // types inherited from Traits
    typedef typename Base::Point_2                      Point_2;
    typedef typename Base::Direction_2                  Direction_2;
    // predicates and constructions inherited from Traits
    typedef typename Base::Has_on_negative_side_2       HONS;
    typedef typename Base::Construct_vector_2           CV2;
    typedef typename Base::Construct_direction_2        CD2;
    typedef typename Base::Construct_line_2             Construct_line_2;
    typedef typename Base::Compare_angle_with_x_axis_2  CAWXA;

    using Traits::has_on_negative_side_2_object;
    using Traits::construct_line_2_object;
    using Traits::construct_vector_2_object;
    using Traits::compare_angle_with_x_axis_2_object;

    Min_quadrilateral_traits_wrapper(const Traits& bt) : Base(bt) {}

    // ---------------------------------------------------------------
    // Right_of_implicit_line_2
    // ---------------------------------------------------------------
    typedef boost::function3<bool,Point_2,Point_2,Direction_2> 
      Right_of_implicit_line_2;
    
    Right_of_implicit_line_2 right_of_implicit_line_2_object() const {
      return boost::bind(has_on_negative_side_2_object(),
			 boost::bind(construct_line_2_object(), _2, _3),
			 _1);
    }
    
    typedef boost::function2<Direction_2,Point_2,Point_2> 
      Construct_direction_2;
    
    Construct_direction_2 construct_direction_2_object() const {
      return boost::bind(Base::construct_direction_2_object(),
			 boost::bind(construct_vector_2_object(), _1, _2));
    }
    
    template < class Kernel >
    class Rdbmop
    : public std::binary_function< Direction_2, int, Direction_2 >
    {
      typename Kernel::Construct_perpendicular_vector_2   cperpvec;
      typename Kernel::Construct_vector_from_direction_2  cvec;
      typename Kernel::Construct_direction_2              dir;
      typename Kernel::Construct_opposite_direction_2     oppdir;
    public:
    
      Rdbmop() {}
    
      Rdbmop(const Kernel& k)
      : cperpvec(k.construct_perpendicular_vector_2_object()),
        cvec(k.construct_vector_from_direction_2_object()),
        dir(k.construct_direction_2_object()),
        oppdir(k.construct_opposite_direction_2_object())
      {}
    
      Direction_2
      operator()(const Direction_2& d, int i) const
      {
        // FIXME: here I would like to construct a vector from a
        // direction, but this is not in the kernel concept
        // maybe, we can get rid of directions soon...
        CGAL_precondition(i >= 0 && i < 4);
        if (i == 0) return d;
        if (i == 1) return dir(cperpvec(cvec(d), CLOCKWISE));
        if (i == 2) return oppdir(d);
        return dir(cperpvec(cvec(d), COUNTERCLOCKWISE));
      }
    };
    
    typedef Rdbmop<Traits> Rotate_direction_by_multiple_of_pi_2;
    
    Rotate_direction_by_multiple_of_pi_2
    rotate_direction_by_multiple_of_pi_2_object() const
    { return Rotate_direction_by_multiple_of_pi_2(*this); }
    
    typedef boost::function2<bool,Direction_2,Direction_2>
      Less_angle_with_x_axis_2;
    Less_angle_with_x_axis_2 less_angle_with_x_axis_2_object() const {
      return boost::bind(std::equal_to<Comparison_result>(), 
                         boost::bind(compare_angle_with_x_axis_2_object(),
                                     _1, _2),
                         SMALLER);
    }

  };
} // namespace Optimisation

template < class ForwardIterator, class OutputIterator, class BTraits >
OutputIterator
min_rectangle_2(
  ForwardIterator f,
  ForwardIterator l,
  OutputIterator o,
  BTraits& bt)
{
  typedef Optimisation::Min_quadrilateral_traits_wrapper<BTraits> Traits;
  Traits t(bt);
  CGAL_optimisation_expensive_precondition(is_convex_2(f, l, t));
  CGAL_optimisation_expensive_precondition(
    orientation_2(f, l, t) == COUNTERCLOCKWISE);

  // check for trivial cases
  if (f == l) return o;
  ForwardIterator tst = f;
  if (++tst == l) {
    // all points are equal
    for (int i = 0; i < 4; ++i) *o++ = *f;
    return o;
  }

  // types from the traits class
  typedef typename Traits::Rectangle_2            Rectangle_2;
  typedef typename Traits::Direction_2            Direction_2;
  typedef typename Traits::Construct_direction_2  Construct_direction_2;
  typedef typename Traits::Construct_rectangle_2  Construct_rectangle_2;

  Construct_direction_2 direction = t.construct_direction_2_object();
  Construct_rectangle_2 rectangle = t.construct_rectangle_2_object();
  typename Traits::Rotate_direction_by_multiple_of_pi_2
    rotate = t.rotate_direction_by_multiple_of_pi_2_object();
  typename Traits::Less_angle_with_x_axis_2
    less_angle = t.less_angle_with_x_axis_2_object();
  typename Traits::Area_less_rectangle_2
    area_less = t.area_less_rectangle_2_object();

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
    dir[i] = rotate(direction(*(curr[i]), *cp), i);
  }

  int yet_to_finish = 0;
  for (int i1 = 0; i1 < 4; ++i1) {
    CGAL_optimisation_assertion(limit[i1] != l);
    if (curr[i1] != limit[i1])
      ++yet_to_finish;
  }

  int low = less_angle(dir[0], dir[1]) ? 0 : 1;
  int upp = less_angle(dir[2], dir[3]) ? 2 : 3;

  int event = less_angle(dir[low], dir[upp]) ? low : upp;

  Rectangle_2 rect_so_far =
    rectangle(*(curr[0]), dir[event], *(curr[1]), *(curr[2]), *(curr[3]));

  for (;;) {
    if (++curr[event] == l)
      curr[event] = f;
    ForwardIterator cp = curr[event];
    if (++cp == l)
      cp = f;

    dir[event] = rotate(direction(*(curr[event]), *cp), event);

    if (curr[event] == limit[event])
      if (--yet_to_finish <= 0)
        break;

    if (event < 2)
      low = less_angle(dir[0], dir[1]) ? 0 : 1;
    else
      upp = less_angle(dir[2], dir[3]) ? 2 : 3;

    event = less_angle(dir[low], dir[upp]) ? low : upp;

    Rectangle_2 test_rect = rectangle(*(curr[0]), dir[event],
                                      *(curr[1]), *(curr[2]), *(curr[3]));
    if (area_less(test_rect, rect_so_far))
      rect_so_far = test_rect;

  } // for (;;)

  return t.copy_rectangle_vertices_2(rect_so_far, o);

} // min_rectangle_2( f, l, o , t)

template < class ForwardIterator, class OutputIterator, class BTraits >
OutputIterator
min_parallelogram_2(ForwardIterator f,
                    ForwardIterator l,
                    OutputIterator o,
                    BTraits& bt)
{
  typedef Optimisation::Min_quadrilateral_traits_wrapper<BTraits> Traits;
  Traits t(bt);
  CGAL_optimisation_expensive_precondition(is_convex_2(f, l, t));

  // types from the traits class
  typedef typename Traits::Direction_2            Direction_2;
  typedef typename Traits::Parallelogram_2        Parallelogram_2;
  typedef typename Traits::Construct_direction_2  Construct_direction_2;
  typedef typename Traits::Equal_2                Equal_2;

  Equal_2 equal = t.equal_2_object();
  Construct_direction_2 direction = t.construct_direction_2_object();
  typename Traits::Construct_parallelogram_2
    parallelogram = t.construct_parallelogram_2_object();
  typename Traits::Less_angle_with_x_axis_2
    less_angle = t.less_angle_with_x_axis_2_object();
  typename Traits::Area_less_parallelogram_2
    area_less = t.area_less_parallelogram_2_object();
  typename Traits::Right_of_implicit_line_2
    right_of_line = t.right_of_implicit_line_2_object();

  // check for trivial cases
  if (f == l) return o;
  
  ForwardIterator first;
  do {
    first = f;
    if (++f == l) {
      // all points are equal
      for (int i = 0; i < 4; ++i) *o++ = *first;
      return o;
    }
  } while (equal(*first, *f));

  // quadruple of points defining the bounding box
  ForwardIterator curr[4];
  // initialised to the points defining the bounding box
  convex_bounding_box_2(first, l, curr, t);


  ForwardIterator low   = curr[0];
  ForwardIterator upp   = curr[2];
  ForwardIterator right = low;
  ForwardIterator left  = upp;

  int yet_to_finish = 2;

  // initialize parallelogram
  ForwardIterator ln = low;
  do
    if (++ln == l)
      ln = first;
  while (equal(*ln, *low));
  Direction_2 d_low = direction(*low, *ln);
  ForwardIterator un = upp;
  do
    if (++un == l)
      un = first;
  while (equal(*un, *upp));
  Direction_2 d_upp = direction(*un, *upp);

  bool low_goes_next = less_angle(d_low, d_upp);
  Direction_2 next_dir = low_goes_next ? d_low : d_upp;

  Direction_2 d_leftright = next_dir;
  for (;;) {
    // compute the next left/right candidate and store it to d_leftright
    ForwardIterator rig = right;
    do
      if (++rig == l)
        rig = first;
    while (equal(*rig, *right));
    Direction_2 d_right = direction(*right, *rig);
  
    ForwardIterator len = left;
    do
      if (++len == l)
        len = first;
    while (equal(*len, *left));
    Direction_2 d_left = direction(*len, *left);
  
    if (less_angle(d_right, d_left))
      if (right_of_line(*rig, *left, next_dir))
        right = rig;
      else {
        d_leftright = d_right;
        break;
      }
    else
      if (right_of_line(*right, *len, next_dir))
        left = len;
      else {
        d_leftright = d_left;
        break;
      }
  } // for (;;)

  Parallelogram_2 para_so_far =
    parallelogram(*low, next_dir, *right, d_leftright, *upp, *left);

  for (;;) {
    if (low_goes_next) {
      low = ln;
      if (low == curr[2])
        if (--yet_to_finish <= 0)
          break;
    } else {
      upp = un;
      if (upp == curr[0])
        if (--yet_to_finish <= 0)
          break;
    }

    // compute the next lower/upper candidate
    ln = low;
    do
      if (++ln == l)
        ln = first;
    while (equal(*ln, *low));
    d_low = direction(*low, *ln);
    un = upp;
    do
      if (++un == l)
        un = first;
    while (equal(*un, *upp));
    d_upp = direction(*un, *upp);

    low_goes_next = less_angle(d_low, d_upp);
    next_dir = low_goes_next ? d_low : d_upp;

    for (;;) {
      // compute the next left/right candidate and store it to d_leftright
      ForwardIterator rig = right;
      do
        if (++rig == l)
          rig = first;
      while (equal(*rig, *right));
      Direction_2 d_right = direction(*right, *rig);
    
      ForwardIterator len = left;
      do
        if (++len == l)
          len = first;
      while (equal(*len, *left));
      Direction_2 d_left = direction(*len, *left);
    
      if (less_angle(d_right, d_left))
        if (right_of_line(*rig, *left, next_dir))
          right = rig;
        else {
          d_leftright = d_right;
          break;
        }
      else
        if (right_of_line(*right, *len, next_dir))
          left = len;
        else {
          d_leftright = d_left;
          break;
        }
    } // for (;;)

    // check whether we found a smaller parallelogram
    Parallelogram_2 test_para =
      parallelogram(*low, next_dir, *right, d_leftright, *upp, *left);


    if (area_less(test_para, para_so_far))
      para_so_far = test_para;

  } // for (;;)


   return t.copy_parallelogram_vertices_2(para_so_far, o);
 } // min_parallelogram_2(f, l, o , t)
template < class ForwardIterator, class OutputIterator, class BTraits >
OutputIterator
min_strip_2(ForwardIterator f,
            ForwardIterator l,
            OutputIterator o,
            BTraits& bt)
{
  typedef Optimisation::Min_quadrilateral_traits_wrapper<BTraits> Traits;
  Traits t(bt);
  CGAL_optimisation_expensive_precondition(is_convex_2(f, l, t));

  // types from the traits class
  typedef typename Traits::Direction_2            Direction_2;
  typedef typename Traits::Strip_2                Strip_2;
  typedef typename Traits::Equal_2                Equal_2;
  typedef typename Traits::Construct_direction_2  Construct_direction_2;
  typedef typename Traits::Construct_strip_2      Construct_strip_2;
  typedef typename Traits::Width_less_strip_2     Width_less_strip_2;

  Equal_2 equal = t.equal_2_object();
  Construct_direction_2 direction = t.construct_direction_2_object();
  Construct_strip_2 strip = t.construct_strip_2_object();
  Width_less_strip_2 width_less = t.width_less_strip_2_object();
  typename Traits::Less_angle_with_x_axis_2
    less_angle = t.less_angle_with_x_axis_2_object();

  // check for trivial cases
  if (f == l) return o;
  ForwardIterator first;
  do {
    first = f;
    if (++f == l)
      // strip undefined, if no two distinct points exist
      return o;
  } while (equal(*first, *f));

  // quadruple of points defining the bounding box
  ForwardIterator curr[4];
  // initialised to the points defining the bounding box
  convex_bounding_box_2(first, l, curr, t);

  ForwardIterator low = curr[0];
  ForwardIterator upp = curr[2];

  int yet_to_finish = 2;

  ForwardIterator nlow = low;
  if (++nlow == l)
    nlow = first;
  Direction_2 low_dir = direction(*low, *nlow);
  ForwardIterator nupp = upp;
  if (++nupp == l)
    nupp = first;
  Direction_2 upp_dir = direction(*nupp, *upp);

  bool low_goes_next = less_angle(low_dir, upp_dir);
  Strip_2 strip_so_far = low_goes_next ?
    strip(*low, low_dir, *upp) : strip(*low, upp_dir, *upp);

  for (;;) {
    // compute next direction
    if (low_goes_next) {
      low = nlow;
      if (low == curr[2])
        if (--yet_to_finish <= 0)
          break;
      if (++nlow == l)
        nlow = first;
      low_dir = direction(*low, *nlow);
    } else {
      upp = nupp;
      if (upp == curr[0])
        if (--yet_to_finish <= 0)
          break;
      if (++nupp == l)
        nupp = first;
      upp_dir = direction(*nupp, *upp);
    }

    low_goes_next = less_angle(low_dir, upp_dir);
    Strip_2 test_strip = low_goes_next ?
      strip(*low, low_dir, *upp) : strip(*low, upp_dir, *upp);
    if (width_less(test_strip, strip_so_far))
      strip_so_far = test_strip;

  } // for (;;)

  // return the result
  return t.copy_strip_lines_2(strip_so_far, o);

} // min_strip_2(f, l, o, t)


} //namespace CGAL
#include <CGAL/Min_quadrilateral_traits_2.h>
namespace CGAL {

template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
min_rectangle_2(ForwardIterator f,
         ForwardIterator l,
         OutputIterator o)
{
  typedef typename std::iterator_traits< ForwardIterator >::value_type VT;
  typedef typename Kernel_traits<VT>::Kernel Kernel;
  Min_quadrilateral_default_traits_2<Kernel> t;
  return min_rectangle_2(f, l, o, t);
} // min_rectangle_2(f, l, o)

#ifndef CGAL_NO_DEPRECATED_CODE
// backwards compatibility
template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
minimum_enclosing_rectangle_2(ForwardIterator f,
                       ForwardIterator l,
                       OutputIterator o)
{ return min_rectangle_2(f, l, o); }
#endif // CGAL_NO_DEPRECATED_CODE
template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
min_parallelogram_2(ForwardIterator f,
         ForwardIterator l,
         OutputIterator o)
{
  typedef typename std::iterator_traits< ForwardIterator >::value_type VT;
  typedef typename Kernel_traits<VT>::Kernel Kernel;
  Min_quadrilateral_default_traits_2<Kernel> t;
  return min_parallelogram_2(f, l, o, t);
} // min_parallelogram_2(f, l, o)

#ifndef CGAL_NO_DEPRECATED_CODE
// backwards compatibility
template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
minimum_enclosing_parallelogram_2(ForwardIterator f,
                       ForwardIterator l,
                       OutputIterator o)
{ return min_parallelogram_2(f, l, o); }
#endif // CGAL_NO_DEPRECATED_CODE
template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
min_strip_2(ForwardIterator f,
         ForwardIterator l,
         OutputIterator o)
{
  typedef typename std::iterator_traits< ForwardIterator >::value_type VT;
  typedef typename Kernel_traits<VT>::Kernel Kernel;
  Min_quadrilateral_default_traits_2<Kernel> t;
  return min_strip_2(f, l, o, t);
} // min_strip_2(f, l, o)

#ifndef CGAL_NO_DEPRECATED_CODE
// backwards compatibility
template < class ForwardIterator, class OutputIterator >
inline
OutputIterator
minimum_enclosing_strip_2(ForwardIterator f,
                       ForwardIterator l,
                       OutputIterator o)
{ return min_strip_2(f, l, o); }
#endif // CGAL_NO_DEPRECATED_CODE

} //namespace CGAL

#endif // ! (CGAL_MIN_QUADRILATERAL_2_H)
