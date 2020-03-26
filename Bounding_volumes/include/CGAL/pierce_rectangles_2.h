// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_PIERCE_RECTANGLES_2_H
#define CGAL_PIERCE_RECTANGLES_2_H 1

#include <CGAL/license/Bounding_volumes.h>


#include <CGAL/Optimisation/assertions.h>
#include <CGAL/circulator.h>
#include <CGAL/algorithm.h>
#include <algorithm>
#include <iterator>
#include <vector>
#include <boost/bind.hpp>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4355) // complaint about using 'this' to
#endif                          // initialize a member

namespace CGAL {

//!!! STL-extensions
template < class T >
struct Wastebasket {
  typedef std::output_iterator_tag iterator_category;
  typedef Wastebasket< T >         iterator;

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

template < class Traits_ >
struct Loc_domain {
  // ---------------------------------------------
  // types:

  typedef Traits_                               Traits;
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
        if (traits.less_x_2_object()(*i, minx)) minx = *i;
        if (traits.less_y_2_object()(*i, miny)) miny = *i;
      }
      else {
        if (traits.less_y_2_object()(*i, miny)) miny = *i;
        if (traits.less_x_2_object()(maxx, *i)) maxx = *i;
      }
    else
      if (j == 2) {
        if (traits.less_x_2_object()(maxx, *i)) maxx = *i;
        if (traits.less_y_2_object()(maxy, *i)) maxy = *i;
      }
      else {
        if (traits.less_y_2_object()(maxy, *i)) maxy = *i;
        if (traits.less_x_2_object()(*i, minx)) minx = *i;
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
      if (traits.less_x_2_object()(*i, minx)) minx = *i;
      if (traits.less_x_2_object()(maxx, *i)) maxx = *i;
      if (traits.less_y_2_object()(*i, miny)) miny = *i;
      if (traits.less_y_2_object()(maxy, *i)) maxy = *i;
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
      return traits.construct_point_2_above_right_implicit_point_2_object()(
        minx, miny, r);
    else if (i == 1)
      return traits.construct_point_2_above_left_implicit_point_2_object()(
        maxx, miny, r);
    else if (i == 2)
      return traits.construct_point_2_below_left_implicit_point_2_object()(
        maxx, maxy, r);
    return traits.construct_point_2_below_right_implicit_point_2_object()(
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
        CGAL_optimisation_assertion(!traits.less_x_2_object()(*i, minx));
        CGAL_optimisation_assertion(!traits.less_x_2_object()(maxx, *i));
        CGAL_optimisation_assertion(!traits.less_y_2_object()(*i, miny));
        CGAL_optimisation_assertion(!traits.less_y_2_object()(maxy, *i));
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
template < class Traits_ >
struct Staircases : public Loc_domain< Traits_ > {
  typedef Traits_                           Traits;
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
  : Base(b, e, t),
    sorted(this->pts),
    xgy(t.signed_x_distance_2_object()(this->maxx, this->minx) >
        t.signed_y_distance_2_object()(this->maxy, this->miny))
  {
    using std::sort;
    using std::find_if;

    Container& xsort = xgy ? sorted : this->pts;
    Container& ysort = xgy ? this->pts : sorted;

    // build top-left and bottom-right staircases
    sort(ysort.begin(), ysort.end(), this->traits.less_y_2_object());
    // bottom-right
    Iterator i = ysort.begin();
    do {
      brstc.push_back(*i++);
      i = find_if(i, ysort.end(),
                  boost::bind(this->traits.less_x_2_object(), brstc.back(), _1));
    } while (i != ysort.end());
    // top-left
    Riterator j = ysort.rbegin();
    do {
      tlstc.push_back(*j++);
      j = find_if(j, ysort.rend(),
                  boost::bind(this->traits.less_x_2_object(), _1, tlstc.back()));
    } while (j != ysort.rend());

    // build left-bottom and right-top staircases
    sort(xsort.begin(), xsort.end(), this->traits.less_x_2_object());
    // left-bottom
    i = xsort.begin();
    do {
      lbstc.push_back(*i++);
      i = find_if(i, xsort.end(),
                  boost::bind(this->traits.less_y_2_object(), _1, lbstc.back()));
    } while (i != xsort.end());
    // right-top
    j = xsort.rbegin();
    do {
      rtstc.push_back(*j++);
      j = find_if(j, xsort.rend(),
                  boost::bind(this->traits.less_y_2_object(), rtstc.back(), _1));
    } while (j != xsort.rend());
  } // Staircases(b, e, t)

  bool is_middle_empty() const {
    //!!! the "middle" point could be precomputed in advance
    Citerator i = this->pts.begin();
    FT rr = FT(2) * this->r;
    do
      if (this->traits.signed_x_distance_2_object()(this->maxx, *i) > rr &&
          this->traits.signed_x_distance_2_object()(*i, this->minx) > rr &&
          this->traits.signed_y_distance_2_object()(*i, this->miny) > rr &&
          this->traits.signed_y_distance_2_object()(this->maxy, *i) > rr)
        return false;
    while (++i != this->pts.end());
    return true;
  } // is_middle()

  bool is_x_greater_y() const { return xgy; }

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
      this->traits.construct_point_2_above_right_implicit_point_2_object()(
        this->minx, this->miny, FT(2) * this->r);
    Point_2 q =
      this->traits.construct_point_2_above_left_implicit_point_2_object()(
        this->maxx, this->miny, FT(2) * this->r);

    Citerator i =
      min_element_if(
        this->pts.begin(), this->pts.end(),
        this->traits.less_x_2_object(),
        boost::bind(std::logical_and< bool >(),
                    boost::bind(this->traits.less_x_2_object(), p, _1),
                    boost::bind(this->traits.less_y_2_object(), p, _1)));
    Citerator j =
      max_element_if(
        this->pts.begin(), this->pts.end(),
        this->traits.less_x_2_object(),
        boost::bind(std::logical_and< bool >(),
                    boost::bind(this->traits.less_x_2_object(), _1, q),
                    boost::bind(this->traits.less_y_2_object(), q, _1)));
    return Intervall(i == this->pts.end() ? this->maxx : *i,
                     j == this->pts.end() ? this->minx : *j);
  } // top_intervall()

  Intervall bottom_intervall() const {
    Point_2 p =
      this->traits.construct_point_2_below_right_implicit_point_2_object()(
        this->minx, this->maxy, FT(2) * this->r);
    Point_2 q =
      this->traits.construct_point_2_below_left_implicit_point_2_object()(
        this->maxx, this->maxy, FT(2) * this->r);

    Citerator i =
      min_element_if(
        this->pts.begin(), this->pts.end(),
        this->traits.less_x_2_object(),
        boost::bind(std::logical_and< bool >(),
                    boost::bind(this->traits.less_x_2_object(), p, _1),
                    boost::bind(this->traits.less_y_2_object(), _1, p)));
    Citerator j =
      max_element_if(
        this->pts.begin(), this->pts.end(),
        this->traits.less_x_2_object(),
        boost::bind(std::logical_and< bool >(),
                    boost::bind(this->traits.less_x_2_object(), _1, q),
                    boost::bind(this->traits.less_y_2_object(), _1, q)));
    return Intervall(i == this->pts.end() ? this->maxx : *i,
                     j == this->pts.end() ? this->minx : *j);
  } // bottom_intervall()

  Intervall left_intervall() const {
    Point_2 p =
      this->traits.construct_point_2_above_left_implicit_point_2_object()(
        this->maxx, this->miny, FT(2) * this->r);
    Point_2 q =
      this->traits.construct_point_2_below_left_implicit_point_2_object()(
        this->maxx, this->maxy, FT(2) * this->r);

    Citerator i =
      min_element_if(
        this->pts.begin(), this->pts.end(),
        this->traits.less_y_2_object(),
        boost::bind(std::logical_and< bool >(),
                    boost::bind(this->traits.less_x_2_object(), _1, p),
                    boost::bind(this->traits.less_y_2_object(), p, _1)));
    Citerator j =
      max_element_if(
        this->pts.begin(), this->pts.end(),
        this->traits.less_y_2_object(),
        boost::bind(std::logical_and< bool >(),
                    boost::bind(this->traits.less_x_2_object(), _1, q),
                    boost::bind(this->traits.less_y_2_object(), _1, q)));
    return Intervall(i == this->pts.end() ? this->maxy : *i,
                     j == this->pts.end() ? this->miny : *j);
  } // left_intervall()

  Intervall right_intervall() const {
    Point_2 p =
      this->traits.construct_point_2_above_right_implicit_point_2_object()(
        this->minx, this->miny, FT(2) * this->r);
    Point_2 q =
      this->traits.construct_point_2_below_right_implicit_point_2_object()(
        this->minx, this->maxy, FT(2) * this->r);

    Citerator i =
      min_element_if(
        this->pts.begin(), this->pts.end(),
        this->traits.less_y_2_object(),
        boost::bind(std::logical_and< bool >(),
                    boost::bind(this->traits.less_x_2_object(), p, _1),
                    boost::bind(this->traits.less_y_2_object(), p, _1)));
    Citerator j =
      max_element_if(
        this->pts.begin(), this->pts.end(),
        this->traits.less_y_2_object(),
        boost::bind(std::logical_and< bool >(),
                    boost::bind(this->traits.less_x_2_object(), q, _1),
                    boost::bind(this->traits.less_y_2_object(), _1, q)));
    return Intervall(i == this->pts.end() ? this->maxy : *i,
                     j == this->pts.end() ? this->miny : *j);
  } // right_intervall()

  template < class OutputIterator >
  OutputIterator shared_intervall(OutputIterator o) const {
    if (xgy) {
      if (this->traits.signed_y_distance_2_object()(
          this->maxy, this->miny) > FT(4) * this->r)
        return o;
      Point_2 p =
        this->traits.construct_point_2_below_right_implicit_point_2_object()(
          this->minx, this->maxy, FT(2) * this->r);
      Point_2 q =
        this->traits.construct_point_2_above_left_implicit_point_2_object()(
          this->maxx, this->miny, FT(2) * this->r);
      //!!! start with binary search
      for (Citerator i = sorted.begin(); i != sorted.end(); ++i)
        if (this->traits.less_x_2_object()(p, *i) &&
            this->traits.less_x_2_object()(*i, q) &&
            !this->traits.less_y_2_object()(*i, p) &&
            !this->traits.less_y_2_object()(q, *i))
          *o++ = *i;
    } else {
      if (this->traits.signed_x_distance_2_object()(
          this->maxx, this->minx) > FT(4) * this->r)
        return o;
      Point_2 p =
        this->traits.construct_point_2_above_left_implicit_point_2_object()(
          this->maxx, this->miny, FT(2) * this->r);
      Point_2 q =
        this->traits.construct_point_2_below_right_implicit_point_2_object()(
          this->minx, this->maxy, FT(2) * this->r);
      //!!! start with binary search
      for (Citerator i = sorted.begin(); i != sorted.end(); ++i)
        if (!this->traits.less_x_2_object()(*i, p) &&
            !this->traits.less_x_2_object()(q, *i) &&
            this->traits.less_y_2_object()(p, *i) &&
            this->traits.less_y_2_object()(*i, q))
          *o++ = *i;
    }
    return o;
  } // shared_intervall(o)

private:
  Container tlstc, lbstc, brstc, rtstc, sorted;
  // exceeds the x-dimension of the location domain its y-dimension?
  bool xgy;
};


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
  using std::find_if;
  using std::less;

  typedef typename Traits::FT           FT;
  typename Traits::Infinity_distance_2 dist =
    d.traits.infinity_distance_2_object();
  typename Traits::Signed_infinity_distance_2 sdist =
    d.traits.signed_infinity_distance_2_object();

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
                boost::bind(less<FT>(),
                     d.r,
                     boost::bind(Min<FT>(),
                          boost::bind(dist, d[0], _1),
                          boost::bind(dist, d[2], _1)))))
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
                boost::bind(less<FT>(),
                     d.r,
                     boost::bind(Min<FT>(),
                          boost::bind(dist, d[1], _1),
                          boost::bind(dist, d[3], _1)))))
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
  using std::find_if;
  using std::less;
  using std::iter_swap;

  CGAL_optimisation_precondition(!d.empty());

  // typedefs:
  typedef typename Traits::FT                      FT;
  typedef typename Traits::Point_2                 Point_2;
  typedef typename Loc_domain< Traits >::Iterator  Iterator;
  typename Traits::Infinity_distance_2 dist =
    d.traits.infinity_distance_2_object();

  // test the four corners:
  for (int k = 0; k < 4; ++k) {

    // extract all points which are close enough to d[k]
    Point_2 corner = d[k];

    // find first point not covered by the rectangle at d[k]
    Iterator i = find_if(d.begin(), d.end(),
                         boost::bind(less<FT>(), d.r, boost::bind(dist, corner, _1)));

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

    CGAL_optimisation_expensive_assertion(
      save_end == find_if(d.end(), save_end,
                          boost::bind(less<FT>(), d.r, boost::bind(dist, corner, _1))));
    CGAL_optimisation_expensive_assertion(
      d.end() == find_if(d.begin(), d.end(),
                         boost::bind(std::greater_equal<FT>(),
                              d.r,
                              boost::bind(dist, corner, _1))));


    two_cover_points(d, o, ok);

    // restore saved sides of d:
    d.extreme(k) = save_side1;
    d.extreme((k+1) % 4) = save_side2;

    if (ok) {
      // does any rectangle contain the corner?
      if (d.end() != save_end) {
        *o++ = corner;
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
} //namespace CGAL
namespace CGAL {
template < class OutputIterator, class Traits >
OutputIterator
four_cover_points(Staircases< Traits >& d, OutputIterator o, bool& ok)
{

  using std::less;
  using std::iter_swap;
  using std::find_if;
  using std::back_inserter;

  typedef typename Traits::Point_2                  Point_2;
  typedef typename Traits::FT                       FT;
  typedef typename Traits::Less_x_2                 Less_x_2;
  typedef typename Traits::Less_y_2                 Less_y_2;
  typedef typename Traits::Signed_x_distance_2      Signed_x_distance_2;
  typedef typename Traits::Signed_y_distance_2      Signed_y_distance_2;
  typedef typename Traits::Infinity_distance_2      Infinity_distance_2;
  typedef typename Staircases< Traits >::Container  Container;
  typedef typename Staircases< Traits >::Iterator   Iterator;
  typedef typename Staircases< Traits >::Citerator  Citerator;
  typedef typename Staircases< Traits >::Intervall  Intervall;

  Infinity_distance_2 dist   = d.traits.infinity_distance_2_object();
  Less_x_2 lessx             = d.traits.less_x_2_object();
  Less_y_2 lessy             = d.traits.less_y_2_object();
  Signed_x_distance_2 sdistx = d.traits.signed_x_distance_2_object();
  Signed_y_distance_2 sdisty = d.traits.signed_y_distance_2_object();

  typename Traits::Construct_point_2_above_right_implicit_point_2
  cparip =
    d.traits.construct_point_2_above_right_implicit_point_2_object();
  typename Traits::Construct_point_2_above_left_implicit_point_2
  cpalip =
    d.traits.construct_point_2_above_left_implicit_point_2_object();
  typename Traits::Construct_point_2_below_right_implicit_point_2
  cpbrip =
    d.traits.construct_point_2_below_right_implicit_point_2_object();



  // test the four corners:
  for (int j = 0; j < 5; ++j) {
    const int k = j < 4 ? j : 3;

    // extract all points which are close enough to this point
    Point_2 corner = d[k];
    if (j >= 3) {
      if (j == 3) {
        Citerator i = d.tlstc_begin();
        while (sdistx(*i, d.minx) > FT(2) * d.r)
          ++i;
        corner = cpbrip(d.minx, *i, d.r);
      } else {
        Citerator i = d.tlstc_end();
        while (sdisty(d.maxy, *--i) > FT(2) * d.r) {}
        corner = cpbrip(*i, d.maxy, d.r);
      }
    }

    // find first point not covered by the rectangle at d[k]
    Iterator i = find_if(d.begin(), d.end(),
                         boost::bind(less<FT>(), d.r, boost::bind(dist, corner, _1)));

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

    CGAL_optimisation_expensive_assertion(
      save_end == find_if(d.end(), save_end,
                          boost::bind(less<FT>(), d.r, boost::bind(dist, corner, _1))));
    CGAL_optimisation_expensive_assertion(
      d.end() == find_if(d.begin(), d.end(),
                         boost::bind(std::greater_equal<FT>(),
                              d.r,
                              boost::bind(dist, corner, _1))));


    three_cover_points(d, o, ok);

    // restore saved sides of d:
    d.extreme(k) = save_side1;
    d.extreme((k+1) % 4) = save_side2;

    if (ok) {
      // does any rectangle contain the corner?
      if (d.end() != save_end) {
        *o++ = corner;
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
    Container share;
    d.shared_intervall(back_inserter(share));
    Citerator shf = share.end();
    Citerator shl = share.end();

    Citerator tl = d.tlstc_begin();
    Citerator lb = d.lbstc_begin();
    Citerator br = d.brstc_begin();
    Citerator rt = d.rtstc_begin();

    // make sure the top intervall is covered (left endpoint)
    // (it might be that top_i.first determines the placement of
    //  the top square)
    Point_2 top = top_i.first;
    if (lessx(top, *tl))
      while (++tl != d.tlstc_end() && !lessx(*tl, top)) {}
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

          // check the shared x-intervall
          if (!share.empty() && d.is_x_greater_y()) {
            // compute position of top in share
    #ifndef _MSC_VER
            while (shf != share.begin() && !lessx(*(shf - 1), top))
    #else
            while (shf != Citerator(share.begin()) &&
                   !lessx(*(shf - 1), top))
    #endif
              --shf;
    #ifndef _MSC_VER
            while (shl != share.begin() &&
    #else
            while (shl != Citerator(share.begin()) &&
    #endif
                   sdistx(*(shl - 1), top) > FT(2) * d.r)
              --shl;

            // make sure shared intervall is covered (left endpoint)
    #ifndef _MSC_VER
           if ((shf != share.begin() || shl == share.end()) &&
               lessx(share.front(), bottom))
             bottom = share.front();
           else if (shl != share.end() && lessx(*shl, bottom))
             bottom = *shl;
    #else
           if ((shf != Citerator(share.begin()) ||
                shl == Citerator(share.end())) &&
               lessx(share.front(), bottom))
             bottom = share.front();
           else if (shl != Citerator(share.end()) && lessx(*shl, bottom))
             bottom = *shl;
    #endif

          }


          // make sure the bottom and the shared intervall (right endpoint)
          // are covered
    #ifndef _MSC_VER
          if (sdistx(bottom_i.second, bottom) <= FT(2) * d.r &&
              (!d.is_x_greater_y() ||
               ((shl == share.end() ||
                 sdistx(share.back(), bottom) <= FT(2) * d.r) &&
                (shf == share.begin() ||
                 sdistx(*(shf - 1), bottom) <= FT(2) * d.r))))
    #else
          if (sdistx(bottom_i.second, bottom) <= FT(2) * d.r &&
              ((!d.is_x_greater_y() ||
               (shl == Citerator(share.end()) ||
                sdistx(share.back(), bottom) <= FT(2) * d.r) &&
               (shf == Citerator(share.begin()) ||
                sdistx(*(shf - 1), bottom) <= FT(2) * d.r))))
    #endif
            {
              // compute position of right square
              while (br != d.brstc_end() && sdistx(*br, bottom) <= FT(2) * d.r)
                ++br;
              // has the bottom square reached the bottom-right corner?
              if (br == d.brstc_end())
                break;
              Point_2 right = lessy(right_i.first, *br) ? right_i.first : *br;

              // check the shared y-intervall
              if (!share.empty() && !d.is_x_greater_y()) {
                // compute position of left in share
    #ifndef _MSC_VER
                while (shf != share.begin() &&
                       sdisty(left, *(shf - 1)) <= FT(2) * d.r)
    #else
                while (shf != Citerator(share.begin()) &&
                       sdisty(left, *(shf - 1)) <= FT(2) * d.r)
    #endif
                  --shf;
    #ifndef _MSC_VER
                while (shl != share.begin() &&
    #else
                while (shl != Citerator(share.begin()) &&
    #endif
                       lessy(left, *(shl - 1)))
                  --shl;

                // make sure shared intervall is covered (bottom endpoint)
    #ifndef _MSC_VER
                if ((shf != share.begin() || shl == share.end()) &&
                    lessy(share.front(), right))
                  right = share.front();
                else if (shl != share.end() && lessy(*shl, right))
                  right = *shl;
    #else
                if ((shf != Citerator(share.begin()) ||
                     shl == Citerator(share.end())) &&
                    lessy(share.front(), right))
                  right = share.front();
                else if (shl != Citerator(share.end())  && lessy(*shl, right))
                  right = *shl;
    #endif

              }


              // make sure the right intervall and the shared intervall
              // (top endpoint) are covered
    #ifndef _MSC_VER
              if (sdisty(right_i.second, right) <= FT(2) * d.r &&
                  (d.is_x_greater_y() ||
                   ((shl == share.end() ||
                     sdisty(share.back(), right) <= FT(2) * d.r) &&
                    (shf == share.begin() ||
                     sdisty(*(shf - 1), right) <= FT(2) * d.r))))
    #else
              if (sdisty(right_i.second, right) <= FT(2) * d.r &&
                  (d.is_x_greater_y() ||
                   ((shl == Citerator(share.end()) ||
                     sdisty(share.back(), right) <= FT(2) * d.r) &&
                    (shf == Citerator(share.begin()) ||
                     sdisty(*(shf - 1), right) <= FT(2) * d.r))))
    #endif
                {
                  // compute right bound for top square
                  while (rt != d.rtstc_end() &&
                         sdisty(*rt, right) <= FT(2) * d.r)
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

            } // if (bottom and shared intervall are covered)

        } // if (sdisty(left, left_i.first) <= FT(2) * d.r)

        top = *tl;
        if (!lessy(left_i.second, *tl) || ++tl == d.tlstc_end() ||
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

} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // ! (CGAL_PIERCE_RECTANGLES_2_H)
