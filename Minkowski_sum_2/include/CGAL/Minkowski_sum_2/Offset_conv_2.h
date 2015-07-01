// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein   <wein_r@yahoo.com>

#ifndef CGAL_OFFSET_CONV_H
#define CGAL_OFFSET_CONV_H

#include <CGAL/Minkowski_sum_2/Union_of_curve_cycles_2.h>
#include <list>

namespace CGAL {

/*! \class
 * A class for computing the offset of a given polygon by a given radius
 * by constructing a single convolution cycle and computing its interior.
 */
template <class Base_>
class Offset_by_convolution_2 : private Base_
{
private:

  typedef Base_                                          Base;

public:

  typedef typename Base::Basic_kernel                    Kernel;
  typedef typename Base::Basic_NT                        NT;
  typedef typename Base::Polygon_2                       Polygon_2;
  typedef typename Base::Polygon_with_holes_2            Polygon_with_holes_2;
  typedef typename Base::Offset_polygon_2                Offset_polygon_2;

private:

  typedef typename Base::Labeled_traits_2                Labeled_traits_2;
  typedef typename Base::Labeled_curve_2                 Labeled_curve_2;
  typedef std::list<Labeled_curve_2>                     Curves_list;

  typedef Union_of_curve_cycles_2<Labeled_traits_2,
                                  Offset_polygon_2>      Union_2;

  using Base::_offset_polygon;
public:

  /*! Constructor. */
  Offset_by_convolution_2 (const Base_& base) :
    Base (base)
  {}

  /*!
   * Compute the offset of a simple polygon by a given radius.
   * Note that as the input polygon may not be convex, its offset may not be
   * simply connected. The result is therefore represented as the outer
   * boundary of the Minkowski sum (which is always a simple offset polygon)
   * and a container of offset polygons, representing the holes in this "outer"
   * polygon.
   * \param pgn The polygon.
   * \param r The offset radius.
   * \param off_bound Output: The outer boundary of the offset polygon.
   * \param off_holes Output: An output iterator for the holes in the offset.
   * \pre The polygon is simple.
   * \return A past-the-end iterator for the holes container.
   */
  template <class OutputIterator>
  OutputIterator operator() (const Polygon_2& pgn,
                             const NT& r,
                             Offset_polygon_2& off_bound,
                             OutputIterator off_holes) const
  {
    CGAL_precondition (pgn.is_simple());

    // Compute the curves that form the single convolution cycle for the
    // given polygon.
    Curves_list                     cycle;

    _offset_polygon (pgn,
                     CGAL::COUNTERCLOCKWISE,
                     r,
                     1,                       // The ID of the single cycle.
                     std::back_inserter (cycle));

    // Compute the union of the cycles that represent the offset polygon.
    Union_2     unite;

    off_holes = unite (cycle.begin(), cycle.end(),
                       off_bound, off_holes);

    return (off_holes);
  }

  /*!
   * Compute the offset of a polygon with holes by a given radius.
   * The result is represented as the outer boundary of the Minkowski sum
   * (which is always a simple offset polygon) and a container of offset
   * polygons, representing the holes in this "outer" polygon.
   * \param pwh The polygon with holes.
   * \param r The offset radius.
   * \param off_bound Output: The outer boundary of the offset polygon.
   * \param off_holes Output: An output iterator for the holes in the offset.
   * \pre The polygon is bounded (has an outer boundary).
   * \return A past-the-end iterator for the holes container.
   */
  template <class OutputIterator>
  OutputIterator operator() (const Polygon_with_holes_2& pwh,
                             const NT& r,
                             Offset_polygon_2& off_bound,
                             OutputIterator off_holes) const
  {
    CGAL_precondition (! pwh.is_unbounded());

    // Compute the curves that form the convolution cycle for the polygon
    // that forms the outer boundary.
    Curves_list                     cycle;
    unsigned int                    cycle_id = 1;

    _offset_polygon (pwh.outer_boundary(),
                     CGAL::COUNTERCLOCKWISE,
                     r,
                     cycle_id,
                     std::back_inserter (cycle));

    // Go over the polygon holes and compute the convolution cycle for each
    // hole. Note that in this case we traverse the holes in clockwise
    // orientation.
    typename Polygon_with_holes_2::Hole_const_iterator  hoit;

    for (hoit = pwh.holes_begin(); hoit != pwh.holes_end(); ++hoit)
    {
      cycle_id++;
      _offset_polygon (*hoit,
                       CGAL::CLOCKWISE,
                       r,
                       cycle_id,
                       std::back_inserter (cycle));
    }

    // Compute the union of the cycles that represent the offset polygon.
    Union_2     unite;

    off_holes = unite (cycle.begin(), cycle.end(),
                       off_bound, off_holes);

    return (off_holes);
  }

  /*!
   * Compute the inset of a simple polygon by a given radius.
   * Note that as the input polygon may not be convex, its offset may not be
   * simply connected. The result is therefore represented as a sequence of
   * polygons (which may also be empty).
   * \param pgn The polygon.
   * \param r The inset radius.
   * \param oi Output: An output iterator for the inset polygons.
   * \pre The polygon is simple.
   * \return A past-the-end iterator for the polygons container.
   */
  template <class OutputIterator>
  OutputIterator inset (const Polygon_2& pgn,
                        const NT& r,
                        OutputIterator oi) const
  {
    CGAL_precondition (pgn.is_simple());

    // Compute the curves that form the single convolution cycle for the
    // given polygon. Note that we traverse the polygon in clockwise direction,
    // as we treat it as a hole in the unbounded plane.
    Curves_list                     cycle;

    _offset_polygon (pgn,
                     CGAL::CLOCKWISE,
                     r,
                     1,                       // The ID of the single cycle.
                     std::back_inserter (cycle));

    // Compute the union of the cycles that represent the offset polygon.
    Union_2             unite;

    oi = unite.inverse (cycle.begin(), cycle.end(),
                        oi);

    return (oi);
  }
};

} //namespace CGAL

#endif
