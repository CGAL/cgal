// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ron Wein   <wein_r@yahoo.com>

#ifndef CGAL_OFFSET_DECOMP_H
#define CGAL_OFFSET_DECOMP_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Minkowski_sum_2/Union_of_curve_cycles_2.h>
#include <list>

namespace CGAL {

/*! \class
 * A class for computing the offset of a given polygon by a given radius,
 * by decomposing the polygon into convex sub-polygons and computing the union
 * off their offsets.
 */
template <class Base_, class DecompStrategy_>
class Offset_by_decomposition_2 : private Base_
{
private:

  typedef Base_                                        Base;

  using Base::_offset_polygon;

public:

  typedef typename Base::Basic_kernel                  Kernel;
  typedef typename Base::Basic_NT                      NT;
  typedef typename Base::Polygon_2                     Polygon_2;
  typedef typename Base::Offset_polygon_2              Offset_polygon_2;
  typedef DecompStrategy_                              Decomposition_strategy;

private:

  typedef std::list<Polygon_2>                         Polygons_list;
  typedef typename Polygons_list::iterator             Polygons_iterator;

  typedef typename Base::Labeled_traits_2              Labeled_traits_2;
  typedef typename Base::Labeled_curve_2               Labeled_curve_2;
  typedef std::list<Labeled_curve_2>                   Curves_list;

  typedef Union_of_curve_cycles_2<Labeled_traits_2,
                                  Offset_polygon_2>    Union_2;

public:

  /*! Constructor. */
  Offset_by_decomposition_2 (const Base_& base) :
    Base (base)
  {}

  /*!
   * Compute the offset of a simple polygon by a given radius.
   * Note that as the input polygon may not be convex, its offset may not be
   * simply connected. The result is therefore represented as the outer
   * boundary of the Minkowski sum (which is always a simple offset polygon)
   * and a container of offset polygons, representing the holes in this "outer"
   * polygon.
   * \param traits Arrangement traits that can deal with line segments and
   *               circular arcs.
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

    // Decompose the input polygon into convex sub-polygons.
    Decomposition_strategy  decomp_strat;
    Polygons_list           sub_pgns;
    Polygons_iterator       iter;

    decomp_strat (pgn, std::back_inserter(sub_pgns));

    // Compute the offset of each polygon separately.
    Curves_list                     boundary_curves;
    unsigned int                    pgn_id = 1;

    for (iter = sub_pgns.begin(); iter != sub_pgns.end(); ++iter)
    {
      _offset_polygon (*iter,
                       CGAL::COUNTERCLOCKWISE,
                       r,
                       pgn_id,
                       std::back_inserter(boundary_curves));
      pgn_id++;
    }

    // Compute the union of the cycles that represent the offset polygon.
    Union_2     unite;

    off_holes = unite (boundary_curves.begin(), boundary_curves.end(),
                       off_bound, off_holes);

    return (off_holes);
  }

};

} //namespace CGAL

#endif
