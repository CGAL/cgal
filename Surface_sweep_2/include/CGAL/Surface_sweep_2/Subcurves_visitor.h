// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Baruch Zukerman <baruchzu@post.tau.ac.il>
//            Efi Fogel       <efif@post.tau.ac.il>
//            (based on old version by Tali Zvi)

#ifndef CGAL_SURFACE_SWEEP_2_SUBCURVES_VISITOR_H
#define CGAL_SURFACE_SWEEP_2_SUBCURVES_VISITOR_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Definition of a surface-sweep visitor that reports all maximal x-monotone
 * non-intersecting x-monotone curves induced by a set of input curves.
 */

#include <vector>

#include <CGAL/Surface_sweep_2/Default_visitor.h>
#include <CGAL/Surface_sweep_2/Surface_sweep_2_utils.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class Subcurves_visitor
 *
 * A simple sweep-line visitor that reports all maximal x-monotone
 * non-intersecting x-monotone curves induced by a set of input curves. Used by
 * compute_subcurves().
 */
template <typename GeometryTraits_2, typename OutputIerator,
          typename Allocator_ = CGAL_ALLOCATOR(int)>
class Subcurves_visitor :
  public Default_visitor<Subcurves_visitor<GeometryTraits_2, OutputIerator,
                                           Allocator_>,
                         GeometryTraits_2, Allocator_>
{
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef OutputIerator                                 Output_ierator;
  typedef Allocator_                                    Allocator;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Subcurves_visitor<Gt2, Output_ierator, Allocator>
                                                        Self;
  typedef Default_visitor<Self, Gt2, Allocator>         Base;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  typedef typename Base::Event                          Event;
  typedef typename Base::Subcurve                       Subcurve;

  typedef typename Subcurve::Status_line_iterator       Status_line_iterator;

  typedef typename Base::Surface_sweep_2                Surface_sweep_2;

protected:
  // Data members:
  Output_ierator m_out;           // The output curves.
  bool m_overlapping;             // Should we report overlapping curves twice.

public:
  Subcurves_visitor(Output_ierator out, bool overlapping) :
    m_out(out),
    m_overlapping(overlapping)
  {}

  template <typename CurveIterator>
  void sweep(CurveIterator begin, CurveIterator end)
  {
    std::vector<X_monotone_curve_2> curves_vec;
    std::vector<Point_2> points_vec;

    curves_vec.reserve(std::distance(begin,end));
    make_x_monotone(begin, end,
                    std::back_inserter(curves_vec),
                    std::back_inserter(points_vec),
                    this->traits());

    // Perform the sweep.
    Surface_sweep_2* sl = this->surface_sweep();
    sl->sweep(curves_vec.begin(), curves_vec.end(),
              points_vec.begin(), points_vec.end());
  }

  void add_subcurve(const X_monotone_curve_2& cv, Subcurve *sc)
  {
    if (!m_overlapping) {
      // Report the curve once, whether it represents an overlap or not.
      *m_out = cv;
      ++m_out;
    }
    else {
      unsigned int overlap_depth = sc->overlap_depth();
      for (unsigned int i = 0; i < overlap_depth; ++i) {
        *m_out = cv;
        ++m_out;
      }
    }
  }

  Output_ierator output_iterator() { return m_out; }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
