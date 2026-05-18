// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Baruch Zukerman <baruchzu@post.tau.ac.il>
//             Efi Fogel       <efif@post.tau.ac.il>
//             (based on old version by Tali Zvi)

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
template <typename GeometryTraits_2, typename OutputIerator, typename Allocator_ = CGAL_ALLOCATOR(int)>
class Subcurves_visitor :
  public Default_visitor<Subcurves_visitor<GeometryTraits_2, OutputIerator, Allocator_>, GeometryTraits_2, Allocator_> {
public:
  using Geometry_traits_2 = GeometryTraits_2;
  using Output_ierator = OutputIerator;
  using Allocator = Allocator_;

private:
  using Gt2 = Geometry_traits_2;
  using Self = Subcurves_visitor<Gt2, Output_ierator, Allocator>;
  using Base = Default_visitor<Self, Gt2, Allocator>;

public:
  using X_monotone_curve_2 = typename Gt2::X_monotone_curve_2;
  using Point_2 = typename Gt2::Point_2;

  using Event = typename Base::Event;
  using Subcurve = typename Base::Subcurve;

  using Status_line_iterator = typename Subcurve::Status_line_iterator;

  using Surface_sweep_2 = typename Base::Surface_sweep_2;

protected:
  // Data members:
  Output_ierator m_out;           // the output curves.
  bool m_overlapping;             // should we report overlapping curves twice.

public:
  /*! constructs.
   */
  Subcurves_visitor(Output_ierator out, bool overlapping) :
    m_out(out),
    m_overlapping(overlapping)
  {}

  /*!
   */
  void add_subcurve(const X_monotone_curve_2& cv, Subcurve* sc) {
    if (! m_overlapping) {
      // Report the curve once, whether it represents an overlap or not.
      *m_out = cv;
      ++m_out;
      return;
    }
    std::size_t overlap_depth = sc->overlap_depth();
    for (std::size_t i = 0; i < overlap_depth; ++i) {
      *m_out = cv;
      ++m_out;
    }
  }

  /*!
   */
  Output_ierator output_iterator() { return m_out; }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
