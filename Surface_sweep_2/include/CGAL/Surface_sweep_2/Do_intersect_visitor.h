// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Baruch Zukerman  <baruchzu@post.tau.ac.il>
//             Ron Wein         <wein@post.tau.ac.il>
//             Efi Fogel        <efif@post.tau.ac.il>

#ifndef CGAL_SURFACE_SWEEP_2_DO_INTERSECT_VISITOR_H
#define CGAL_SURFACE_SWEEP_2_DO_INTERSECT_VISITOR_H

#include <CGAL/license/Surface_sweep_2.h>

/*! \file
 *
 * Definition of the basic sweep-line visitors, for the usage of the global
 * sweep-line functions.
 */

#include <vector>
#include <type_traits>

#include <CGAL/Surface_sweep_2/Default_visitor.h>

namespace CGAL {
namespace Surface_sweep_2 {

/*! \class Do_intersect_visitor
 *
 * A plane-sweep visitor that determines whether the curves in a given set intersect.
 */
template <typename GeometryTraits_2, typename Allocator_ = CGAL_ALLOCATOR(int)>
class Do_intersect_visitor :
  public Default_visitor<Do_intersect_visitor<GeometryTraits_2, Allocator_>, GeometryTraits_2, Allocator_> {
protected:
  bool m_found_x;               // have we found an intersection so far.

public:
  /*! constructs.
   */
  Do_intersect_visitor() : m_found_x(false) {}

  /*! destructs.
   */
  virtual ~Do_intersect_visitor() {}

  /*!
   */
  void found_intersection() {
    m_found_x = true;
    this->stop_sweep();
  }

  /*!
   */
  bool do_intersect() { return m_found_x; }
};

} // namespace Surface_sweep_2
} // namespace CGAL

#endif
