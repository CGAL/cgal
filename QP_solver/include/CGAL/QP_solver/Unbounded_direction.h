// Copyright (c) 1997-2007  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Schoenherr
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp
//                 Kaspar Fischer

#ifndef CGAL_QP_SOLVER_UNBOUNDED_DIRECTION_H
#define CGAL_QP_SOLVER_UNBOUNDED_DIRECTION_H

#include <CGAL/license/QP_solver.h>


namespace CGAL {
template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::unbounded_direction_value(int i) const
{
 if (is_basic(i)) {                  // basic variable?
   return direction == 1 ? -q_x_O[in_B[i]] : q_x_O[in_B[i]];
 } else {                            // non-basic variable?
   if (i == j)                       // most recent entering variable?
     return direction == 1 ? d : -d;
   return et0;
 }
}
} //namespace CGAL

#endif // CGAL_QP_SOLVER_UNBOUNDED_DIRECTION_H

// ===== EOF ==================================================================
