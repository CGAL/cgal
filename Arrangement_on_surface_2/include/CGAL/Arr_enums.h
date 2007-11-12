// Copyright (c) 2006 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s): Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_ENUM_H
#define CGAL_ARR_ENUM_H

/*! \file
 * Definition of the enumeration types for the arrangement package.
 */

#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

/*! \enum
 * Selection of a curve end.
 */
enum Curve_end
{
  MIN_END,
  MAX_END
};

/*! \enum
 * Indicator whether a halfedge is directed from left to right (from the
 * xy-lexicographically smaller vertex to the larger one), or from right to
 * left.
 */
enum Halfedge_direction
{
  LEFT_TO_RIGHT = -1,
  RIGHT_TO_LEFT = 1
};

/*! \enum
 * The various boundary conditions.
 */
enum Boundary_type
{
  NO_BOUNDARY = 0,
  MINUS_INFINITY = -1,
  PLUS_INFINITY = 1,
  AFTER_DISCONTINUITY = -2,
  BEFORE_DISCONTINUITY = 2,
  AFTER_SINGULARITY = -3,
  BEFORE_SINGULARITY = 3
};

inline CGAL::Sign sign (Boundary_type bt)
{
  if (static_cast<int> (bt) < 0)
    return (CGAL::NEGATIVE);
  else if (static_cast<int> (bt) > 0)
    return (CGAL::POSITIVE);

  return (CGAL::ZERO);
}

/*! \enum
 * 
 */
enum Arr_halfedge_direction {
  ARR_LEFT_TO_RIGHT = -1,
  ARR_RIGHT_TO_LEFT = 1
};

/*! \enum
 * 
 */
enum Arr_curve_end {
  ARR_MIN_END,
  ARR_MAX_END
};

/*! \enum The various surface boundary types
 * For example
 * - The plain has unbounded boundaries.
 * - The sphere has 2 contraction points and one identification curve.
 */
enum Arr_boundary_type {
  ARR_UNBOUNDED = 0,
  ARR_BOUNDED,
  ARR_CONTRACTION,
  ARR_IDENTIFICATION,
  ARR_NUM_BOUNDARY_TYPES
};

/*! \enum The various surface parameter space options categorizing the
 * surface range according to the parameter domain.
 */
enum Arr_parameter_space {
  ARR_LEFT_BOUNDARY = 0,
  ARR_RIGHT_BOUNDARY,
  ARR_BOTTOM_BOUNDARY,
  ARR_TOP_BOUNDARY,
  ARR_INTERIOR
};

CGAL_END_NAMESPACE

#endif

