// Copyright (c) 2013-06  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : 

#ifndef CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H
#define CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>


#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace internal {




} /* namespace internal */

/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

template <typename ConcurrencyTag,
          typename PointRange,
          typename NamedParameters = parameters::Default_named_parameters
>
void
constraint_based_smooth_point_set(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np = parameters::default_values())
{

}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H
