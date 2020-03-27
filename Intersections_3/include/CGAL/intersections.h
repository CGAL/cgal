// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman



#ifndef CGAL_INTERSECTIONS_H
#define CGAL_INTERSECTIONS_H

// intersection and do_intersect pick the Kernel of their arguments
// based on Kernel_traits and Intersection_traits and use the functors
// Kernel::Do_intersect_{2,3} or Kernel::Intersect_{2,3} and the
// functors invoke on of the overloads of
// CGAL::internal::intersection/do_intersect.

// To define a new intersection or do_intersect function a
// specialization of Intersection_traits_{2,3} must be present and
// overloads internal::intersection(A, B) and internal::intersection(B, A) must
// be provided.

// If only a do_intersect function without a corresponding
// intersection function is added CGAL::do_intersect must also be
// overloaded in addition to CGAL::internal::do_intersect.

#include <CGAL/intersection_2.h>
#include <CGAL/intersection_3.h>

#endif // CGAL_INTERSECTIONS_H
