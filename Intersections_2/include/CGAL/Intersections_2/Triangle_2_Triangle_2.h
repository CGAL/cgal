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


#ifndef CGAL_INTERSECTIONS_2_TRIANGLE_2_TRIANGLE_2_H
#define CGAL_INTERSECTIONS_2_TRIANGLE_2_TRIANGLE_2_H

#include <CGAL/Intersections_2/internal/Triangle_2_Triangle_2_do_intersect_impl.h>
#include <CGAL/Intersections_2/internal/Triangle_2_Triangle_2_intersection_impl.h>

namespace CGAL {
CGAL_DO_INTERSECT_FUNCTION_SELF(Triangle_2, 2)
CGAL_INTERSECTION_FUNCTION_SELF(Triangle_2, 2)
} // namespace CGAL

#endif
