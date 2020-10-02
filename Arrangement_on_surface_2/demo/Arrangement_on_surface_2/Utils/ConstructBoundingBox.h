// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#ifndef ARRANGEMENT_DEMO_CONSTRUCT_BOUNDING_BOX
#define ARRANGEMENT_DEMO_CONSTRUCT_BOUNDING_BOX

#include <CGAL/Bbox_2.h>

// bounding box utility for arrangements
// doesn't have to be exact, only good enough for rendering
template <typename Traits_>
class ConstructBoundingBox
{
public:
  using Traits = Traits_;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Point_2 = typename Traits::Point_2;

  CGAL::Bbox_2 operator()(const X_monotone_curve_2& curve);
  CGAL::Bbox_2 operator()(const Point_2& point);
};

#endif
