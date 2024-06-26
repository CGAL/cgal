// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   cgal_kernel.h
 * @author Gernot Walzl
 * @date   2011-11-10
 */

#ifndef CGAL_KERNEL_H
#define CGAL_KERNEL_H

#include "config.h"

#include <CGAL/double.h>
#include <CGAL/float.h>
#include <CGAL/number_utils.h>

#ifdef USE_CGAL

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

namespace CGAL {

// typedef Cartesian<double> CK;
typedef CGAL::Simple_cartesian<CGAL::Exact_rational> SC_ER;
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt EPECK_w_sqrt;

typedef EPECK K;

typedef K::FT             FT;

typedef K::Point_2        Point2;
typedef K::Vector_2       Vector2;
typedef K::Line_2         Line2;
typedef K::Segment_2      Segment2;

typedef K::Point_3        Point3;
typedef K::Vector_3       Vector3;
typedef K::Ray_3          Ray3;
typedef K::Line_3         Line3;
typedef K::Segment_3      Segment3;
typedef K::Triangle_3     Triangle3;
typedef K::Plane_3        Plane3;
typedef K::Sphere_3       Sphere3;

template <typename FT>
decltype(auto) sqrt_with_warning(const FT& v)
{
  return CGAL::approximate_sqrt(v);
}

template <typename FT>
FT disallowed_sqrt(const FT& v)
{
  CGAL_assertion(false);
  std::exit(1);
  return CGAL::approximate_sqrt(v);
}

} // namespace CGAL

#else /* USE_CGAL */

namespace CGAL {

typedef void              K;

typedef double            FT;

} // namespace CGAL

#endif  /* USE_CGAL */

#endif /* CGAL_KERNEL_H */
