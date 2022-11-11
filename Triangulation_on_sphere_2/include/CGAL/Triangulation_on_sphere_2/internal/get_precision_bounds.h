// Copyright (c) 1997-2021 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//                 Sébastien Loriot

#ifndef CGAL_TRIANGULATION_ON_SPHERE_GET_PRECISION_BOUNDS_H
#define CGAL_TRIANGULATION_ON_SPHERE_GET_PRECISION_BOUNDS_H

#include <CGAL/license/Triangulation_on_sphere_2.h>

#include <CGAL/double.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/number_utils.h>

#include <type_traits>

namespace CGAL {
namespace Triangulations_on_sphere_2 {
namespace internal {

// @todo could do something more suble than requiring exact SQRT representation (Root_of_2 etc.)
template <typename FT,
          bool has_exact_rep =
            is_same_or_derived<
              Field_with_sqrt_tag,
              typename Algebraic_structure_traits<FT>::Algebraic_category>::value &&
            !std::is_floating_point<FT>::value>
struct ToS2_precision_bound
{
  ToS2_precision_bound(const FT /*radius*/) { }

  // If you get a link error here, your number type is not supported and needs to implement its own specific values.
  // See the manual of the package for more information: https://doc.cgal.org/latest/Triangulation_on_sphere_2/
  FT get_squared_min_dist() const;
  FT get_squared_min_radius() const;
  FT get_squared_max_radius() const;
};

template <typename FT>
struct ToS2_precision_bound<FT, true> // exact representation of points on the sphere
{
  ToS2_precision_bound(const FT radius) : _sq_radius(CGAL::square(radius)) { }

  FT get_squared_min_dist() const { return 0; }
  FT get_squared_min_radius() const { return _sq_radius; }
  FT get_squared_max_radius() const { return _sq_radius; }

private:
  FT _sq_radius;
};

template <>
struct ToS2_precision_bound<double, false>
{
  ToS2_precision_bound(const double radius) : _radius(radius) { }

  double get_squared_min_dist() const { return CGAL::square(_radius * std::pow(2, -23)); }
  double get_squared_min_radius() const { return CGAL::square(_radius * (1 - std::pow(2, -50))); }
  double get_squared_max_radius() const { return CGAL::square(_radius * (1 + std::pow(2, -50))); }

private:
  double _radius;
};

} // namespace internal
} // namespace Triangulations_on_sphere_2
} // namespace CGAL

#endif // CGAL_TRIANGULATION_ON_SPHERE_GET_PRECISION_BOUNDS_H
