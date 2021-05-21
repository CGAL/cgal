// Copyright (c) 2020 GeometryFactory SARL (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_WEIGHTS_INTERNAL_PROJECTION_TRAITS_BASE_3_H
#define CGAL_WEIGHTS_INTERNAL_PROJECTION_TRAITS_BASE_3_H

// #include <CGAL/license/Weights.h>

namespace CGAL {
namespace Weights {
namespace internal {
namespace ProjectionTraitsCartesianFunctors {

template <class Traits>
class Compare_along_axis {

  using Vector_3 = typename Traits::Vector_3;
  using Point = typename Traits::Point_2;
  Vector_3 base;

public:
  Compare_along_axis(const Vector_3& base) : base(base) { }
  using result_type = Comparison_result;
  Comparison_result operator() (const Point &p, const Point &q) const {
    return compare(base * (p - q), 0);
  }
};

template <class Traits>
class Less_along_axis {

  using Vector_3 = typename Traits::Vector_3;
  using Point = typename Traits::Point_2;
  Vector_3 base;

public:
  Less_along_axis(const Vector_3& base) : base(base) { }
  using result_type = bool;
  bool operator() (const Point &p, const Point &q) const {
    return base * (p - q) < 0;
  }
};

template<class Traits>
class Less_xy_along_axis {

  using Vector_3 = typename Traits::Vector_3;
  using Point = typename Traits::Point_2;
  Vector_3 base1, base2;

public:
  Less_xy_along_axis(const Vector_3& base1, const Vector_3& base2) :
  base1(base1), base2(base2)
  { }

  using result_type = bool;
  bool operator() (const Point &p, const Point &q) const {

    const Compare_along_axis<Traits> cx(base1);
    const Comparison_result crx = cx(p, q);
    if (crx == SMALLER) {
      return true;
    }
    if (crx == LARGER) {
      return false;
    }
    const Less_along_axis<Traits> ly(base2);
    return ly(p, q);
  }
};

template<class Traits>
class Equal_along_axis {

  using Vector_3 = typename Traits::Vector_3;
  using Point = typename Traits::Point_2;
  Vector_3 base1, base2;

public:
  Equal_along_axis(const Vector_3& base1, const Vector_3& base2) :
  base1(base1), base2(base2)
  { }

  using result_type = bool;
  bool operator() (const Point &p, const Point &q) const {

    const Compare_along_axis<Traits> cx(base1);
    const Compare_along_axis<Traits> cy(base2);

    const Comparison_result crx = cx(p, q);
    const Comparison_result cry = cy(p, q);

    if (crx == EQUAL && cry == EQUAL) {
      return true;
    }
    return false;
  }
};

template<class Traits>
class Projected_orientation_with_normal_3 {

  using K        = typename Traits::K;
  using Point    = typename Traits::Point_2;
  using Vector_3 = typename Traits::Vector_3;
  Vector_3 normal;

public:
  using Orientation = typename K::Orientation;
  using result_type = Orientation;

  Projected_orientation_with_normal_3(const Vector_3& normal_) : normal(normal_) { }
  Orientation operator()(const Point& p, const Point& q, const Point& r) const {
    return orientation(q - p, r - p, normal);
  }
};

template<class Traits>
class Projected_collinear_with_normal_3 {

  using K        = typename Traits::K;
  using Point    = typename Traits::Point_2;
  using Vector_3 = typename Traits::Vector_3;
  Vector_3 normal;

public:
  using Collinear   = typename K::Collinear_3;
  using result_type = bool;

  Projected_collinear_with_normal_3(const Vector_3& normal_) : normal(normal_) { }
  bool operator()(const Point& p, const Point& q, const Point& r) const {
    return collinear(p, q, r);
  }
};

} // namespace ProjectionTraitsCartesianFunctors

template<class Kernel>
class Projection_traits_base_3 {

  using Self = Projection_traits_base_3<Kernel>;
  typename Kernel::Vector_3 n, b1, b2;

public:
  using K        = Kernel;
  using FT       = typename K::FT;
  using Point_2  = typename K::Point_3;
  using Vector_2 = typename K::Vector_3;
  using Vector_3 = typename K::Vector_3;

  explicit Projection_traits_base_3(const Vector_3& n_) : n(n_) {

    CGAL_precondition(n != Vector_3(0, 0, 0));
    const FT nx = n.x();
    const FT ny = n.y();
    const FT nz = n.z();
    if (CGAL::abs(nz) >= CGAL::abs(ny)) {
      b1 = Vector_3(nz, FT(0), -nx);
    } else {
      b1 = Vector_3(ny, -nx, FT(0));
    }
    b2 = cross_product(n, b1);
  }

  const Vector_3& normal() const {
    return n;
  }

  const Vector_3& base1() const {
    return b1;
  }

  const Vector_3& base2() const {
    return b2;
  }

  using Less_xy_2
    = ProjectionTraitsCartesianFunctors::Less_xy_along_axis<Self>;

  using Equal_2
    = ProjectionTraitsCartesianFunctors::Equal_along_axis<Self>;

  using Orientation_2
    = ProjectionTraitsCartesianFunctors::Projected_orientation_with_normal_3<Self>;

  using Collinear_2
    = ProjectionTraitsCartesianFunctors::Projected_collinear_with_normal_3<Self>;

  using Compute_squared_distance_2
    = typename K::Compute_squared_distance_3;

  using Compute_squared_length_2
    = typename K::Compute_squared_length_3;

  using Compute_scalar_product_2
    = typename K::Compute_scalar_product_3;

  using Compute_area_2
    = typename K::Compute_area_3;

  using Construct_vector_2
    = typename K::Construct_vector_3;

  using Construct_circumcenter_2
    = typename K::Construct_circumcenter_3;

  Less_xy_2
  less_xy_2_object() const {
    return Less_xy_2(this->base1(), this->base2());
  }

  Equal_2
  equal_2_object() const {
    return Equal_2(this->base1(), this->base2());
  }

  Orientation_2
  orientation_2_object() const {
    return Orientation_2(this->normal());
  }

  Collinear_2
  collinear_2_object() const {
    return Collinear_2(this->normal());
  }

  Compute_squared_distance_2
  compute_squared_distance_2_object() const {
    return Compute_squared_distance_2();
  }

  Compute_squared_length_2
  compute_squared_length_2_object() const {
    return Compute_squared_length_2();
  }

  Compute_scalar_product_2
  compute_scalar_product_2_object() const {
    return Compute_scalar_product_2();
  }

  Compute_area_2
  compute_area_2_object() const {
    return Compute_area_2();
  }

  Construct_vector_2
  construct_vector_2_object() const {
    return Construct_vector_2();
  }

  Construct_circumcenter_2
  construct_circumcenter_2_object() const {
    return Construct_circumcenter_2();
  }
};

} // namespace internal
} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHTS_INTERNAL_PROJECTION_TRAITS_BASE_3_H
