// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_SEGMENT_TRAITS_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_SEGMENT_TRAITS_1_H

#include <CGAL/Arrangement_on_curve_1/Ft_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Parameter_space.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Geometry traits for a 1D arrangement on a closed segment \f$[a, b]\f$.  Both
// sides are closed.  Points are represented as `Kernel::FT`.  Pair with
// `Bounded_topology_traits`.
template <typename Kernel_>
class Segment_traits_1 : public Ft_traits_1<Kernel_> {
  using Kernel = Kernel_;
  using Base = Ft_traits_1<Kernel_>;

public:
  using Point_1 = typename Base::Point_1;

  using Left_side_category = Closed_side_tag;
  using Right_side_category = Closed_side_tag;

  /*! constructs from the boundary points
   */
  Segment_traits_1(const Point_1& left_bound, const Point_1& right_bound) :
    m_left_bound(left_bound),
    m_right_bound(right_bound)
  {}

  const Point_1& left_bound() const { return m_left_bound; }
  const Point_1& right_bound() const { return m_right_bound; }

  /*! \class Parameter_space_in_x_1
   * A functor that obtains the location of a point in the interval parameter space.
   */
  class Parameter_space_in_x_1 {
    const Segment_traits_1& m_tr;

  public:
    explicit Parameter_space_in_x_1(const Segment_traits_1& tr) : m_tr(tr) {}

    /*! obtains the location of `p` in the interval parameter space;
     * if `p` coincides with the leftmost point, returns LEFT_BOUNDARY;
     * otherwise, if `p` coincides with the rightmost point, returns RIGHT_BOUNDARY;
     * otherwise, returns INTERIOR
     */
    Parameter_space operator()(const Point_1& p) const {
      if (p == m_tr.m_left_bound) return Parameter_space::LEFT_BOUNDARY;
      if (p == m_tr.m_right_bound) return Parameter_space::RIGHT_BOUNDARY;
      return Parameter_space::INTERIOR;
    }
  };

  /*! obtains a functor object of type `Parameter_space_in_x_1`.
   */
  Parameter_space_in_x_1 parameter_space_in_x_1_object() const { return Parameter_space_in_x_1(*this); }

private:
  Point_1 m_left_bound;         //! the left bound
  Point_1 m_right_bound;        //! the right bound
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
