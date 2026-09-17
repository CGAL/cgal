// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_DIRECTIONAL_NON_CACHING_SEGMENT_BASIC_TRAITS_2_H
#define CGAL_ARR_DIRECTIONAL_NON_CACHING_SEGMENT_BASIC_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Arr_non_caching_segment_basic_traits_2.h>

namespace CGAL {

/*! \class
 * A model of the following concepts:
 * 1. AosBasicTraits_2,
 * 2. AosDirectionalXMonotoneTraits_2,
 * 4. AosConstructXMonotoneCurveTraits_2, and
 * 3. AosOpenBoundaryTraits_2
 * It handles linear curves.
 */
template <class Kernel_T>
class Arr_directional_non_caching_segment_basic_traits_2 :
  public Arr_non_caching_segment_basic_traits_2<Kernel_T> {
public:
  using Kernel = Kernel_T;

  using Base = Arr_non_caching_segment_basic_traits_2<Kernel>;
  using Segment_assertions = typename Base::Segment_assertions;
  using Has_exact_division = typename Base::Has_exact_division;

  /*! constructs default. */
  Arr_directional_non_caching_segment_basic_traits_2() : Base() {}

  /// \name Types and functors inherited from the base, required by the
  // AosBasicTraits_2 concept.
  //@{

  // Traits types:
  using Has_left_category = typename Base::Has_left_category;

  using Left_side_category = typename Base::Left_side_category;
  using Bottom_side_category = typename Base::Bottom_side_category;
  using Top_side_category = typename Base::Top_side_category;
  using Right_side_category = typename Base::Right_side_category;

  using Point_2 = typename Base::Point_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  /*! Compare the x-coordinates of two points */
  using Compare_x_2 = typename Base::Compare_x_2;

  /*! Compare two points lexigoraphically; by x, then by y */
  using Compare_xy_2 = typename Base::Compare_xy_2;

  /*! Obtain the left endpoint of a given segment */
  using Construct_min_vertex_2 = typename Base::Construct_min_vertex_2;

  /*! Obtain the right endpoint of a given segment */
  using Construct_max_vertex_2 = typename Base::Construct_max_vertex_2;

  /*! Check whether a given segment is vertical */
  using Is_vertical_2 = typename Base::Is_vertical_2;

  /*! Return the location of a given point with respect to an input segment */
  using Compare_y_at_x_2 = typename Base::Compare_y_at_x_2;

  /*! Check if two segments or if two points are identical */
  using Equal_2 = typename Base::Equal_2;

  /*! Compare the y value of two segments immediately to the left of their intersection point */
  using Compare_y_at_x_left_2 = typename Base::Compare_y_at_x_left_2;

  /*! Compare the y value of two segments immediately to the right of their intersection point */
  using Compare_y_at_x_right_2 = typename Base::Compare_y_at_x_right_2;

  /*! Construct a segment. */
  using Construct_x_monotone_curve_2 = typename Base::Construct_x_monotone_curve_2;

  //@}

  /// \name Types and functors introduced here, required by the
  // ArrangementDirectionalXMonotoneTraits concept.
  //@{
  using Construct_opposite_2 = typename Kernel::Construct_opposite_segment_2;

  /*! obtains a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const { return Construct_opposite_2(); }

  class Compare_endpoints_xy_2 {
  protected:
    typedef Arr_directional_non_caching_segment_basic_traits_2<Kernel> Traits;

    /*! The traits (in case it has state). */
    const Traits& m_traits;

    /*! Constructs
     * \param traits the traits (in case it has state).
     */
    Compare_endpoints_xy_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_directional_non_caching_segment_basic_traits_2<Kernel>;

  public:
    /*! compares the two endpoints of a given curve lexigoraphically.
     * \param cv The curve.
     * \return SMALLER if cv is directed from left to right and LARGER
     * otherwise.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv) const {
      using Construct_vertex_2 = typename Kernel::Construct_vertex_2;

      const Kernel& kernel = m_traits;
      Construct_vertex_2 ctr_v = kernel.construct_vertex_2_object();
      Compare_xy_2 cmp_xy = m_traits.compare_xy_2_object();
      return(cmp_xy(ctr_v(cv,0), ctr_v(cv,1)));
    }
  };

  /*! obtains a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const { return Compare_endpoints_xy_2(*this); }
  //@}
};

}

#include <CGAL/enable_warnings.h>

#endif
