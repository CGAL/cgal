// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_DIRECTIONAL_NON_CACHING_SEGMENT_BASIC_TRAITS_2_H
#define CGAL_ARR_DIRECTIONAL_NON_CACHING_SEGMENT_BASIC_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Arr_non_caching_segment_basic_traits_2.h>

namespace CGAL {

/*! \class
 * A model of the following concepts:
 * 1. ArrangementBasicTraits_2,
 * 2. ArrangementDirectionalXMonotoneTraits_2,
 * 4. ArrangementConstructXMonotoneCurveTraits_2, and
 * 3. ArrangementOpenBoundaryTraits_2
 * It handles linear curves.
 */
template <class Kernel_T>
class Arr_directional_non_caching_segment_basic_traits_2 :
  public Arr_non_caching_segment_basic_traits_2<Kernel_T>
{
public:
  typedef Kernel_T                                         Kernel;

  typedef Arr_non_caching_segment_basic_traits_2<Kernel>   Base;
  typedef typename Base::Segment_assertions                Segment_assertions;
  typedef typename Base::Has_exact_division                Has_exact_division;

  /*! Default constructor. */
  Arr_directional_non_caching_segment_basic_traits_2() : Base() {}

  /// \name Types and functors inherited from the base, required by the
  // ArrangementBasicTraits_2 concept.
  //@{

  // Traits types:
  typedef typename Base::Has_left_category         Has_left_category;
  typedef typename Base::Has_do_intersect_category Has_do_intersect_category;

  typedef typename Base::Left_side_category        Left_side_category;
  typedef typename Base::Bottom_side_category      Bottom_side_category;
  typedef typename Base::Top_side_category         Top_side_category;
  typedef typename Base::Right_side_category       Right_side_category;

  typedef typename Base::Point_2                   Point_2;
  typedef typename Base::X_monotone_curve_2        X_monotone_curve_2;

  /*! Compare the x-coordinates of two points */
  typedef typename Base::Compare_x_2             Compare_x_2;

  /*! Compare two points lexigoraphically; by x, then by y */
  typedef typename Base::Compare_xy_2            Compare_xy_2;

  /*! Obtain the left endpoint of a given segment */
  typedef typename Base::Construct_min_vertex_2  Construct_min_vertex_2;

  /*! Obtain the right endpoint of a given segment */
  typedef typename Base::Construct_max_vertex_2  Construct_max_vertex_2;

  /*! Check whether a given segment is vertical */
  typedef typename Base::Is_vertical_2           Is_vertical_2;

  /*! Return the location of a given point with respect to an input segment */
  typedef typename Base::Compare_y_at_x_2        Compare_y_at_x_2;

  /*! Check if two segments or if two points are identical */
  typedef typename Base::Equal_2                 Equal_2;

  /*! Compare the y value of two segments immediately to the left of their
   * intersection point
   */
  typedef typename Base::Compare_y_at_x_left_2   Compare_y_at_x_left_2;

  /*! Compare the y value of two segments immediately to the right of their
   * intersection point
   */
  typedef typename Base::Compare_y_at_x_right_2  Compare_y_at_x_right_2;

  /*! Construct a segment. */
  typedef typename Base::Construct_x_monotone_curve_2
                                                 Construct_x_monotone_curve_2;

  /*! Obtain an approximation of a point coordinate. */
  typedef typename Base::Approximate_number_type Approximate_number_type;

  //@}

  /// \name Types and functors introduced here, required by the
  // ArrangementDirectionalXMonotoneTraits concept.
  //@{
  typedef typename Kernel::Construct_opposite_segment_2  Construct_opposite_2;

  /*! Obtain a Construct_opposite_2 functor object. */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }

  class Compare_endpoints_xy_2 {
  protected:
    typedef Arr_directional_non_caching_segment_basic_traits_2<Kernel> Traits;

    /*! The traits (in case it has state). */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state).
     */
    Compare_endpoints_xy_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_directional_non_caching_segment_basic_traits_2<Kernel>;

  public:
    /*! Compare the two endpoints of a given curve lexigoraphically.
     * \param cv The curve.
     * \return SMALLER if cv is directed from left to right and LARGER
     * otherwise.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv) const
    {
      typedef typename Kernel::Construct_vertex_2     Construct_vertex_2;

      const Kernel& kernel = m_traits;
      Construct_vertex_2 ctr_v = kernel.construct_vertex_2_object();
      Compare_xy_2 cmp_xy = m_traits.compare_xy_2_object();
      return(cmp_xy(ctr_v(cv,0), ctr_v(cv,1)));
    }
  };

  /*! Obtain a Compare_endpoints_xy_2 functor object. */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(*this); }
  //@}
};

}

#include <CGAL/enable_warnings.h>

#endif
