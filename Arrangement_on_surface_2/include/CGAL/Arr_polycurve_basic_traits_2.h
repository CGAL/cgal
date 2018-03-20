// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University(Israel).
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
//                 Ron Wein  <wein@post.tau.ac.il>
//                 Dror Atariah <dror.atariah@fu-berlin.de>
//                 Waqar Khan <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_ARR_POLYCURVE_BASIC_TRAITS_2_H
#define CGAL_ARR_POLYCURVE_BASIC_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The traits-class for the general piece-wise (polycurve) type of curves of the
 * arrangement package.
 */

#include <iterator>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_non_caching_segment_traits_2.h>
#include <CGAL/Arr_geometry_traits/Polycurve_2.h>
#include <CGAL/Arr_geometry_traits/IO/Polycurve_2_iostream.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>

namespace CGAL {

template <typename SubcurveTraits_2 = Arr_non_caching_segment_traits_2<> >
class Arr_polycurve_basic_traits_2 {
public:
  typedef SubcurveTraits_2                                 Subcurve_traits_2;

  /// \name Types and functors inherited from the subcurve geometry traits.
  //@{

  typedef typename Subcurve_traits_2::Has_left_category    Has_left_category;
  typedef typename Subcurve_traits_2::Has_do_intersect_category
    Has_do_intersect_category;

  typedef typename Subcurve_traits_2::Left_side_category   Left_side_category;
  typedef typename Subcurve_traits_2::Bottom_side_category Bottom_side_category;
  typedef typename Subcurve_traits_2::Top_side_category    Top_side_category;
  typedef typename Subcurve_traits_2::Right_side_category  Right_side_category;

  typedef typename Arr_are_all_sides_oblivious_tag<Left_side_category,
                                                   Bottom_side_category,
                                                   Top_side_category,
                                                   Right_side_category>::result
    Are_all_sides_oblivious_tag;

  typedef typename Subcurve_traits_2::Point_2              Point_2;
  typedef typename Subcurve_traits_2::X_monotone_curve_2   X_monotone_subcurve_2;

  //@}

  // Backward compatibility:
  typedef X_monotone_subcurve_2                            X_monotone_segment_2;

private:
  typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>  Self;

  // Data members:
  const Subcurve_traits_2* m_subcurve_traits;  // The base segment-traits class.
  bool m_own_traits;

protected:
  enum { INVALID_INDEX = 0xffffffff };

public:
  /*! Construct default. */
  Arr_polycurve_basic_traits_2() :
    m_subcurve_traits(new Subcurve_traits_2()),
    m_own_traits(true)
  {}

  /*! Construct from a subcurve traits.
   * \param seg_traits an already existing subcurve tarits which is passed will
   *        be used by the class.
   */
  Arr_polycurve_basic_traits_2(const Subcurve_traits_2* geom_traits) :
    m_subcurve_traits(geom_traits), m_own_traits(false) {}

  /*! Construct copy.
   * If the 'other' polycurve traits owns its subcurve traits, then make
   * this polycurve traits own its subcurve traits as well
   * \param other the other traits.
   */
  Arr_polycurve_basic_traits_2(const Arr_polycurve_basic_traits_2& other)
  {
    m_subcurve_traits = (other.m_own_traits) ?
      new Subcurve_traits_2() : other.m_subcurve_traits;
    m_own_traits = other.m_own_traits;
  }

  /* Destructor
   * Deletes the subcurve tarits class in case it was constructed during the
   * construction of this.
   */
  ~Arr_polycurve_basic_traits_2()
  { if (m_own_traits) delete m_subcurve_traits; }

  /*! Obtain the subcurve traits.
   * \return the subcurve traits.
   */
  const Subcurve_traits_2* subcurve_traits_2() const
  { return m_subcurve_traits; }

  /// \name Types and functors defined here, required by the
  // ArrangementBasicTraits concept.
  //@{

  /*! An x monotone polycurve represents a continuous piecewise-linear
   * curve which is either strongly x-monotone or vertical. Again,
   * the polycurve is without degenerated subcurves.
   */
  typedef internal::X_monotone_polycurve_2<X_monotone_subcurve_2, Point_2>
                                                           X_monotone_curve_2;
  typedef typename X_monotone_curve_2::Size                Size;
  typedef typename X_monotone_curve_2::size_type           size_type;


  /*! Compare the x-coordinates of two points. */
  class Compare_x_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Compare_x_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Compare the x-coordinates of two directional points.
     * \param p1 the first directional point.
     * \param p2 the second directional point.
     * \return SMALLER - x(p1) < x(p2);
     *         EQUAL   - x(p1) = x(p2);
     *         LARGER  - x(p1) > x(p2).
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { return m_poly_traits.subcurve_traits_2()->compare_x_2_object()(p1, p2); }

    /*! Compare two ends of x-monotone curves in x.
     * \param xs1 the first curve.
     * \param ce1 the curve-end indicator of the first x-monotone curve xs1:
     *            ARR_MIN_END - the minimal end of xs1 or
     *            ARR_MAX_END - the maximal end of xs1.
     * \param p2 the second curve end.
     */
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const Point_2& p2)
    { return operator()(xs1, ce1, p2, Are_all_sides_oblivious_tag()); }

    /*! Compare two ends of x-monotone curves in x.
     * \param xs1 the first curve.
     * \param ce1 the curve-end indicator of the first x-monotone curve xs1:
     *            ARR_MIN_END - the minimal end of xs1 or
     *            ARR_MAX_END - the maximal end of xs1.
     * \param xs2 the second curve.
     * \param ce2 the curve-end indicator of the second x-monoton curve xs2:
     *            ARR_MIN_END - the minimal end of xs2 or
     *            ARR_MAX_END - the maximal end of xs2.
     */
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2)
    { return operator()(xs1, ce1, xs2, ce2, Are_all_sides_oblivious_tag()); }

  private:
    // Oblivious implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const Point_2& p2,
                                 Arr_all_sides_oblivious_tag) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      const Point_2& p1 = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      return geom_traits->compare_x_2_object()(p1, p2);
    }

    // Boundary implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const Point_2& p2,
                                 Arr_not_all_sides_oblivious_tag) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Parameter_space_in_x_2 ps_x =
        geom_traits->parameter_space_in_x_2_object();
      const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);

      if (ps_x1 != ARR_INTERIOR) {
        if (ps_x1 == ARR_LEFT_BOUNDARY) return SMALLER;
        if (ps_x1 == ARR_RIGHT_BOUNDARY) return LARGER;
      }

      typename Subcurve_traits_2::Parameter_space_in_y_2 ps_y =
        geom_traits->parameter_space_in_y_2_object();
      const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);
      if (ps_y1 == ARR_INTERIOR) {
        const Point_2& p1 = (ce1 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs1) :
          geom_traits->construct_min_vertex_2_object()(xs1);
        return geom_traits->compare_x_2_object()(p1, p2);
      }
      typename Subcurve_traits_2::Compare_x_on_boundary_2 cmp_x_on_bnd =
        geom_traits->compare_x_on_boundary_2_object();
      return opposite(cmp_x_on_bnd(p2, xs1, ce1));
    }

    // Oblivious implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2,
                                 Arr_all_sides_oblivious_tag) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      const Point_2& p1 = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      const Point_2& p2 = (ce2 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs2) :
        geom_traits->construct_min_vertex_2_object()(xs2);
      return geom_traits->compare_x_2_object()(p1, p2);
    }

    // Boundary implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2,
                                 Arr_not_all_sides_oblivious_tag) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Parameter_space_in_x_2 ps_x =
        geom_traits->parameter_space_in_x_2_object();
      const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);
      const Arr_parameter_space ps_x2 = ps_x(xs2, ce2);

      if (ps_x1 != ps_x2) {
        if (ps_x1 == ARR_LEFT_BOUNDARY) return SMALLER;
        if (ps_x1 == ARR_RIGHT_BOUNDARY) return LARGER;
        if (ps_x2 == ARR_LEFT_BOUNDARY) return LARGER;
        if (ps_x2 == ARR_RIGHT_BOUNDARY) return SMALLER;
      }

      // ps_x1 == ps_x2
      if (ps_x1 != ARR_INTERIOR) return EQUAL;

      typename Subcurve_traits_2::Parameter_space_in_y_2 ps_y =
        geom_traits->parameter_space_in_y_2_object();
      const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);
      const Arr_parameter_space ps_y2 = ps_y(xs2, ce2);
      if (ps_y1 == ARR_INTERIOR) {
        const Point_2& p1 = (ce1 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs1) :
          geom_traits->construct_min_vertex_2_object()(xs1);
        if (ps_y2 == ARR_INTERIOR) {
          const Point_2& p2 = (ce2 == ARR_MAX_END) ?
            geom_traits->construct_max_vertex_2_object()(xs2) :
            geom_traits->construct_min_vertex_2_object()(xs2);
          return geom_traits->compare_x_2_object()(p1, p2);
        }
        typename Subcurve_traits_2::Compare_x_on_boundary_2 cmp_x_on_bnd =
          geom_traits->compare_x_on_boundary_2_object();
        return cmp_x_on_bnd(p1, xs2, ce2);
      }
      if (ps_y2 == ARR_INTERIOR) {
        const Point_2& p2 = (ce2 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs2) :
          geom_traits->construct_min_vertex_2_object()(xs2);
        typename Subcurve_traits_2::Compare_x_on_boundary_2 cmp_x_on_bnd =
          geom_traits->compare_x_on_boundary_2_object();
        return opposite(cmp_x_on_bnd(p2, xs1, ce1));
      }
      return geom_traits->compare_x_on_boundary_2_object()(xs1, ce1, xs2, ce2);
    }
  };

  /*! Obtain a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object() const
  { return Compare_x_2(*this); }

  /*! Compare two curve-ends or points lexigoraphically: by x, then by y. */
  class Compare_xy_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Compare_xy_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Compare two directional points lexigoraphically: by x, then by y.
     * \param p1 the first enpoint directional point.
     * \param p2 the second endpoint directional point.
     * \return SMALLER - x(p1) < x(p2);
     *         SMALLER - x(p1) = x(p2) and y(p1) < y(p2);
     *         EQUAL   - x(p1) = x(p2) and y(p1) = y(p2);
     *         LARGER  - x(p1) = x(p2) and y(p1) > y(p2);
     *         LARGER  - x(p1) > x(p2).
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { return m_poly_traits.subcurve_traits_2()->compare_xy_2_object()(p1, p2); }

    /*! Compare two ends of x-monotone curves lexicographically.
     * \param xs1 the first curve.
     * \param ce1 the curve-end indicator of the first x-monotone curve xs1:
     *            ARR_MIN_END - the minimal end of xs1 or
     *            ARR_MAX_END - the maximal end of xs1.
     * \param p2 the second curve end.
     */
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const Point_2& p2)
    { return operator()(xs1, ce1, p2, Are_all_sides_oblivious_tag()); }

    /*! Compare two ends of x-monotone curves lexicographically.
     * \param xs1 the first curve.
     * \param ce1 the curve-end indicator of the first x-monotone curve xs1:
     *            ARR_MIN_END - the minimal end of xs1 or
     *            ARR_MAX_END - the maximal end of xs1.
     * \param xs2 the second curve.
     * \param ce2 the curve-end indicator of the second x-monoton curve xs2:
     *            ARR_MIN_END - the minimal end of xs2 or
     *            ARR_MAX_END - the maximal end of xs2.
     */
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2)
    { return operator()(xs1, ce1, xs2, ce2, Are_all_sides_oblivious_tag()); }

  private:
    // Oblivious implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const Point_2& p2,
                                 Arr_all_sides_oblivious_tag) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      const Point_2& p1 = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      return geom_traits->compare_xy_2_object()(p1, p2);
    }

    // Boundary implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const Point_2& p2,
                                 Arr_not_all_sides_oblivious_tag) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Parameter_space_in_x_2 ps_x =
        geom_traits->parameter_space_in_x_2_object();
      typename Subcurve_traits_2::Parameter_space_in_y_2 ps_y =
        geom_traits->parameter_space_in_y_2_object();
      const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);
      const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);

      if (ps_x1 != ARR_INTERIOR) {
        if (ps_x1 == ARR_LEFT_BOUNDARY) return SMALLER;
        if (ps_x1 == ARR_RIGHT_BOUNDARY) return LARGER;
      }

      if (ps_y1 == ARR_INTERIOR) {
        const Point_2& p1 = (ce1 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs1) :
          geom_traits->construct_min_vertex_2_object()(xs1);
        return geom_traits->compare_xy_2_object()(p1, p2);
      }

      // EFEF: missing implementation for open boundary.
      typename Subcurve_traits_2::Compare_x_on_boundary_2 cmp_x_on_bnd =
        geom_traits->compare_x_on_boundary_2_object();
      Comparison_result res = opposite(cmp_x_on_bnd(p2, xs1, ce1));
      if (res != EQUAL) return res;
      if (ps_y1 == ARR_TOP_BOUNDARY) return LARGER;
      CGAL_assertion(ps_y1 == ARR_BOTTOM_BOUNDARY);
      return SMALLER;
    }

    // Oblivious implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2,
                                 Arr_all_sides_oblivious_tag) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      const Point_2& p1 = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      const Point_2& p2 = (ce2 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs2) :
        geom_traits->construct_min_vertex_2_object()(xs2);
      return geom_traits->compare_xy_2_object()(p1, p2);
    }

    // Boundary implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2,
                                 Arr_not_all_sides_oblivious_tag) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Parameter_space_in_x_2 ps_x =
        geom_traits->parameter_space_in_x_2_object();
      typename Subcurve_traits_2::Parameter_space_in_y_2 ps_y =
        geom_traits->parameter_space_in_y_2_object();
      const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);
      const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);
      const Arr_parameter_space ps_x2 = ps_x(xs2, ce2);
      const Arr_parameter_space ps_y2 = ps_y(xs2, ce2);

      if (ps_x1 != ps_x2) {
        if (ps_x1 == ARR_LEFT_BOUNDARY) return SMALLER;
        if (ps_x1 == ARR_RIGHT_BOUNDARY) return LARGER;
        if (ps_x2 == ARR_LEFT_BOUNDARY) return LARGER;
        if (ps_x2 == ARR_RIGHT_BOUNDARY) return SMALLER;
      }

      if ((ps_x1 == ARR_INTERIOR) && (ps_y1 == ARR_INTERIOR)) {
        const Point_2& p1 = (ce1 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs1) :
          geom_traits->construct_min_vertex_2_object()(xs1);
        // ps1 == ARR_INTERIOR

        if ((ps_x2 == ARR_INTERIOR) && (ps_y2 == ARR_INTERIOR)) {
          const Point_2& p2 = (ce2 == ARR_MAX_END) ?
            geom_traits->construct_max_vertex_2_object()(xs2) :
            geom_traits->construct_min_vertex_2_object()(xs2);

          // ps1 == ARR_INTERIOR
          // ps2 == ARR_INTERIOR
          return geom_traits->compare_xy_2_object()(p1, p2);
        }

        // The cases ps_x2 == ARR_{LEFT,RIGHT}_BOUNDARY are handled above

        // ps1 == ARR_INTERIOR
        // ps_x2 == ARR_INTERIOR
        // ps_y2 != ARR_INTERIOR
        CGAL_assertion(ps_x2 == ARR_INTERIOR);
        // EFEF: missing implementation for open boundary.
        typename Subcurve_traits_2::Compare_x_on_boundary_2 cmp_x_on_bnd =
          geom_traits->compare_x_on_boundary_2_object();
        Comparison_result res = cmp_x_on_bnd(p1, xs2, ce2);
        if (res != EQUAL) return res;
        if (ps_y2 == ARR_TOP_BOUNDARY) return SMALLER;
        CGAL_assertion(ps_y2 == ARR_BOTTOM_BOUNDARY);
        return LARGER;
      }

      // ps1 != ARR_INTERIOR
      if ((ps_x2 == ARR_INTERIOR) && (ps_y2 == ARR_INTERIOR)) {
        const Point_2& p2 = (ce2 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs2) :
          geom_traits->construct_min_vertex_2_object()(xs2);

        // The cases ps_x1 == ARR_{LEFT,RIGHT}_BOUNDARY are handled above

        // ps_x1 == ARR_INTERIOR
        // ps_y1 != ARR_INTERIOR
        // ps2 == ARR_INTERIOR
        CGAL_assertion(ps_x1 == ARR_INTERIOR);
        typename Subcurve_traits_2::Compare_x_on_boundary_2 cmp_x_on_bnd =
          geom_traits->compare_x_on_boundary_2_object();
        Comparison_result res = cmp_x_on_bnd(p2, xs1, ce1);
        if (res != EQUAL) return opposite(res);
        if (ps_y1 == ARR_TOP_BOUNDARY) return LARGER;
        CGAL_assertion(ps_y1 == ARR_BOTTOM_BOUNDARY);
        return SMALLER;
      }

      // ps1 != ARR_INTERIOR
      // ps2 != ARR_INTERIOR
      // ps_x1 == ps_x2
      if (ps_x1 == ARR_INTERIOR) {
        // ps_y1 != ARR_INTERIOR
        // ps_y2 != ARR_INTERIOR
        Comparison_result res =
          geom_traits->compare_x_on_boundary_2_object()(xs1, ce1, xs2, ce2);
        if (res != EQUAL) return res;
        if (ps_y1 == ps_y2) return EQUAL;
        return (ps_y1 == ARR_BOTTOM_BOUNDARY) ? SMALLER : LARGER;
      }

      CGAL_assertion(ce1 == ce2);
      const Point_2& p1 = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      const Point_2& p2 = (ce2 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs2) :
        geom_traits->construct_min_vertex_2_object()(xs2);
      return geom_traits->compare_y_on_boundary_2_object()(p1, p2);
    }
  };

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(*this); }

  class Construct_min_vertex_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /* Constructor. */
    Construct_min_vertex_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Obtain the left endpoint of the x-monotone polycurve.
     * \todo: is it possible to make the return type const reference if the
     *        return type of the subcurve traits is const reference?
     * \param cv The polycurve curve.
     * \return The left endpoint.
     */
    Point_2 operator()(const X_monotone_curve_2& cv) const
    {
      CGAL_assertion(cv.number_of_subcurves() > 0);

      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();

      if (geom_traits->compare_endpoints_xy_2_object()(cv[0]) == SMALLER)
        return geom_traits->construct_min_vertex_2_object()(cv[0]);
      else
        return geom_traits->
          construct_min_vertex_2_object()(cv[cv.number_of_subcurves()-1]);
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(*this); }

  class Construct_max_vertex_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Construct_max_vertex_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Obtain the right endpoint of the x-monotone polycurve.
     * \todo: is it possible to make the return type const reference if the
     *        return type of the subcurve traits is const reference?
     * \param cv The polycurve.
     * \return The right endpoint.
     */
    Point_2 operator()(const X_monotone_curve_2& cv) const
    {
      CGAL_assertion(cv.number_of_subcurves() > 0);

      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();

      if (geom_traits->compare_endpoints_xy_2_object()(cv[0]) == SMALLER)
        return geom_traits->
          construct_max_vertex_2_object()(cv[cv.number_of_subcurves()-1]);
      else
        return geom_traits->construct_max_vertex_2_object()(cv[0]);
    }
  };

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(*this); }

  class Is_vertical_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
    Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Is_vertical_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv) const
    {
      // An x-monotone polycurve can represent a vertical segment only if it
      // is comprised of vertical segments. If the first subcurve is vertical,
      // all subcurves are vertical in an x-monotone polycurve
      return m_poly_traits.subcurve_traits_2()->is_vertical_2_object()(cv[0]);
    }
  };

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(*this); }

  class Compare_y_at_x_2 {
  private:
    // Oblivious implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_all_sides_oblivious_tag) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      const Point_2& p = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      return geom_traits->compare_y_at_x_2_object()(p, xs2);
    }

    // Boundary implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_not_all_sides_oblivious_tag) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Parameter_space_in_x_2 ps_x =
        geom_traits->parameter_space_in_x_2_object();
      typename Subcurve_traits_2::Parameter_space_in_y_2 ps_y =
        geom_traits->parameter_space_in_y_2_object();
      typename Subcurve_traits_2::Construct_min_vertex_2 min_vertex =
        geom_traits->construct_min_vertex_2_object();
      typename Subcurve_traits_2::Construct_max_vertex_2 max_vertex =
        geom_traits->construct_max_vertex_2_object();

      const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);
      const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);

      CGAL_assertion(((ce1 == ARR_MAX_END) && (ps_x1 != ARR_LEFT_BOUNDARY)) ||
                     ((ce1 == ARR_MIN_END) && (ps_x1 != ARR_RIGHT_BOUNDARY)));

      if (ps_x1 == ARR_INTERIOR) {
        if (ps_y1 == ARR_TOP_BOUNDARY) {
          typename Subcurve_traits_2::Equal_2 equal =
            geom_traits->equal_2_object();
          const Point_2& p = (ce1 == ARR_MAX_END) ?
            max_vertex(xs1) : min_vertex(xs1);
          if (equal(p, max_vertex(xs2))) return EQUAL;
          if (equal(p, min_vertex(xs2))) return EQUAL;
          return LARGER;
        }
        if (ps_y1 == ARR_BOTTOM_BOUNDARY) {
          typename Subcurve_traits_2::Equal_2 equal =
            geom_traits->equal_2_object();
          const Point_2& p = (ce1 == ARR_MAX_END) ?
            max_vertex(xs1) : min_vertex(xs1);
          if (equal(p, max_vertex(xs2))) return EQUAL;
          if (equal(p, min_vertex(xs2))) return EQUAL;
          return SMALLER;
        }
        // ps_y1 == ARR_INTERIOR
        const Point_2& p = (ce1 == ARR_MAX_END) ?
          max_vertex(xs1) : min_vertex(xs1);
        return geom_traits->compare_y_at_x_2_object()(p, xs2);
      }
      // ps_x1 == ARR_RIGHT_BOUNDARY || ARR_LEFT_BOUNDARY
      const Point_2& p1 = (ce1 == ARR_MAX_END) ?
        max_vertex(xs1) : min_vertex(xs1);
      const Point_2& p2 = (ce1 == ARR_MAX_END) ?
        max_vertex(xs2) : min_vertex(xs2);
      return geom_traits->compare_y_on_boundary_2_object()(p1, p2);
    }

  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Compare_y_at_x_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits) {}

    /*! Obtain the location of the given point with respect to the input curve.
     * \param p The point.
     * \param xcv The polycurve curve.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) < cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      if (! m_poly_traits.is_vertical_2_object()(xcv)) {
        // Get the index of the subcurve in xcv containing p.
        std::size_t i =
          m_poly_traits.locate_impl(xcv, p, Are_all_sides_oblivious_tag());
        CGAL_precondition(i != INVALID_INDEX);

        // Compare the subcurve xcv[i] and p.
        return geom_traits->compare_y_at_x_2_object()(p, xcv[i]);
      }
      // The curve is vertical
      #ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      const Comparison_result SMLLR = SMALLER;
      const Comparison_result LRGR = LARGER;
      #else
      const bool is_left_to_right = m_poly_traits.subcurve_traits_2()->
                                      compare_endpoints_xy_2_object()(xcv[0])
                                        == SMALLER;
      const Comparison_result SMLLR = is_left_to_right?SMALLER:LARGER;
      const Comparison_result LRGR = is_left_to_right?LARGER:SMALLER;
      #endif

      Comparison_result rc = geom_traits->compare_y_at_x_2_object()(p, xcv[0]);
      if (rc == SMLLR) return SMLLR;
      std::size_t n = xcv.number_of_subcurves();
      rc = geom_traits->compare_y_at_x_2_object()(p, xcv[n-1]);
      if (rc == LRGR) return LRGR;
      return EQUAL;
    }

    /*! Obtain the location of the given curve_end with respect to the input
     * curve.
     * \param xcv The polycurve curve.
     * \param ce the curve-end indicator of the x-monotone subcurve xl:
     *            ARR_MIN_END - the minimal end of xl or
     *            ARR_MAX_END - the maximal end of xl.
     * \param xcv The polycurve curve.
     * \pre the curve-end is in the x-range of xcv.
     * \return SMALLER if if y(xs, ce) < cv(x(xs, ce)), i.e.  the curve-end
     *           is below the curve xcv;
     *         LARGER if y(xs, ce) > cv(x(xs, ce)), i.e. the curve-end is
     *           above the curve xcv;
     *         EQUAL if the curve-end lies on the curve xcv.
     */
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2) const
    { return operator()(xs1, ce1, xs2, Are_all_sides_oblivious_tag()); }
  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(*this); }

  class Compare_y_at_x_left_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Compare_y_at_x_left_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Compare the y value of two x-monotone curves immediately to the left
     * of their intersection point.
     * \param cv1 The first polycurve curve.
     * \param cv2 The second polycurve curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined(lexicographically) to its left.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& p) const
    {
      // Get the indices of the subcurves in cv1 and cv2 containing p and
      // defined to its left.
      std::size_t i1 = m_poly_traits.locate_side(cv1, p, false);
      std::size_t i2 = m_poly_traits.locate_side(cv2, p, false);

      CGAL_precondition(i1 != INVALID_INDEX);
      CGAL_precondition(i2 != INVALID_INDEX);

      // Compare cv1[i1] and cv2[i2] at p's left.
      return m_poly_traits.subcurve_traits_2()->
        compare_y_at_x_left_2_object()(cv1[i1], cv2[i2], p);
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(*this); }

  class Compare_y_at_x_right_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Compare_y_at_x_right_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Compare the y value of two x-monotone curves immediately to the right
     * of their intersection point.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined(lexicographically) to its right.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the right of p: SMALLER, LARGER or EQUAL.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& p) const
    {
      // Get the indices of the subcurves in cv1 and cv2 containing p and
      // defined to its right.
      std::size_t i1=m_poly_traits.locate_side(cv1, p, true);
      std::size_t i2=m_poly_traits.locate_side(cv2, p, true);

      CGAL_precondition(i1 != INVALID_INDEX);
      CGAL_precondition(i2 != INVALID_INDEX);

      // Compare cv1[i1] and cv2[i2] at p's right.
      return m_poly_traits.subcurve_traits_2()->
        compare_y_at_x_right_2_object()(cv1[i1], cv2[i2], p);
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 functor object.
   */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(*this); }

  class Equal_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Equal_2(const Polycurve_basic_traits_2& poly_tr) : m_poly_traits(poly_tr) {}

    /*! Check whether the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same;(false) otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    { return m_poly_traits.subcurve_traits_2()->equal_2_object()(p1, p2); }

    /*! Check whether the two x-monotone curves are the same (have the same
     * graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return(true) if the two curves are the same;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      // Check the pairwise equality of the contained subcurves.
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Equal_2 equal =
        geom_traits->equal_2_object();
      typename Subcurve_traits_2::Compare_x_2 compare_x =
        geom_traits->compare_x_2_object();
      typename Subcurve_traits_2::Compare_y_at_x_2 compare_y_at_x =
        geom_traits->compare_y_at_x_2_object();
      typename Subcurve_traits_2::Construct_max_vertex_2 max_vertex =
        geom_traits->construct_max_vertex_2_object();
      typename Subcurve_traits_2::Compare_endpoints_xy_2 comp_endpt =
        geom_traits->compare_endpoints_xy_2_object();
      Is_vertical_2 is_vertical = m_poly_traits.is_vertical_2_object();
      Construct_min_vertex_2 xpoly_min_v =
        m_poly_traits.construct_min_vertex_2_object();
      Construct_max_vertex_2 xpoly_max_v =
        m_poly_traits.construct_max_vertex_2_object();

      // The first and last points of the subcurves should be equal.
      bool res = equal(xpoly_min_v(cv1), xpoly_min_v(cv2));
      if (!res) return false;
      res = equal(xpoly_max_v(cv1), xpoly_max_v(cv2));
      if (!res) return false;

      // If the first and last points are equal and the curves are vertical,
      // it means that it is equal.
      bool ver1 = is_vertical(cv1);
      bool ver2 = is_vertical(cv2);
      // both curves are vertical and therefore equal.
      if (ver1 && ver2) return true;
      // one is vertical and the other is not - hence not equal.
      if (ver1 || ver2) return false;

      // If we arrived here it means that the first and last point of the
      // curve are equal.
      std::size_t i = 0;
      std::size_t j = 0;
      std::size_t n1 = cv1.number_of_subcurves();
      std::size_t n2 = cv2.number_of_subcurves();
      Comparison_result is_cv1_left_to_right = comp_endpt(cv1[0]);
      Comparison_result is_cv2_left_to_right = comp_endpt(cv2[0]);

      while ((i < (n1-1)) || (j < (n2-1))) {
        Point_2 point1, point2;
        std::size_t cv1_seg_ind, cv2_seg_ind;
        if (SMALLER == is_cv1_left_to_right) {
          cv1_seg_ind = i;
          point1 = max_vertex(cv1[cv1_seg_ind]);
        }
        else {
          cv1_seg_ind = n1 - 1 - i;
          point1 = max_vertex(cv1[cv1_seg_ind]);
        }
        if (SMALLER == is_cv2_left_to_right) {
          cv2_seg_ind = j;
          point2 = max_vertex(cv2[cv2_seg_ind]);
        }
        else {
          cv2_seg_ind = n2 - 1 - j;
          point2 = max_vertex(cv2[cv2_seg_ind]);
        }

        bool res = equal(point1, point2);
        // Easy case - the two points are equal
        if (res) {
          ++i;
          ++j;
        }
        else {
          Comparison_result res_x = compare_x(point1,point2);
          // Check if the different point is a collinear point situated on
          // the line between its two neighbors.
          if (SMALLER == res_x) {
            Comparison_result res_y_at_x =
              compare_y_at_x(point1, cv2[cv2_seg_ind]);
            if (EQUAL == res_y_at_x) ++i;
            else return false;
          }
          else if (LARGER == res_x) {
            Comparison_result res_y_at_x =
              compare_y_at_x(point2,cv1[cv1_seg_ind]);
            if (EQUAL == res_y_at_x) ++j;
            else return false;
          }
          else return false;
        }
      }
      return true;
    }
  };

  /*! Obtain an Equal_2 functor object. */
  Equal_2 equal_2_object() const { return Equal_2(*this); }

  class Compare_endpoints_xy_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Compare_endpoints_xy_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Compare the endpoints of an \(x\)-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param cv The curve.
     * \return SMALLER if the curve is oriented left-to-right;
     *         LARGER if the curve is oriented right-to-left.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv) const
    {
      return (m_poly_traits.subcurve_traits_2()->
              compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) ?
        (SMALLER) : (LARGER);
    }
  };
  //@}

  /// \name Types and functors defined here, required by the
  // ArrangementDirectionalXMonotoneTraits_2 concept.
  //@{

  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(*this); }

  class Construct_opposite_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor */
    Construct_opposite_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Construct the reversed \(x\)-monotone polycurve of the input.
     * Note that the functor constructs the opposites of _all_ subcurves
     * constituting xcv.
     * \param xcv the \(x\)-monotone polycurve to be reveres
     * \pre xcv contains at least one subcurve
     * \return An \(x\)-monotone polycurve with the same graph as the input xcv
     *         only with a reverse orientation.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Construct_opposite_2 const_op =
        geom_traits->construct_opposite_2_object();
      std::vector<X_monotone_subcurve_2> rev_segs(xcv.number_of_subcurves());;
      typename X_monotone_curve_2::Subcurve_const_iterator sit;
      typename X_monotone_curve_2::Subcurve_iterator tit = rev_segs.begin();
      for (sit = xcv.subcurves_begin(); sit != xcv.subcurves_end(); ++sit)
        *tit++ = const_op(*sit);
      return X_monotone_curve_2(rev_segs.rbegin(), rev_segs.rend());
    }
  };

  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(*this); }

  ///@}

  /// \name Types and functors defined here, required by the
  // ArrangementLandmarkTraits concept.
  //@{

#if 0
  // The following block assumes that the subcurve traits template parameter
  // is a model of the ArrangementLandmarkTraits concept; in other words, it
  // defines the nested types Approximate_number_type and Approximate_2 and
  // the member function approximate_2_object(). It cannot be used as is if
  // the subcurve traits does not model the ArrangementLandmarkTraits concept.
  // The functor Construct_x_monotone_curve_2 is provided regardless of the
  // subcurve traits.

  typedef typename Subcurve_traits_2::Approximate_number_type
  Approximate_number_type;
  typedef typename Subcurve_traits_2::Approximate_2    Approximate_2;

  /*! Obtain an Approximate_2 functor object. */
  Approximate_2 approximate_2_object() const
  { return subcurve_traits_2()->approximate_2_object(); }
#else
  // The following block defines the nested types Approximate_number_type and
  // Approximate_2 and the member function approximate_2_object() based on the
  // corresponding types and function definitions of the subcurve traits. If
  // the subcurve traits does not provide these definitions, they are defined
  // as dummies. Essentially, the polycurve traits becomes a practical model of
  // the ArrangementLandmarkTraits concept only if the subcurve traits is a
  // model of this concept.
  //
  // The following implementation is inspired by
  // http://stackoverflow.com/a/11816999/1915421

  template <typename T>
  struct Void {
    typedef void type;
  };

  template <typename T, typename _ = void>
  struct has_approximate_2 {
    // Generic implementation
    typedef void                        Approximate_number_type;
    typedef void                        Approximate_2;
  };

  template <typename T>
  struct has_approximate_2<T, typename Void<typename T::Approximate_2>::type>
  {
    // Specialization for types holding a nested type T::Approximate_2
    typedef typename T::Approximate_number_type
                                        Approximate_number_type;
    typedef typename T::Approximate_2   Approximate_2;
  };

  typedef typename has_approximate_2<Subcurve_traits_2>::Approximate_number_type
                                        Approximate_number_type;
  typedef typename has_approximate_2<Subcurve_traits_2>::Approximate_2
                                        Approximate_2;

  /*! Obtain an Approximate_2 functor object. */
  Approximate_2 approximate_2_object_impl(boost::false_type) const
  { return subcurve_traits_2()->approximate_2_object(); }

  Approximate_2 approximate_2_object_impl(boost::true_type) const { }

  Approximate_2 approximate_2_object() const
  {
    typedef typename boost::is_same<void, Approximate_2>::type      Is_void;
    return approximate_2_object_impl(Is_void());
  }
#endif

  class Construct_x_monotone_curve_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Construct_x_monotone_curve_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Obtain an x-monotone polycurve that consists of one given subcurve.
     * \param seg input subcurve.
     * \pre seg is not degenerated.
     * \return An x-monotone polycurve with one subcurve.
     */
    X_monotone_curve_2 operator()(const X_monotone_subcurve_2& seg) const
    {
      CGAL_precondition_code
        (
         /* Test that the subcurve is not degenerated. We do this test
          * independently from the subcurve traits in use, as we do not
          * allow a polycurve with degenerated subcurves.
          */
         const Subcurve_traits_2* geom_traits =
           m_poly_traits.subcurve_traits_2();
         typename Subcurve_traits_2::Construct_min_vertex_2 get_min_v =
           geom_traits->construct_min_vertex_2_object();
         typename Subcurve_traits_2::Construct_max_vertex_2 get_max_v =
           geom_traits->construct_max_vertex_2_object();
         typename Subcurve_traits_2::Equal_2 equal =
           geom_traits->equal_2_object();

         CGAL_precondition_msg(!equal(get_min_v(seg),get_max_v(seg)),
                               "Cannot construct a degenerated subcurve");
         );

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      if (m_poly_traits.subcurve_traits_2()->
          compare_endpoints_xy_2_object()(seg) == LARGER)
        return X_monotone_subcurve_2(m_poly_traits.subcurve_traits_2()->
                                    construct_opposite_2_object()(seg));
#endif

      return X_monotone_curve_2(seg);
    }

    /*! Construct an x-monotone polycurve which is well-oriented from a range of
     * elements.
     * \pre Range should from a continuous well-oriented x-monotone polycurve.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 operator()(ForwardIterator begin,
                                  ForwardIterator end) const
    {
      typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
      typedef typename boost::is_same<VT,Point_2>::type Is_point;

      // Dispatch the range to the appropriate implementation.
      return constructor_impl(begin, end, Is_point());
    }

    /*! Construct an x-monotone polycurve from a range of points.
     * The polycurve may be oriented left-to-right or right-to-left
     * depending on the lexicographical order of the points in the
     * input.
     * \pre Range contains at least two points.
     * \pre No two consecutive points are the same.
     * \pre The points form an continuous well-oriented x-monotone polycurve.
     * \post By the construction the returned polycurve is well-oriented.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 constructor_impl(ForwardIterator /* begin */,
                                        ForwardIterator /* end */,
                                        boost::true_type) const
    { CGAL_error_msg("Cannot construct a polycurve from a range of points!"); }

    /*! Obtain an x-monotone polycurve from a range of subcurves.
     * \param begin An iterator pointing to the first subcurve in the range.
     * \param end An iterator pointing to the past-the-end subcurve
     * in the range.
     * \pre The range contains at least one subcurve.
     * \pre Subcurves correspond to a well-oriented polycurve. That
     *      is, the target of the i-th subcurve is an source of the
     *      (i+1)th subcurve.
     * \pre The sequence of subcurves in the range forms a weak x-monotone
     *      polycurve.
     * \pre The container should support bidirectional iteration.
     * \return A continuous, well-oriented x-monotone polycurve which
     *         is directed either left-to-right or right-to-left
     *         depending on the subcurves in the input.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 constructor_impl(ForwardIterator begin,
                                        ForwardIterator end,
                                        boost::false_type) const
    {
      CGAL_precondition_msg
        (
         begin != end,
         "Input range of subcurves has to contain at least one subcurve"
         );

      CGAL_precondition_code
        (
         const Subcurve_traits_2* geom_traits =
           m_poly_traits.subcurve_traits_2();
         typename Subcurve_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
           geom_traits->compare_endpoints_xy_2_object();
         typename Subcurve_traits_2::Construct_min_vertex_2 get_min_v =
           geom_traits->construct_min_vertex_2_object();
         typename Subcurve_traits_2::Construct_max_vertex_2 get_max_v =
           geom_traits->construct_max_vertex_2_object();
         typename Subcurve_traits_2::Equal_2 equal =
           geom_traits->equal_2_object();

         ForwardIterator curr = begin;
         ForwardIterator next = begin;
         ++next;
         );


      CGAL_precondition_msg
        (
         (next != end) || !equal(get_max_v(*curr),get_min_v(*curr)),
         "Cannot construct a polycurve with degenerated subcurve"
         );

      CGAL_precondition_code
        (
         // Range contains at least two subcurves

         Comparison_result init_dir = cmp_seg_endpts(*curr);
         while (next != end){
           CGAL_precondition_msg
             (!equal(get_min_v(*next),get_max_v(*next)),
              "Cannot construct a polycurve with degenerated subcurve"
              );
           CGAL_precondition_msg
             (
              init_dir == cmp_seg_endpts(*next),
              "Subcurves must form x-monotone polycurve"
              );
           if (init_dir == SMALLER){
             CGAL_precondition_msg
               (
                equal(get_max_v(*curr),get_min_v(*next)),
                "Subcurves should concatenate in source->target manner"
                );
           }
           else{
             CGAL_precondition_msg
               (
                equal(get_min_v(*curr),get_max_v(*next)),
                "Subcurves should concatenate in source->target manner"
                );
           }
           ++curr;
           ++next;
         }
         );

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      if (m_poly_traits.subcurve_traits_2()->
          compare_endpoints_xy_2_object()(*begin) == LARGER)
      {
        X_monotone_curve_2 xcv(begin, end);
        return m_poly_traits.construct_opposite_2_object()(xcv);
      }
#endif

      return X_monotone_curve_2(begin, end);
    }
  };

  /*! Obtain a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(*this); }

  //@}

  /// \name Types and functors defined here, required by the
  // ArrangementOpenBoundaryTraits_2 concept.
  //@{

  /*! A function object that obtains the parameter space of a geometric
   * entity along the x-axis
   */
  class Parameter_space_in_x_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Parameter_space_in_x_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Obtains the parameter space at the end of a curve along the x-axis .
     * Note that if the curve-end coincides with a pole, then unless the curve
     * coincides with the identification curve, the curve-end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the curve coincides with the identification curve, it is assumed to
     * be smaller than any other object.
     * \param xcv the curve
     * \param ce the curve-end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the curve xcv.
     *   ARR_LEFT_BOUNDARY  - the curve approaches the identification curve
     *                        from the right at the curve left end.
     *   ARR_INTERIOR       - the curve does not approache the identification
     *                        curve.
     *   ARR_RIGHT_BOUNDARY - the curve approaches the identification curve
     *                        from the left at the curve right end.
     * \pre xcv does not coincide with the vertical identification curve.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      Comparison_result direction =
        geom_traits->compare_endpoints_xy_2_object()(xcv[0]);
      const X_monotone_subcurve_2& xs =
        (((direction == SMALLER) && (ce == ARR_MAX_END)) ||
         ((direction == LARGER) && (ce == ARR_MIN_END))) ?
        xcv[xcv.number_of_subcurves()-1] : xcv[0];
      return geom_traits->parameter_space_in_x_2_object()(xs, ce);
    }

    /*! Obtains the parameter space at a point along the x-axis.
     * \param p the point.
     * \return the parameter space at p.
     * \pre p does not lie on the vertical identification curve.
     */
    Arr_parameter_space operator()(const Point_2 p) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->parameter_space_in_x_2_object()(p);
    }
  };

  /*! Obtain a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(*this); }

  /*! A function object that obtains the parameter space of a geometric
   * entity along the y-axis
   */
  class Parameter_space_in_y_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Parameter_space_in_y_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Obtains the parameter space at the end of an curve along the y-axis .
     * Note that if the curve-end coincides with a pole, then unless the curve
     * coincides with the identification curve, the curve-end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the curve coincides with the identification curve, it is assumed to
     * be smaller than any other object.
     * \param xcv the curve
     * \param ce the curve-end indicator:
     *     ARR_MIN_END - the minimal end of xcv or
     *     ARR_MAX_END - the maximal end of xcv
     * \return the parameter space at the ce end of the curve xcv.
     *   ARR_BOTTOM_BOUNDARY  - the curve approaches the south pole at the
     *                          curve left end.
     *   ARR_INTERIOR         - the curve does not approache a contraction
     *                          point.
     *   ARR_TOP_BOUNDARY     - the curve approaches the north pole at the
     *                          curve right end.
     * There are no horizontal identification curves!
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      Comparison_result direction =
        geom_traits->compare_endpoints_xy_2_object()(xcv[0]);
      const X_monotone_subcurve_2& xs =
        (((direction == SMALLER) && (ce == ARR_MAX_END)) ||
         ((direction == LARGER) && (ce == ARR_MIN_END))) ?
        xcv[xcv.number_of_subcurves()-1] : xcv[0];
      return geom_traits->parameter_space_in_y_2_object()(xs, ce);
    }

    /*! Obtains the parameter space at a point along the y-axis.
     * \param p the point.
     * \return the parameter space at p.
     * \pre p does not lie on the horizontal identification curve.
     * There are no horizontal identification curves!
     */
    Arr_parameter_space operator()(const Point_2 p) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->parameter_space_in_y_2_object()(p);
    }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(*this); }

  /*! A functor that compares the x-coordinate of curve-ends and points on the
   * boundary of the parameter space.
   */
  class Compare_x_on_boundary_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Compare_x_on_boundary_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Compare the x-limit of a point with the x-coordinate of an
     * x-curve-end on the boundary.
     * \param point the point.
     * \param xcv the x-curve, the endpoint of which is compared.
     * \param ce the x-curve-end indicator:
     *            ARR_MIN_END - the minimal end of xcv or
     *            ARR_MAX_END - the maximal end of xcv.
     * \return the comparison result:
     *         SMALLER - x(p) < x(xcv, ce);
     *         EQUAL   - x(p) = x(xcv, ce);
     *         LARGER  - x(p) > x(xcv, ce).
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the x-curve xcv lies on the top boundary.
     * \pre xcv does not coincide with the vertical identification curve.
     */
    Comparison_result operator()(const Point_2& point,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    {
      const Subcurve_traits_2* geom_traits =
        m_poly_traits.subcurve_traits_2();
      Comparison_result direction =
        geom_traits->compare_endpoints_xy_2_object()(xcv[0]);
      const X_monotone_subcurve_2& xs =
        (((direction == SMALLER) && (ce == ARR_MAX_END)) ||
         ((direction == LARGER) && (ce == ARR_MIN_END))) ?
        xcv[0] : xcv[xcv.number_of_subcurves()-1];
      return geom_traits->compare_x_on_boundary_2_object()(point, xs, ce);
    }

    /*! Compare the x-coordinates of 2 curve-ends near the boundary of the
     * parameter space.
     * \param xcv1 the first curve.
     * \param ce1 the first curve-end indicator:
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second curve.
     * \param ce2 the second  curve-end indicator:
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) < x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce1 end of the curve xcv1 lies on a pole (implying ce1 is
     *      vertical).
     * \pre the ce2 end of the curve xcv2 lies on a pole (implying ce2 is
     *      vertical).
     * \pre xcv1 does not coincide with the vertical identification curve.
     * \pre xcv2 does not coincide with the vertical identification curve.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      Comparison_result direction1 =
        geom_traits->compare_endpoints_xy_2_object()(xcv1[0]);
      const X_monotone_subcurve_2& xs1 =
        (((direction1 == SMALLER) && (ce1 == ARR_MAX_END)) ||
         ((direction1 == LARGER) && (ce1 == ARR_MIN_END))) ?
        xcv1[0] : xcv1[xcv1.number_of_subcurves()-1];
      Comparison_result direction2 =
        geom_traits->compare_endpoints_xy_2_object()(xcv2[0]);
      const X_monotone_subcurve_2& xs2 =
        (((direction2 == SMALLER) && (ce2 == ARR_MAX_END)) ||
         ((direction2 == LARGER) && (ce2 == ARR_MIN_END))) ?
        xcv2[0] : xcv2[xcv2.number_of_subcurves()-1];
      return geom_traits->compare_x_on_boundary_2_object()(xs1, ce1, xs2, ce2);
    }
  };

  /*! Obtain a Compare_x_on_boundary_2 function object. */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
  { return Compare_x_on_boundary_2(*this); }

  class Compare_x_at_limit_2{
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    Compare_x_at_limit_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    unsigned int get_curve_index (const X_monotone_curve_2& xcv,
                                  const Arr_curve_end ce) const
    {
      //waqar:: dont know why it is opposite in Parameter_space_in_x...
      // I think this is because of the way the subcurves are stored in the
      // curve_vector.
      // I am assuming that min end depends upon the direction and not the
      // x-value.
      // and also that min end subcurve is always placed at position 0 of the
      // vector.
      // Comfirm with Eric.
      return (ce == ARR_MIN_END) ? 0 : xcv.number_of_subcurves() - 1;
    }

    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Compare_x_at_limit_2 compare_x_at_limit =
        geom_traits->compare_x_at_limit_2_object();

      unsigned int index = this->get_curve_index(xcv, ce);
      return compare_x_at_limit(p, xcv[index], ce );
    }

    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1/* for xcv1 */,
                                 const X_monotone_curve_2 & xcv2,
                                 Arr_curve_end ce2/*! for xcv2 */) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Compare_x_at_limit_2 compare_x_at_limit =
        geom_traits->compare_x_at_limit_2_object();

      unsigned int index_1 = this->get_curve_index(xcv1, ce1);
      unsigned int index_2 = this->get_curve_index(xcv2, ce2);

      return compare_x_at_limit(xcv1[index_1], ce1, xcv2[index_2], ce2);
    }

    Comparison_result operator()(const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce1/* for xcv */,
                                 const X_monotone_subcurve_2& xseg,
                                 Arr_curve_end ce2/*! for xseg */) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Compare_x_at_limit_2
        compare_x_at_limit = geom_traits->compare_x_at_limit_2_object();

      unsigned int index = this->get_curve_index(xcv, ce1 );

      return compare_x_at_limit(xcv[index], ce1, xseg, ce2 );
    }
  };

  Compare_x_at_limit_2 compare_x_at_limit_2_object() const
  { return Compare_x_at_limit_2(*this); }

  class Compare_x_near_limit_2{
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    Compare_x_near_limit_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    unsigned int get_curve_index (const X_monotone_curve_2& xcv,
                                  const Arr_curve_end ce) const
    {
      //waqar:: dont know why it is opposite in Parameter_space_in_x...
      // I think this is because of the way the subcurves are stored in the
      // curve_vector.
      // I am assuming that min end depends upon the direction and not the
      // x-value.
      // and also that min end subcurve is always placed at position 0 of the
      // vector.
      // Comfirm with Eric.
      unsigned int index =
        (ce == ARR_MIN_END) ? 0 : xcv.number_of_subcurves() - 1;
      return index;
    }

    Comparison_result operator()(const X_monotone_curve_2 xcv1,
                                 const X_monotone_curve_2 xcv2,
                                 Arr_curve_end ce) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Compare_x_near_limit_2
        cmp_x_near_limit = geom_traits->compare_x_near_limit_2_object();

      unsigned int index_1 = this->get_curve_index(xcv1, ce);
      unsigned int index_2 = this->get_curve_index(xcv2, ce);

      return cmp_x_near_limit(xcv1[index_1], xcv2[index_2], ce);
    }
  };

  Compare_x_near_limit_2 compare_x_near_limit_2_object() const
  { return Compare_x_near_limit_2(*this); }

  /*! A functor that compares the y-coordinate of two given points
   * that lie on the vertical identification curve.
   */
  class Compare_y_on_boundary_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Compare_y_on_boundary_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Compare the y-coordinate of two given points that lie on the vertical
     * identification curve.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return SMALLER - p1 is lexicographically smaller than p2;
     *         EQUAL   - p1 and p2 coincides;
     *         LARGER  - p1 is lexicographically larger than p2;
     * \pre p1 lies on the vertical identification curve.
     * \pre p2 lies on the vertical identification curve.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->compare_y_on_boundary_2_object()(p1, p2);
    }
  };

  /*! Obtain a Compare_y_on_boundary_2 function object */
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const
  { return Compare_y_on_boundary_2(*this); }

  /*! A functor that compares the y-coordinates of curve-ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Compare_y_near_boundary_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Compare the y-coordinates of 2 curves at their ends near the boundary
     * of the parameter space.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param ce the curve-end indicator:
     *     ARR_MIN_END - the minimal end or
     *     ARR_MAX_END - the maximal end
     * \return the second comparison result.
     * \pre the ce ends of the curves xcv1 and xcv2 lie either on the left
     *      boundary or on the right boundary of the parameter space (implying
     *      that they cannot be vertical).
     * There is no horizontal identification curve!
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      Comparison_result direction1 =
        geom_traits->compare_endpoints_xy_2_object()(xcv1[0]);
      const X_monotone_subcurve_2& xs1 =
        (((direction1 == SMALLER) && (ce == ARR_MAX_END)) ||
         ((direction1 == LARGER) && (ce == ARR_MIN_END))) ?
        xcv1[0] : xcv1[xcv1.number_of_subcurves()-1];
      Comparison_result direction2 =
        geom_traits->compare_endpoints_xy_2_object()(xcv2[0]);
      const X_monotone_subcurve_2& xs2 =
        (((direction2 == SMALLER) && (ce == ARR_MAX_END)) ||
         ((direction2 == LARGER) && (ce == ARR_MIN_END))) ?
        xcv2[0] : xcv2[xcv2.number_of_subcurves()-1];
      return geom_traits->compare_y_near_boundary_2_object()(xs1, xs2, ce);
    }
  };

  /*! Obtain a Compare_y_near_boundary_2 function object */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const
  { return Compare_y_near_boundary_2(*this); }

  /*! A functor that indicates whether a geometric object lies on the
   * vertical identification curve.
   */
  class Is_on_y_identification_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Is_on_y_identification_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Determine whether a point lies in the vertical boundary.
     * \param p the point.
     * \return a Boolean indicating whether p lies in the vertical boundary.
     */
    bool operator()(const Point_2& p) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->is_on_y_identification_2_object()(p);
    }

    /*! Determine whether an x-monotone curve lies in the vertical boundary.
     * \param xcv the x-monotone curve.
     * \return a Boolean indicating whether xcv lies in the vertical boundary.
     */
    bool operator()(const X_monotone_curve_2& xcv) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename X_monotone_curve_2::Subcurve_const_iterator it;
      for (it = xcv.subcurves_begin(); it != xcv.subcurves_end(); ++it)
        if (! geom_traits->is_on_y_identification_2_object()(*it)) return false;
      return true;
    }
  };

  /*! Obtain a Is_on_y_identification_2 function object */
  Is_on_y_identification_2 is_on_y_identification_2_object() const
  { return Is_on_y_identification_2(*this); }

  /*! A functor that indicates whether a geometric object lies on the
   * horizontal identification curve.
   */
  class Is_on_x_identification_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Is_on_x_identification_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Determine whether a point lies in the vertical boundary.
     * \param p the point.
     * \return a Boolean indicating whether p lies in the vertical boundary.
     */
    bool operator()(const Point_2& p) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->is_on_x_identification_2_object()(p);
    }

    /*! Determine whether an x-monotone curve lies in the vertical boundary.
     * \param xcv the x-monotone curve.
     * \return a Boolean indicating whether xcv lies in the vertical boundary.
     */
    bool operator()(const X_monotone_curve_2& xcv) const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename X_monotone_curve_2::Subcurve_const_iterator it;
      for (it = xcv.subcurves_begin(); it != xcv.subcurves_end(); ++it)
        if (! geom_traits->is_on_x_identification_2_object()(*it)) return false;
      return true;
    }
  };

  /*! Obtain a Is_on_x_identification_2 function object */
  Is_on_x_identification_2 is_on_x_identification_2_object() const
  { return Is_on_x_identification_2(*this); }

  //@}

  /// \name Types and functors defined here not required by any concept.
  //@{

  /*! \class
   * A functor that obtains the number of points of a polycurve.
   */
  class Number_of_points_2 {
  public:
    size_type operator()(const X_monotone_curve_2& cv) const
    {
      size_type num_seg = cv.number_of_subcurves();
      return (num_seg == 0) ? 0 : num_seg + 1;
    }
  };

  Number_of_points_2 number_of_points_2_object() const
  { return Number_of_points_2(); }

  /* Functor to augment a polycurve by adding a subcurve at the back.
   * TODO: Test all the operator()'s. (Don't forget vertical cases!)
   */
  class Push_back_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Push_back_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! Append a subcurve to an existing x-monotone polycurve at the back.
     */
    void operator()(X_monotone_curve_2& xcv, const X_monotone_subcurve_2& seg)
      const
    { push_back_2_impl<void*>(xcv, seg, Are_all_sides_oblivious_tag()); }

  private:
    // Oblivious implementation
    template<typename>
    void push_back_2_impl(X_monotone_curve_2& xcv,
                          const X_monotone_subcurve_2& seg,
                          Arr_all_sides_oblivious_tag) const
    {
      CGAL_precondition_code
        (
         typedef typename X_monotone_curve_2::size_type size_type;
         size_type num_seg = xcv.number_of_subcurves();
         const Subcurve_traits_2* geom_traits =
           m_poly_traits.subcurve_traits_2();
         typename Subcurve_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
           geom_traits->compare_endpoints_xy_2_object();
         Comparison_result dir = cmp_seg_endpts(seg);
         typename Subcurve_traits_2::Construct_max_vertex_2 get_max_v =
           geom_traits->construct_max_vertex_2_object();
         typename Subcurve_traits_2::Construct_min_vertex_2 get_min_v =
           geom_traits->construct_min_vertex_2_object();
         typename Subcurve_traits_2::Equal_2 equal =
           geom_traits->equal_2_object();
         typename Subcurve_traits_2::Is_vertical_2 is_vertical =
           geom_traits->is_vertical_2_object();

         CGAL_precondition_msg((num_seg == 0) ||
                               ((is_vertical(xcv[0]) && is_vertical(seg)) ||
                                (!is_vertical(xcv[0]) && !is_vertical(seg))),
                               "xcv is vertical and seg is not or vice versa!");

         CGAL_precondition_msg((num_seg == 0) ||
                               (cmp_seg_endpts(xcv[0]) == dir),
                               "xcv and seg do not have the same orientation!");

           CGAL_precondition_msg((num_seg == 0) ||
                                 !equal(get_min_v(seg), get_max_v(seg)),
                                 "Seg degenerates to a point!");

           CGAL_precondition_msg((num_seg == 0) ||
                                 (((dir != SMALLER) ||
                                   equal(get_max_v(xcv[num_seg-1]),
                                         get_min_v(seg)))),
                                 "Seg does not connect to the right!");

           CGAL_precondition_msg((num_seg == 0) ||
                                 (((dir != LARGER) ||
                                   equal(get_min_v(xcv[num_seg-1]),
                                         get_max_v(seg)))),
                                 "Seg does not connect to the left!");
         ); // precondition code ends

      xcv.push_back(seg);
    }

    // Boundary implementation
    template<typename>
    void push_back_2_impl(X_monotone_curve_2& xcv,
                          const X_monotone_subcurve_2& seg,
                          Arr_not_all_sides_oblivious_tag) const
    {
      CGAL_precondition_code
        (
         typedef typename X_monotone_curve_2::size_type size_type;
         size_type num_seg = xcv.number_of_subcurves();
         const Subcurve_traits_2* geom_traits =
           m_poly_traits.subcurve_traits_2();
         typename Subcurve_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
           geom_traits->compare_endpoints_xy_2_object();
         Comparison_result dir = cmp_seg_endpts(seg);
         typename Subcurve_traits_2::Construct_max_vertex_2 get_max_v =
           geom_traits->construct_max_vertex_2_object();
         typename Subcurve_traits_2::Construct_min_vertex_2 get_min_v =
           geom_traits->construct_min_vertex_2_object();
         typename Subcurve_traits_2::Equal_2 equal =
           geom_traits->equal_2_object();
         typename Subcurve_traits_2::Is_vertical_2 is_vertical =
           geom_traits->is_vertical_2_object();
         typename Subcurve_traits_2::Parameter_space_in_x_2 ps_x =
           geom_traits->parameter_space_in_x_2_object();
         typename Subcurve_traits_2::Parameter_space_in_y_2 ps_y =
           geom_traits->parameter_space_in_y_2_object();

         CGAL_precondition_msg((num_seg == 0) ||
                               ((is_vertical(xcv[0]) && is_vertical(seg)) ||
                                (!is_vertical(xcv[0]) && !is_vertical(seg))),
                               "xcv is vertical and seg is not or vice versa!");

         CGAL_precondition_msg((num_seg == 0) ||
                               (cmp_seg_endpts(xcv[0]) == dir),
                               "xcv and seg do not have the same orientation!");

         const Arr_parameter_space min_x_seg = ps_x(seg, ARR_MIN_END);
         const Arr_parameter_space min_y_seg = ps_y(seg, ARR_MIN_END);
         const Arr_parameter_space max_x_seg = ps_x(seg, ARR_MAX_END);
         const Arr_parameter_space max_y_seg = ps_y(seg, ARR_MAX_END);

         CGAL_precondition_code(const Arr_parameter_space min_x_cv =
                                ((num_seg>0) ?
                                 ps_x(xcv[num_seg-1], ARR_MIN_END) :
                                 ARR_INTERIOR));
         CGAL_precondition_code(const Arr_parameter_space min_y_cv =
                                ((num_seg>0) ?
                                 ps_y(xcv[num_seg-1], ARR_MIN_END) :
                                 ARR_INTERIOR));
         CGAL_precondition_code(const Arr_parameter_space max_x_cv =
                                ((num_seg>0) ?
                                 ps_x(xcv[num_seg-1], ARR_MAX_END) :
                                 ARR_INTERIOR));
         CGAL_precondition_code(const Arr_parameter_space max_y_cv =
                                ((num_seg>0) ?
                                 ps_y(xcv[num_seg-1], ARR_MAX_END) :
                                 ARR_INTERIOR));

         // A subcurve should not be pushed if the polycurve is directed to
         // the right and reaches the boundary.
         CGAL_precondition_msg(((dir != SMALLER) ||
                                ((max_x_cv == ARR_INTERIOR) &&
                                 (max_y_cv == ARR_INTERIOR))),
                               "Polycurve reaches the boundary to the right."
                               "Can not push back any subcurve further." );

         // A subcurve should not be pushed if the polycurve is directed to
         // the left and reaches the boundary.
         CGAL_precondition_msg(((dir != LARGER) ||
                                ((min_x_cv == ARR_INTERIOR) &&
                                 (min_y_cv == ARR_INTERIOR))),
                               "Polycurve reaches the boundary to the left."
                               "Can not push back any subcurve further." );

         // Something like a line should not be pushed if there is already a
         // subcurve present in the polycurve.
         CGAL_precondition_msg((((min_x_seg == ARR_INTERIOR) &&
                                 (min_y_seg == ARR_INTERIOR)) ||
                                ((max_x_seg == ARR_INTERIOR) &&
                                 (max_y_seg == ARR_INTERIOR)) ||
                                (num_seg == 0) ),
                               "Subcurve reaching the boundary at both ends "
                               "can not be pushed if there is already one or "
                               "more subcurves present in the polycurve.");

         if ((min_x_seg == ARR_INTERIOR) && (min_y_seg == ARR_INTERIOR) &&
             (max_x_seg == ARR_INTERIOR) && (max_y_seg == ARR_INTERIOR))
         {
           CGAL_precondition_msg((num_seg == 0) ||
                                 !equal(get_min_v(seg), get_max_v(seg)),
                                 "Seg degenerates to a point!");
         }

         if ((min_x_seg == ARR_INTERIOR) && (min_y_seg == ARR_INTERIOR)) {
           CGAL_precondition_msg((num_seg == 0) ||
                                 (((dir != SMALLER) ||
                                   equal(get_max_v(xcv[num_seg-1]),
                                         get_min_v(seg)))),
                                 "Seg does not connect to the right!");
         }

         if ((max_x_seg == ARR_INTERIOR) && (max_y_seg == ARR_INTERIOR) ) {
           CGAL_precondition_msg((num_seg == 0) ||
                                 (((dir != LARGER) ||
                                   equal(get_min_v(xcv[num_seg-1]),
                                         get_max_v(seg)))),
                                 "Seg does not connect to the left!");
         }
         ); // precondition code ends

      xcv.push_back(seg);
    }
  };

  /*! Obtain a Push_back_2 functor object. */
  Push_back_2 push_back_2_object() const { return Push_back_2(*this); }

  /* Functor to augment a polycurve by adding a subcurve at the front.
   * TODO: Test all the operator()'s. (Don't forget vertical cases!)
   */
  class Push_front_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

  public:
    /*! Constructor. */
    Push_front_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /* Append a subcurve `seg` to an existing polycurve `xcv` at the front. */
    void operator()(X_monotone_curve_2& xcv, const X_monotone_subcurve_2& seg)
      const
    { push_front_2_impl<void*>(xcv, seg, Are_all_sides_oblivious_tag()); }

  private:
    // Oblivious implementation
    template<typename>
    void push_front_2_impl(X_monotone_curve_2& xcv,
                           const X_monotone_subcurve_2& seg,
                           Arr_all_sides_oblivious_tag)
      const
    {
      CGAL_precondition_code
        (
         typedef typename X_monotone_curve_2::size_type size_type;
         size_type num_seg = xcv.number_of_subcurves();
         const Subcurve_traits_2* geom_traits =
           m_poly_traits.subcurve_traits_2();
         typename Subcurve_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
           geom_traits->compare_endpoints_xy_2_object();
         Comparison_result dir = cmp_seg_endpts(seg);
         typename Subcurve_traits_2::Construct_max_vertex_2 get_max_v =
           geom_traits->construct_max_vertex_2_object();
         typename Subcurve_traits_2::Construct_min_vertex_2 get_min_v =
           geom_traits->construct_min_vertex_2_object();
         typename Subcurve_traits_2::Equal_2 equal =
           geom_traits->equal_2_object();
         typename Subcurve_traits_2::Is_vertical_2 is_vertical =
           geom_traits->is_vertical_2_object();

         CGAL_precondition_msg((num_seg == 0) ||
                               ((is_vertical(xcv[0]) && is_vertical(seg)) ||
                                (!is_vertical(xcv[0]) && !is_vertical(seg))),
                               "xcv is vertical and seg is not or vice versa!");

         CGAL_precondition_msg((num_seg == 0) ||
                               (cmp_seg_endpts(xcv[0]) == dir),
                               "xcv and seg do not have the same orientation!");


           CGAL_precondition_msg((num_seg == 0) ||
                                 !equal(get_min_v(seg), get_max_v(seg)),
                                 "Seg degenerates to a point!");
           CGAL_precondition_msg((num_seg == 0) ||
                                 (((dir != SMALLER) ||
                                   equal(get_min_v(xcv[0]), get_max_v(seg)))),
                                 "Seg does not connect to the left!");

           CGAL_precondition_msg((num_seg == 0) ||
                                 (((dir != LARGER) ||
                                   equal(get_max_v(xcv[0]), get_min_v(seg)))),
                                 "Seg does not connect to the right!");
         ); // precondition code ends
      xcv.push_front(seg);
    }

    // Boundary implementation
    template<typename>
    void push_front_2_impl(X_monotone_curve_2& xcv,
                           const X_monotone_subcurve_2& seg,
                           Arr_not_all_sides_oblivious_tag)
      const
    {
      CGAL_precondition_code
        (
         typedef typename X_monotone_curve_2::size_type size_type;
         size_type num_seg = xcv.number_of_subcurves();
         const Subcurve_traits_2* geom_traits =
           m_poly_traits.subcurve_traits_2();
         typename Subcurve_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
           geom_traits->compare_endpoints_xy_2_object();
         Comparison_result dir = cmp_seg_endpts(seg);
         typename Subcurve_traits_2::Construct_max_vertex_2 get_max_v =
           geom_traits->construct_max_vertex_2_object();
         typename Subcurve_traits_2::Construct_min_vertex_2 get_min_v =
           geom_traits->construct_min_vertex_2_object();
         typename Subcurve_traits_2::Equal_2 equal =
           geom_traits->equal_2_object();
         typename Subcurve_traits_2::Parameter_space_in_x_2 ps_x =
           geom_traits->parameter_space_in_x_2_object();
         typename Subcurve_traits_2::Parameter_space_in_y_2 ps_y =
           geom_traits->parameter_space_in_y_2_object();
         typename Subcurve_traits_2::Is_vertical_2 is_vertical =
           geom_traits->is_vertical_2_object();

         CGAL_precondition_msg((num_seg == 0) ||
                               ((is_vertical(xcv[0]) && is_vertical(seg)) ||
                                (!is_vertical(xcv[0]) && !is_vertical(seg))),
                               "xcv is vertical and seg is not or vice versa!");

         CGAL_precondition_msg((num_seg == 0) ||
                               (cmp_seg_endpts(xcv[0]) == dir),
                               "xcv and seg do not have the same orientation!");

         const Arr_parameter_space min_x_seg = ps_x(seg, ARR_MIN_END);
         const Arr_parameter_space min_y_seg = ps_y(seg, ARR_MIN_END);
         const Arr_parameter_space max_x_seg = ps_x(seg, ARR_MAX_END);
         const Arr_parameter_space max_y_seg = ps_y(seg, ARR_MAX_END);

         // Something like line should not be pushed if there is already a
         // subcurve present in the polycurve.
         CGAL_precondition_msg((((min_x_seg == ARR_INTERIOR) &&
                                 (min_y_seg == ARR_INTERIOR)) ||
                                ((max_x_seg == ARR_INTERIOR) &&
                                 (max_y_seg == ARR_INTERIOR)) ||
                                (num_seg == 0)),
                               "Subcurve reaching the boundary at both ends"
                               "can not be pushed if there is already one "
                               "or more subcurves present in the polycurve.");

         // Something like Ray should not be pushed front if there
         // is already one or more subcurves present in polycurve.
         CGAL_precondition_msg((((dir == SMALLER) &&
                                 (max_x_seg == ARR_INTERIOR) &&
                                 (max_y_seg == ARR_INTERIOR)) ||
                                ((dir == LARGER) &&
                                 (min_x_seg == ARR_INTERIOR) &&
                                 (min_y_seg == ARR_INTERIOR)) ||
                                (num_seg == 0)),
                               "Subcurve reaching the boundary at the "
                               "connecting end can not be pushed ");

         if ((min_x_seg == ARR_INTERIOR) && (min_y_seg == ARR_INTERIOR) &&
             (max_x_seg == ARR_INTERIOR) && (max_y_seg == ARR_INTERIOR))
         {
           CGAL_precondition_msg((num_seg == 0) ||
                                 !equal(get_min_v(seg), get_max_v(seg)),
                                 "Seg degenerates to a point!");
         }

         if ((max_x_seg == ARR_INTERIOR) && (max_y_seg == ARR_INTERIOR)) {
           CGAL_precondition_msg((num_seg == 0) ||
                                 (((dir != SMALLER) ||
                                   equal(get_min_v(xcv[0]), get_max_v(seg)))),
                                 "Seg does not connect to the left!");
         }

         if ((min_x_seg == ARR_INTERIOR) && (min_y_seg == ARR_INTERIOR)) {
           CGAL_precondition_msg((num_seg == 0) ||
                                 (((dir != LARGER) ||
                                   equal(get_max_v(xcv[0]), get_min_v(seg)))),
                                 "Seg does not connect to the right!");
         }
         ); // precondition code ends

      xcv.push_front(seg);
    }
  };

  /*! Obtain a Push_front_2 functor object. */
  Push_front_2 push_front_2_object() const { return Push_front_2(*this); }

  class Trim_2 {
  protected:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /* The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

    /*! \brief returns a trimmed version of the polycurve with src and tgt as
     * end points.
     */
  public:
    /*! Constructor. */
    Trim_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv,
                                  const Point_2& src,
                                  const Point_2& tgt)const
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      typename Subcurve_traits_2::Trim_2 trim = geom_traits->trim_2_object();

      //check whether src and tgt lies on the polycurve/polycurve.
      CGAL_precondition(m_poly_traits.compare_y_at_x_2_object()(src, xcv) ==
                        EQUAL);
      CGAL_precondition(m_poly_traits.compare_y_at_x_2_object()(tgt, xcv) ==
                        EQUAL);

      /* Check whether the source and the target conform with the
       * direction of the polycurve.
       * since the direction of the poly-line/curve should not be changed.
       * we will interchange the source and the target.
       */
      Point_2 source = src;
      Point_2 target = tgt;

      // If curve is oriented from right to left but points are left to right.
      if (m_poly_traits.compare_endpoints_xy_2_object()(xcv) == LARGER &&
          m_poly_traits.compare_x_2_object()(src, tgt) == SMALLER )
      {
        source = tgt;
        target = src;
      }

      /* If curve is oriented from left to right but points are from right
       * to left.
       */
      else if (m_poly_traits.compare_endpoints_xy_2_object()(xcv) == SMALLER &&
               m_poly_traits.compare_x_2_object()(src, tgt) == LARGER )
      {
        source = tgt;
        target = src;
      }

      // std::cout << "**************the new sourc: " << source
      //           << "the new target: " << target << std::endl;
      /*
       * Get the source and target subcurve numbers from the polycurve.
       * The trimmed polycurve will have trimmed end subcurves(containing
       * source and target) along with complete
       * subcurves in between them.
       */
      std::size_t source_id = m_poly_traits.locate(xcv, source);
      std::size_t target_id = m_poly_traits.locate(xcv, target);
      // std::cout << "source number: " << source_id << "  Target number : "
      //           << target_id << std::endl;
      // std::cout << "target subcurve: " << xcv[target_id] << std::endl;

      std::vector<X_monotone_subcurve_2> trimmed_subcurves;

      Comparison_result orientation =
        m_poly_traits.compare_endpoints_xy_2_object()(xcv);

      Point_2 source_max_vertex =
        geom_traits->construct_max_vertex_2_object()(xcv[source_id]);
      Point_2 source_min_vertex =
        geom_traits->construct_min_vertex_2_object()(xcv[source_id]);
      Point_2 target_min_vertex =
        geom_traits->construct_min_vertex_2_object()(xcv[target_id]);
      Point_2 target_max_vertex =
        geom_traits->construct_max_vertex_2_object()(xcv[target_id]);

      //push the trimmed version of the source subcurve.
      // if(sorientation == SMALLER && source != source_max_vertex)
      if ((orientation == SMALLER) &&
          ! geom_traits->equal_2_object()(source, source_max_vertex) )
      {
        if (source_id != target_id )
          trimmed_subcurves.push_back(trim(xcv[source_id],
                                           source, source_max_vertex));
        else trimmed_subcurves.push_back(trim(xcv[source_id], source, target));
      }
      //else if(orientation == LARGER && source != source_min_vertex)
      else if ((orientation == LARGER) &&
               ! geom_traits->equal_2_object()(source, source_min_vertex))
      {
        if (source_id != target_id )
          trimmed_subcurves.push_back(trim(xcv[source_id],
                                          source, source_min_vertex));
        else trimmed_subcurves.push_back(trim(xcv[source_id], source, target));
      }

      //push the middle subcurves as they are.
      for (size_t i = source_id+1; i<target_id; ++i)
        trimmed_subcurves.push_back(xcv[i] );

      //push the appropriately trimmed target subcurve.
      if (source_id != target_id) {
        //if(orientation == SMALLER && target != target_min_vertex)
        if ((orientation == SMALLER) &&
            ! geom_traits->equal_2_object()(target, target_min_vertex))
          trimmed_subcurves.push_back(trim(xcv[target_id],
                                          target_min_vertex, target));

        //else if (orientation == LARGER && target != target_max_vertex)
        else if ((orientation == LARGER) &&
                 ! geom_traits->equal_2_object()(target, target_max_vertex))
          trimmed_subcurves.push_back(trim(xcv[target_id],
                                          target_max_vertex, target));
      }

      return X_monotone_curve_2(trimmed_subcurves.begin(),
                                trimmed_subcurves.end());
    }
  };

  /*! Obtain a Trim_2 functor object. */
  Trim_2 trim_2_object() const { return Trim_2(*this); }

  ///@}

protected:
  /*
   * Roadmap: locate() should return an iterator to the located subcurve
   */

  /*! Obtain the index of the subcurve in the polycurve that contains the
   * point q in its x-range. The function performs a binary search, so if the
   * point q is in the x-range of the polycurve with n subcurves, the subcurve
   * containing it can be located in O(log n) operations.
   * \param cv The polycurve curve.
   * \param q The point.
   * \return An index i such that q is in the x-range of cv[i].
   *         If q is not in the x-range of cv, returns INVALID_INDEX.
   */
  template <typename Compare>
  std::size_t locate_gen(const X_monotone_curve_2& cv, Compare compare) const
  {
    // The direction of cv. SMALLER means left-to-right and
    // otherwise right-to-left
    Comparison_result direction =
      subcurve_traits_2()->compare_endpoints_xy_2_object()(cv[0]);
    std::size_t from, to;
    if (direction == SMALLER) {
      from = 0;
      to = cv.number_of_subcurves() - 1;
    }
    else {
      from = cv.number_of_subcurves() - 1;
      to = 0;
    }

    // Test if q is one of cv's end points
    Comparison_result res_from = compare(cv[from], ARR_MIN_END);
    if (res_from == EQUAL) return from;

    Comparison_result res_to = compare(cv[to], ARR_MAX_END);
    if (res_to == EQUAL) return to;

    if (res_to == res_from) return INVALID_INDEX;

    // Perform a binary search to locate the subcurve that contains q in its
    // range:
    while (((direction == SMALLER) && (to > from)) ||
           ((direction == LARGER) && (to < from)))
    {
      std::size_t mid = (from + to) / 2;
      if (((direction == SMALLER) && (mid > from)) ||
          ((direction == LARGER) && (mid < from)))
      {
        Comparison_result res_mid = compare(cv[mid], ARR_MIN_END);
        if (res_mid == EQUAL) {
          // Ensure that the returned subcurve contains the query point
          // on its right end (if possible)
          if ((direction == SMALLER) && (mid > 0)) --mid;
          else if ((direction == LARGER) &&
                   ((mid + 1) < cv.number_of_subcurves()))
            ++mid;
          return mid;
        }
        if (res_mid == res_from) from = mid;
        else to = (direction == SMALLER) ? mid - 1 : mid + 1;
      }
      else {
        CGAL_assertion(((direction == SMALLER) && (mid < to)) ||
                       ((direction == LARGER) && (mid > to)));
        Comparison_result res_mid = compare(cv[mid], ARR_MAX_END);
        if (res_mid == EQUAL) return mid;
        if (res_mid == res_to) to = mid;
        else from = (direction == SMALLER) ? mid + 1 : mid - 1;
      }
    }
    // In case (from == to), and we know that the polycurve contains the q:
    CGAL_assertion(from == to);
    return from;
  }

  // A utility class that compare a curve-end with a point.
  template <typename Comparer>
  class Compare_points {
  private:
    typedef Arr_polycurve_basic_traits_2<Subcurve_traits_2>
      Polycurve_basic_traits_2;

    /*! The polycurve traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

    const Point_2& m_point;

    Comparer m_compare;

  public:
    // Constructor
    Compare_points(const Polycurve_basic_traits_2& traits, Comparer compare,
                   const Point_2& p) :
      m_poly_traits(traits),
      m_point(p),
      m_compare(compare)
    {}

    // Compare the given curve-end with the stored point.
    Comparison_result operator()(const X_monotone_subcurve_2& xs,
                                 Arr_curve_end ce)
    {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      const Point_2& p = (ce == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs) :
        geom_traits->construct_min_vertex_2_object()(xs);
      return m_compare(p, m_point);
    }
  };

  // A utility class that compares two curve-ends.
  template <typename Comparer>
  class Compare_point_curve_end {
  private:
    const Point_2& m_point;

    Comparer m_compare;

  public:
    // Constructor
    Compare_point_curve_end(Comparer compare, const Point_2& p) :
      m_point(p),
      m_compare(compare)
    {}

    // Compare a given curve-end with the stored point.
    Comparison_result operator()(const X_monotone_subcurve_2& xs,
                                 Arr_curve_end ce)
    { return m_compare(xs, ce, m_point); }
  };

  // A utility class that compare two curve-ends.
  template <typename Comparer>
  class Compare_curve_ends {
  private:
    const X_monotone_subcurve_2& m_x_monotone_subcurve;

    Arr_curve_end m_curve_end;

    Comparer m_compare;

  public:
    // Constructor
    Compare_curve_ends(Comparer compare,
                       const X_monotone_subcurve_2& xs, Arr_curve_end ce) :
      m_x_monotone_subcurve(xs),
      m_curve_end(ce),
      m_compare(compare)
    {}

    // Compare the given curve-end with the stored curve end.
    Comparison_result operator()(const X_monotone_subcurve_2& xs,
                                 Arr_curve_end ce)
    { return m_compare(xs, ce, m_x_monotone_subcurve, m_curve_end); }
  };

  /*! Locate the index of a curve in a polycurve that contains an endpoint
   * of a curve.
   * This implementation is used in the case where at least one side of the
   * parameter space is not oblivious.
   * \param xcv (in) the given polycurve.
   * \param xs (in) the given curve.
   * \param cd (in) the curve-end indicator.
   */
  std::size_t locate_impl(const X_monotone_curve_2& xcv,
                          const X_monotone_subcurve_2& xs,
                          Arr_curve_end ce,
                          Arr_not_all_sides_oblivious_tag) const
  {
    const Subcurve_traits_2* geom_traits = subcurve_traits_2();
    if (geom_traits->is_vertical_2_object()(xcv[0])) {
      // Verify that q has the same x-coord as xcv (which is vertical)
      Compare_x_2 compare_x = compare_x_2_object();
      Comparison_result res = compare_x(xcv[0], ARR_MIN_END, xs, ce);
      if (res != EQUAL) return INVALID_INDEX;

      Compare_curve_ends<Compare_xy_2> compare(compare_xy_2_object(), xs, ce);
      return locate_gen(xcv, compare);
    }

    Compare_curve_ends<Compare_x_2> compare(compare_x_2_object(), xs, ce);
    return locate_gen(xcv, compare);
  }

  /*! Locate the index of a curve in a polycurve that contains an endpoint
   * of a curve.
   * This implementation is used in the case where all sides of the parameter
   * space is oblivious.
   * \param xcv (in) the given polycurve.
   * \param xs (in) the given curve.
   * \param cd (in) the curve-end indicator.
   */
  std::size_t locate_impl(const X_monotone_curve_2& xcv,
                          const X_monotone_subcurve_2& xs,
                          Arr_curve_end ce,
                          Arr_all_sides_oblivious_tag) const
  {
    const Subcurve_traits_2* geom_traits = subcurve_traits_2();
    const Point_2& p = (ce == ARR_MAX_END) ?
      geom_traits->construct_max_vertex_2_object()(xs) :
      geom_traits->construct_min_vertex_2_object()(xs);
    return locate(xcv, p);
  }

  /*! Locate the index of a curve in a polycurve that contains an endpoint
   * of a curve.
   * This implementation is used in the case where at least one side of the
   * parameter space is not oblivious.
   * \param xcv (in) the given polycurve.
   * \param p (in) the endpoint of a curve.
   */
  std::size_t locate_impl(const X_monotone_curve_2& xcv,
                          const Point_2& p,
                          Arr_not_all_sides_oblivious_tag) const
  {
    const Subcurve_traits_2* geom_traits = subcurve_traits_2();
    if (geom_traits->is_vertical_2_object()(xcv[0])) {
      // Verify that q has the same x-coord as xcv (which is vertical)
      Compare_x_2 compare_x = compare_x_2_object();
      Comparison_result res = compare_x(xcv[0], ARR_MIN_END, p);
      if (res != EQUAL) return INVALID_INDEX;

      Compare_point_curve_end<Compare_xy_2> compare(compare_xy_2_object(), p);
      return locate_gen(xcv, compare);
    }

    Compare_point_curve_end<Compare_x_2> compare(compare_x_2_object(), p);
    return locate_gen(xcv, compare);
  }

  /*! Locate the index of a curve in a polycurve that contains an endpoint
   * of a curve.
   * This implementation is used in the case where all sides of the parameter
   * space is oblivious.
   * \param xcv (in) the given polycurve.
   * \param p (in) the endpoint of a curve.
   */
  std::size_t locate_impl(const X_monotone_curve_2& xcv, const Point_2& p,
                          Arr_all_sides_oblivious_tag) const
  { return locate(xcv, p); }

  //
  std::size_t locate(const X_monotone_curve_2& xcv, const Point_2& q) const
  {
    const Subcurve_traits_2* geom_traits = subcurve_traits_2();
    if (geom_traits->is_vertical_2_object()(xcv[0])) {
      // Verify that q has the same x-coord as cv (which is vertical)
      typename Subcurve_traits_2::Construct_min_vertex_2 min_vertex =
        geom_traits->construct_min_vertex_2_object();
      typename Subcurve_traits_2::Compare_x_2 compare_x =
        geom_traits->compare_x_2_object();
      Comparison_result res = compare_x(min_vertex(xcv[0]), q);
      if (res != EQUAL) return INVALID_INDEX;

      Compare_points<Compare_xy_2> compare(geom_traits,
                                           compare_xy_2_object(), q);
      return locate_gen(xcv, compare);
    }

    Compare_points<Compare_x_2> compare(geom_traits, compare_x_2_object(), q);
    return locate_gen(xcv, compare);
  }

  /*! Find the index of the subcurve in the polycurve that is defined to the
   * left(or to the right) of the point q.
   * \param cv The polycurve curve.
   * \param q The point.
   * \param to_right(true) if we wish to locate a subcurve to the right of q,
   *               (false) if we wish to locate a subcurve to its right.
   * \return An index i such that subcurves[i] is defined to the left(or to the
   *         right) of q, or INVALID_INDEX if no such subcurve exists.
   */
  std::size_t locate_side(const X_monotone_curve_2& cv,
                          const Point_2& q, const bool& to_right) const
  {
    // First locate a subcurve subcurves[i] that contains q in its x-range.
    std::size_t i = locate(cv, q);
    if (i == INVALID_INDEX) return INVALID_INDEX;

    typename Subcurve_traits_2::Equal_2 equal =
      subcurve_traits_2()->equal_2_object();
    typename Subcurve_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
      subcurve_traits_2()->compare_endpoints_xy_2_object();
    typename Subcurve_traits_2::Compare_x_2 comp_x =
      subcurve_traits_2()->compare_x_2_object();
    typename Subcurve_traits_2::Is_vertical_2 is_vert =
      subcurve_traits_2()->is_vertical_2_object();
    typename Subcurve_traits_2::Construct_max_vertex_2 get_max_v =
      subcurve_traits_2()->construct_max_vertex_2_object();
    typename Subcurve_traits_2::Construct_min_vertex_2 get_min_v =
      subcurve_traits_2()->construct_min_vertex_2_object();

    Comparison_result direction = cmp_seg_endpts(cv[i]);

    if ((!is_vert(cv[0]) && (comp_x(get_min_v(cv[i]), q) == EQUAL)) ||
        (is_vert(cv[0]) && equal(get_min_v(cv[i]), q))){
      // q is the left endpoint of the i'th subcurve:
      if (to_right) return i;
      else {
        // to_left
        if (direction == SMALLER) {
          if (i == 0) return INVALID_INDEX;
          else return i - 1;
        }
        else {
          if (i == cv.number_of_subcurves()-1) return INVALID_INDEX;
          else return i+1;
        }
      }
    }

    if ((!is_vert(cv[0]) && (comp_x(get_max_v(cv[i]), q) == EQUAL)) ||
        (is_vert(cv[0]) && equal(get_max_v(cv[i]), q)))
    {
      // q is the right endpoint of the i'th subcurve:
      if (!to_right) return i;
      else {
        if (direction == SMALLER) {
          if (i == (cv.number_of_subcurves() - 1)) return INVALID_INDEX;
          else return i + 1;
        }
        else {
          if (i == 0) return INVALID_INDEX;
          else return i-1;
        }
      }
    }

    // In case q is in cv[i]'s interior:
    return i;
  }
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
