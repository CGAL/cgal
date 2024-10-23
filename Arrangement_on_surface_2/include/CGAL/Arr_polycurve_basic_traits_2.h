// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University(Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
#include <type_traits>
#include <tuple>

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
  using Subcurve_traits_2 = SubcurveTraits_2;

  /// \name Types and functors inherited from the subcurve geometry traits.
  //@{

  using Has_left_category = typename Subcurve_traits_2::Has_left_category;
  using Has_do_intersect_category =
    typename Subcurve_traits_2::Has_do_intersect_category;

  using Left_side_category = typename Subcurve_traits_2::Left_side_category;
  using Bottom_side_category = typename Subcurve_traits_2::Bottom_side_category;
  using Top_side_category = typename Subcurve_traits_2::Top_side_category;
  using Right_side_category = typename Subcurve_traits_2::Right_side_category;

  using All_sides_oblivious_category =
    typename Arr_all_sides_oblivious_category<Left_side_category,
                                              Bottom_side_category,
                                              Top_side_category,
                                              Right_side_category>::result;

  using Bottom_or_top_sides_category =
    typename Arr_two_sides_category<Bottom_side_category,
                                    Top_side_category>::result;

  using Point_2 = typename Subcurve_traits_2::Point_2;
  using X_monotone_subcurve_2 = typename Subcurve_traits_2::X_monotone_curve_2;
  using Multiplicity = typename Subcurve_traits_2::Multiplicity;

  //@}

  // Backward compatibility:
  using X_monotone_segment_2 = X_monotone_subcurve_2;

private:
  using Self = Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

  // Data members:
  const Subcurve_traits_2* m_subcurve_traits;  // the base segment-traits class.
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
   * \param seg_traits an already existing subcurve tarits, which is passed in;
   *        it will be used by the class.
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

  /*! An \f$x\f$-monotone polycurve represents a continuous piecewise-linear
   * curve which is either strongly \f$x\f$-monotone or vertical. Again,
   * the polycurve is without degenerated subcurves.
   */
  using X_monotone_curve_2 =
    internal::X_monotone_polycurve_2<X_monotone_subcurve_2, Point_2>;
  using Size = typename X_monotone_curve_2::Size;
  using size_type = typename X_monotone_curve_2::size_type;

  //! Compare the \f$x\f$-coordinates of two points.
  class Compare_x_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Compare_x_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Compare the \f$x\f$-coordinates of two directional points.
     * \param p1 the first directional point.
     * \param p2 the second directional point.
     * \return `SMALLER` - \f$x\f$(`p1`) < \f$x\f$(`p2`);
     *         `EQUAL`   - \f$x\f$(`p1`) = \f$x\f$(`p2`);
     *         `LARGER`  - \f$x\f$(`p1`) > \f$x\f$(`p2`).
     * \pre p1 does not lie on the boundary.
     * \pre p2 does not lie on the boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { return m_poly_traits.subcurve_traits_2()->compare_x_2_object()(p1, p2); }

    /*! Compare two ends of \f$x\f$-monotone curves in \f$x\f$.
     * \param xs1 the first curve.
     * \param ce1 the curve-end indicator of the first \f$x\f$-monotone curve
     *        `xs1`:
     *            `ARR_MIN_END` - the minimal end of `xs1` or
     *            `ARR_MAX_END` - the maximal end of `xs1`.
     * \param p2 the second curve end.
     */
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1, const Point_2& p2)
    { return operator()(xs1, ce1, p2, All_sides_oblivious_category()); }

    /*! Compare two ends of \f$x\f$-monotone curves in x.
     * \param xs1 the first curve.
     * \param ce1 the curve-end indicator of the first \f$x\f$-monotone curve
     *        `xs1`:
     *            `ARR_MIN_END` - the minimal end of `xs1` or
     *            `ARR_MAX_END` - the maximal end of `xs1`.
     * \param xs2 the second curve.
     * \param ce2 the curve-end indicator of the second \f$x\f$-monoton curve
     *        `xs2`:
     *            `ARR_MIN_END` - the minimal end of `xs2` or
     *            `ARR_MAX_END` - the maximal end of `xs2`.
     */
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2)
    { return operator()(xs1, ce1, xs2, ce2, All_sides_oblivious_category()); }

  private:
    // Oblivious implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1, const Point_2& p2,
                                 Arr_all_sides_oblivious_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto p1 = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      return geom_traits->compare_x_2_object()(p1, p2);
    }

    // Boundary implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1, const Point_2& p2,
                                 Arr_not_all_sides_oblivious_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto ps_x = geom_traits->parameter_space_in_x_2_object();
      const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);

      if (ps_x1 != ARR_INTERIOR) {
        if (ps_x1 == ARR_LEFT_BOUNDARY) return SMALLER;
        if (ps_x1 == ARR_RIGHT_BOUNDARY) return LARGER;
      }

      auto ps_y = geom_traits->parameter_space_in_y_2_object();
      const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);
      if (ps_y1 == ARR_INTERIOR) {
        auto p1 = (ce1 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs1) :
          geom_traits->construct_min_vertex_2_object()(xs1);
        return geom_traits->compare_x_2_object()(p1, p2);
      }
      auto cmp_x_on_bnd = geom_traits->compare_x_on_boundary_2_object();
      return opposite(cmp_x_on_bnd(p2, xs1, ce1));
    }

    // Oblivious implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2,
                                 Arr_all_sides_oblivious_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto p1 = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      auto p2 = (ce2 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs2) :
        geom_traits->construct_min_vertex_2_object()(xs2);
      return geom_traits->compare_x_2_object()(p1, p2);
    }

    // Boundary implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2,
                                 Arr_not_all_sides_oblivious_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto ps_x = geom_traits->parameter_space_in_x_2_object();
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

      auto ps_y = geom_traits->parameter_space_in_y_2_object();
      const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);
      const Arr_parameter_space ps_y2 = ps_y(xs2, ce2);
      if (ps_y1 == ARR_INTERIOR) {
        auto p1 = (ce1 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs1) :
          geom_traits->construct_min_vertex_2_object()(xs1);
        if (ps_y2 == ARR_INTERIOR) {
          auto p2 = (ce2 == ARR_MAX_END) ?
            geom_traits->construct_max_vertex_2_object()(xs2) :
            geom_traits->construct_min_vertex_2_object()(xs2);
          return geom_traits->compare_x_2_object()(p1, p2);
        }
        auto cmp_x_on_bnd = geom_traits->compare_x_on_boundary_2_object();
        return cmp_x_on_bnd(p1, xs2, ce2);
      }
      if (ps_y2 == ARR_INTERIOR) {
        auto p2 = (ce2 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs2) :
          geom_traits->construct_min_vertex_2_object()(xs2);
        auto cmp_x_on_bnd = geom_traits->compare_x_on_boundary_2_object();
        return opposite(cmp_x_on_bnd(p2, xs1, ce1));
      }
      return geom_traits->compare_x_on_boundary_2_object()(xs1, ce1, xs2, ce2);
    }
  };

  /*! Obtain a Compare_x_2 functor object. */
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(*this); }

  //! Compare two curve-ends or points lexigoraphically: by x, then by y.
  class Compare_xy_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Compare_xy_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Compare two directional points lexigoraphically: by \f$x\f$, then by
     * \f$y\f$.
     * \param p1 the first endpoint directional point.
     * \param p2 the second endpoint directional point.
     * \return
     *   `SMALLER` - \f$x\f$(p1) < \f$x\f$(p2);
     *   `SMALLER` - \f$x\f$(p1) = \f$x\f$(p2) and \f$y\f$(p1) < \f$y\f$(p2);
     *   `EQUAL`   - \f$x\f$(p1) = \f$x\f$(p2) and \f$y\f$(p1) = \f$y\f$(p2);
     *   `LARGER`  - \f$x\f$(p1) = \f$x\f$(p2) and \f$y\f$(p1) > \f$y\f$(p2);
     *   `LARGER`  - \f$x\f$(p1) > \f$x\f$(p2).
     * \pre `p1` does not lie on the boundary.
     * \pre `p2` does not lie on the boundary.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const
    { return m_poly_traits.subcurve_traits_2()->compare_xy_2_object()(p1, p2); }

    /*! Compare two ends of \f$x\f$-monotone curves lexicographically.
     * \param xs1 the first curve.
     * \param ce1 the curve-end indicator of the first \f$x\f$-monotone curve
     *        `xs1`:
     *            `ARR_MIN_END` - the lexicographically smallest end of `xs1` or
     *            `ARR_MAX_END` - the lexicographically largest end of `xs1`.
     * \param p2 the second curve end.
     */
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1, const Point_2& p2)
    { return operator()(xs1, ce1, p2, All_sides_oblivious_category()); }

    /*! Compare two ends of \f$x\f$-monotone curves lexicographically.
     * \param xs1 the first curve.
     * \param ce1 the curve-end indicator of the first \f$x\f$-monotone curve
     *        `xs1`:
     *            `ARR_MIN_END` - the minimal end of `xs1` or
     *            `ARR_MAX_END` - the maximal end of `xs1`.
     * \param xs2 the second curve.
     * \param ce2 the curve-end indicator of the second \f$x\f$-monoton curve
     *        `xs2`:
     *            `ARR_MIN_END` - the minimal end of `xs2` or
     *            `ARR_MAX_END` - the maximal end of `xs2`.
     */
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2)
    { return operator()(xs1, ce1, xs2, ce2, All_sides_oblivious_category()); }

  private:
    // Oblivious implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1, const Point_2& p2,
                                 Arr_all_sides_oblivious_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto p1 = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      return geom_traits->compare_xy_2_object()(p1, p2);
    }

    // Boundary implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1, const Point_2& p2,
                                 Arr_not_all_sides_oblivious_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto ps_x = geom_traits->parameter_space_in_x_2_object();
      auto ps_y = geom_traits->parameter_space_in_y_2_object();
      const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);
      const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);

      if (ps_x1 != ARR_INTERIOR) {
        if (ps_x1 == ARR_LEFT_BOUNDARY) return SMALLER;
        if (ps_x1 == ARR_RIGHT_BOUNDARY) return LARGER;
      }

      if (ps_y1 == ARR_INTERIOR) {
        auto p1 = (ce1 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs1) :
          geom_traits->construct_min_vertex_2_object()(xs1);
        return geom_traits->compare_xy_2_object()(p1, p2);
      }

      // EFEF: missing implementation for open boundary.
      auto cmp_x_on_bnd = geom_traits->compare_x_on_boundary_2_object();
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
                                 Arr_all_sides_oblivious_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto p1 = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      auto p2 = (ce2 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs2) :
        geom_traits->construct_min_vertex_2_object()(xs2);
      return geom_traits->compare_xy_2_object()(p1, p2);
    }

    // Boundary implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_curve_end ce2,
                                 Arr_not_all_sides_oblivious_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto ps_x = geom_traits->parameter_space_in_x_2_object();
      auto ps_y = geom_traits->parameter_space_in_y_2_object();
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
        auto p1 = (ce1 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs1) :
          geom_traits->construct_min_vertex_2_object()(xs1);
        // ps1 == ARR_INTERIOR

        if ((ps_x2 == ARR_INTERIOR) && (ps_y2 == ARR_INTERIOR)) {
          auto p2 = (ce2 == ARR_MAX_END) ?
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
        auto cmp_x_on_bnd = geom_traits->compare_x_on_boundary_2_object();
        Comparison_result res = cmp_x_on_bnd(p1, xs2, ce2);
        if (res != EQUAL) return res;
        if (ps_y2 == ARR_TOP_BOUNDARY) return SMALLER;
        CGAL_assertion(ps_y2 == ARR_BOTTOM_BOUNDARY);
        return LARGER;
      }

      // ps1 != ARR_INTERIOR
      if ((ps_x2 == ARR_INTERIOR) && (ps_y2 == ARR_INTERIOR)) {
        auto p2 = (ce2 == ARR_MAX_END) ?
          geom_traits->construct_max_vertex_2_object()(xs2) :
          geom_traits->construct_min_vertex_2_object()(xs2);

        // The cases ps_x1 == ARR_{LEFT,RIGHT}_BOUNDARY are handled above

        // ps_x1 == ARR_INTERIOR
        // ps_y1 != ARR_INTERIOR
        // ps2 == ARR_INTERIOR
        CGAL_assertion(ps_x1 == ARR_INTERIOR);
        auto cmp_x_on_bnd = geom_traits->compare_x_on_boundary_2_object();
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
      auto p1 = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      auto p2 = (ce2 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs2) :
        geom_traits->construct_min_vertex_2_object()(xs2);
      return geom_traits->compare_y_on_boundary_2_object()(p1, p2);
    }
  };

  /*! Obtain a Compare_xy_2 functor object. */
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(*this); }

  /*! A functor that obtain the lexicographically smallest endpoint of an
   * \f$x\f$-monotone curve.
   */
  class Construct_min_vertex_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /* Constructor. */
    Construct_min_vertex_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Obtain the left endpoint of the \f$x\f$-monotone polycurve. The return
     * type is the same as the return type of the corresponding operator in the
     * functor of the subtraits, which is either by value or by reference.
     * \param xcv the polycurve curve.
     * \return the lexicographically smallest endpoint.
     */
    using Subcurve_ctr = typename Subcurve_traits_2::Construct_min_vertex_2;
    decltype(std::declval<Subcurve_ctr>().
             operator()(std::declval<X_monotone_subcurve_2>()))
    operator()(const X_monotone_curve_2& xcv) const {
      CGAL_assertion(xcv.number_of_subcurves() > 0);

      const auto* geom_traits = m_poly_traits.subcurve_traits_2();

      if (geom_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER)
        return geom_traits->construct_min_vertex_2_object()(xcv[0]);
      else
        return geom_traits->
          construct_min_vertex_2_object()(xcv[xcv.number_of_subcurves()-1]);
    }
  };

  /*! Obtain a Construct_min_vertex_2 functor object. */
  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(*this); }

  /*! A functor that obtain the lexicographically largest endpoint of an
   * \f$x\f$-monotone curve.
   */
  class Construct_max_vertex_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Construct_max_vertex_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Obtain the right endpoint of the \f$x\f$-monotone polycurve. The return
     * type is the same as the return type of the corresponding operator in the
     * functor of the subtraits, which is either by value or by reference.
     * \param xcv the polycurve.
     * \return the lexicographically largest endpoint.
     */
    using Subcurve_ctr = typename Subcurve_traits_2::Construct_max_vertex_2;
    decltype(std::declval<Subcurve_ctr>().
             operator()(std::declval<X_monotone_subcurve_2>()))
    operator()(const X_monotone_curve_2& xcv) const {
      CGAL_assertion(xcv.number_of_subcurves() > 0);
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto ctr_max_vertex = geom_traits->construct_max_vertex_2_object();
      if (geom_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER)
        return ctr_max_vertex(xcv[xcv.number_of_subcurves()-1]);
      else return ctr_max_vertex(xcv[0]);
    }
  };

  /*! Obtain a Construct_max_vertex_2 functor object. */
  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(*this); }

  /*! A functor that checks whether an \f$x\f$-monotone curve is a vertical. */
  class Is_vertical_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Is_vertical_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Check whether the given \f$x\f$-monotone curve is a vertical segment.
     * \param cv the curve.
     * \return `true` if the curve is a vertical segment; `false` otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv) const {
      /* An \f$x\f$-monotone polycurve can represent a vertical segment only if
       * it comprises vertical segments. If the first subcurve is vertical,
       * all subcurves are vertical in an \f$x\f$-monotone polycurve
       */
      return m_poly_traits.subcurve_traits_2()->is_vertical_2_object()(cv[0]);
    }
  };

  /*! Obtain an Is_vertical_2 functor object. */
  Is_vertical_2 is_vertical_2_object() const { return Is_vertical_2(*this); }

  /*! A functor that compares the \f$y\f$-coordinates of a point and an
   * \f$x\f$-monotone curve at the point \f$x\f$-coordinate
   */
  class Compare_y_at_x_2 {
  private:
    // Oblivious implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_all_sides_oblivious_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto p = (ce1 == ARR_MAX_END) ?
        geom_traits->construct_max_vertex_2_object()(xs1) :
        geom_traits->construct_min_vertex_2_object()(xs1);
      return geom_traits->compare_y_at_x_2_object()(p, xs2);
    }

    // Boundary implementation
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2,
                                 Arr_not_all_sides_oblivious_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto ps_x = geom_traits->parameter_space_in_x_2_object();
      auto ps_y = geom_traits->parameter_space_in_y_2_object();
      auto min_vertex = geom_traits->construct_min_vertex_2_object();
      auto max_vertex = geom_traits->construct_max_vertex_2_object();

      const Arr_parameter_space ps_x1 = ps_x(xs1, ce1);
      const Arr_parameter_space ps_y1 = ps_y(xs1, ce1);

      CGAL_assertion(((ce1 == ARR_MAX_END) && (ps_x1 != ARR_LEFT_BOUNDARY)) ||
                     ((ce1 == ARR_MIN_END) && (ps_x1 != ARR_RIGHT_BOUNDARY)));

      if (ps_x1 == ARR_INTERIOR) {
        auto p = (ce1 == ARR_MAX_END) ? max_vertex(xs1) : min_vertex(xs1);
        if (ps_y1 == ARR_TOP_BOUNDARY) {
          auto equal = geom_traits->equal_2_object();
          if (equal(p, max_vertex(xs2))) return EQUAL;
          if (equal(p, min_vertex(xs2))) return EQUAL;
          return LARGER;
        }
        if (ps_y1 == ARR_BOTTOM_BOUNDARY) {
          auto equal = geom_traits->equal_2_object();
          if (equal(p, max_vertex(xs2))) return EQUAL;
          if (equal(p, min_vertex(xs2))) return EQUAL;
          return SMALLER;
        }
        // ps_y1 == ARR_INTERIOR
        return geom_traits->compare_y_at_x_2_object()(p, xs2);
      }
      // ps_x1 == ARR_RIGHT_BOUNDARY || ARR_LEFT_BOUNDARY
      auto p1 = (ce1 == ARR_MAX_END) ? max_vertex(xs1) : min_vertex(xs1);
      auto p2 = (ce1 == ARR_MAX_END) ? max_vertex(xs2) : min_vertex(xs2);
      return geom_traits->compare_y_on_boundary_2_object()(p1, p2);
    }

    // Compare vertical
    Comparison_result compare_vertical(const Point_2& p,
                                       const X_monotone_curve_2& xcv) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      const Comparison_result l2r_smaller = SMALLER;
      const Comparison_result l2r_larger = LARGER;
#else
      auto cmp_endpints_xy = m_poly_traits.compare_endpoints_xy_2_object();
      const bool l2r = cmp_endpints_xy(xcv[0]) == SMALLER;
      const Comparison_result l2r_smaller = l2r ? SMALLER : LARGER;
      const Comparison_result l2r_larger = l2r ? LARGER : SMALLER;
#endif
      Comparison_result rc = geom_traits->compare_y_at_x_2_object()(p, xcv[0]);
      if (rc == l2r_smaller) return l2r_smaller;
      std::size_t n = xcv.number_of_subcurves();
      rc = geom_traits->compare_y_at_x_2_object()(p, xcv[n-1]);
      return (rc == l2r_larger) ? l2r_larger : EQUAL;
    }

    // Compare non_vertical
    Comparison_result compare_non_vertical(const Point_2& p,
                                           const X_monotone_curve_2& xcv) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      // Get the index of the subcurve in xcv containing p.
      auto i = m_poly_traits.locate_impl(xcv, p, All_sides_oblivious_category());
      CGAL_precondition(i != INVALID_INDEX);

      // Compare the subcurve xcv[i] and p.
      return geom_traits->compare_y_at_x_2_object()(p, xcv[i]);
    }

    // Oblivious implementation
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv,
                                 Arr_all_sides_oblivious_tag) const {
      if (! m_poly_traits.is_vertical_2_object()(xcv))
        return compare_non_vertical(p, xcv);
      return compare_vertical(p, xcv);
    }

    // Boundary implementation
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv,
                                 Arr_not_all_sides_oblivious_tag) const {
      if (! m_poly_traits.is_vertical_2_object()(xcv)) {
        const auto* geom_traits = m_poly_traits.subcurve_traits_2();
        if (geom_traits->is_on_y_identification_2_object()(p)) {
          auto cmp_y = geom_traits->compare_y_on_boundary_2_object();
          return cmp_y(p, m_poly_traits.construct_min_vertex_2_object()(xcv));
        }
        return compare_non_vertical(p, xcv);
      }
      return compare_vertical(p, xcv);
    }

  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Compare_y_at_x_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Obtain the location of the given point with respect to the input curve.
     * \param p the point.
     * \param xcv the polycurve curve.
     * \pre `p` is in the \f$x\f$-range of `xcv`.
     * \return
     *   `SMALLER` if \f$y\f$(p) < cv(x(p)), i.e. the point is below the curve;
     *   `LARGER` if \f$y\f$(p) > cv(x(p)), i.e. the point is above the curve;
     *   `EQUAL` if `p` lies on the curve.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv) const
    { return operator()(p, xcv, All_sides_oblivious_category()); }

    /*! Obtain the location of the given curve_end with respect to the input
     * curve.
     * \param xcv The polycurve curve.
     * \param ce the curve-end indicator of the \f$x\f$-monotone subcurve xl:
     *            `ARR_MIN_END` - the minimal end of xl or
     *            `ARR_MAX_END` - the maximal end of xl.
     * \param xcv The polycurve curve.
     * \pre the curve-end is in the \f$x\f$-range of `xcv`.
     * \return `SMALLER` if if \f$y\f$(xs, ce) < cv(x(xs, ce)), i.e.  the
     *          curve-end is below the curve xcv;
     *         `LARGER` if \f$y\f$(xs, ce) > cv(x(xs, ce)), i.e. the curve-end
     *          is above the curve `xcv`;
     *         `EQUAL` if the curve-end lies on the curve `xcv`.
     */
    Comparison_result operator()(const X_monotone_subcurve_2& xs1,
                                 Arr_curve_end ce1,
                                 const X_monotone_subcurve_2& xs2) const
    { return operator()(xs1, ce1, xs2, All_sides_oblivious_category()); }
  };

  /*! Obtain a Compare_y_at_x_2 functor object. */
  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(*this); }

  /*! A functor that compares the \f$y\f$-coordinates of two \f$x\f$-monotone
   * curves immediately to the left of their intersection point.
   */
  class Compare_y_at_x_left_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Compare_y_at_x_left_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Compare the y value of two \f$x\f$-monotone curves immediately to the
     * left of their intersection point.
     * \param cv1 the first polycurve curve.
     * \param cv2 the second polycurve curve.
     * \param p the intersection point.
     * \pre the point `p` lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return the relative position of `cv1` with respect to `cv2` immdiately
     *         to the left of `p`: `SMALLER`, `LARGER`, or `EQUAL`.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& p) const {
      // Get the indices of the subcurves in cv1 and cv2 containing p and
      // defined to its left.
      std::size_t i1 = m_poly_traits.locate_side(cv1, p, false);
      std::size_t i2 = m_poly_traits.locate_side(cv2, p, false);

      CGAL_precondition(i1 != INVALID_INDEX);
      CGAL_precondition(i2 != INVALID_INDEX);

      // Compare cv1[i1] and cv2[i2] at p's left.
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->compare_y_at_x_left_2_object()(cv1[i1], cv2[i2], p);
    }
  };

  /*! Obtain a Compare_y_at_x_left_2 functor object. */
  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(*this); }

  /*! A functor that compares the \f$y\f$-coordinates of two \f$x\f$-monotone
   * curves immediately to the right of their intersection point.
   */
  class Compare_y_at_x_right_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Compare_y_at_x_right_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Compare the \f$y\f$-value of two \f$x\f$-monotone curves immediately to
     * the right of their intersection point.
     * \param cv1 the first curve.
     * \param cv2 the second curve.
     * \param p the intersection point.
     * \pre the point `p` lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its right.
     * \return the relative position of `cv1` with respect to `cv2` immdiately
     *         to the right of `p`: `SMALLER`, `LARGER`, or `EQUAL`.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2,
                                 const Point_2& p) const {
      // Get the indices of the subcurves in cv1 and cv2 containing p and
      // defined to its right.
      std::size_t i1 = m_poly_traits.locate_side(cv1, p, true);
      std::size_t i2 = m_poly_traits.locate_side(cv2, p, true);

      CGAL_precondition(i1 != INVALID_INDEX);
      CGAL_precondition(i2 != INVALID_INDEX);

      // Compare cv1[i1] and cv2[i2] at p's right.
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->compare_y_at_x_right_2_object()(cv1[i1], cv2[i2], p);
    }
  };

  /*! Obtain a Compare_y_at_x_right_2 functor object.
   */
  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(*this); }

  /*! A functor that checks whether two points and two \f$x\f$-monotone
   * curves are identical.
   */
  class Equal_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Equal_2(const Polycurve_basic_traits_2& poly_tr) : m_poly_traits(poly_tr) {}

  public:
    /*! Check whether the two points are the same.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return `true` if the two point are the same; `false` otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const
    { return m_poly_traits.subcurve_traits_2()->equal_2_object()(p1, p2); }

    /*! Check whether the two \f$x\f$-monotone curves are the same (have the
     * same graph).
     * \param cv1 the first curve.
     * \param cv2 the second curve.
     * \return `true` if the two curves are the same; `false` otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const {
      // Check the pairwise equality of the contained subcurves.
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto equal = geom_traits->equal_2_object();
      auto cmp_x = geom_traits->compare_x_2_object();
      auto cmp_y_at_x = geom_traits->compare_y_at_x_2_object();
      auto max_vertex = geom_traits->construct_max_vertex_2_object();
      auto cmp_endpt = geom_traits->compare_endpoints_xy_2_object();
      auto is_vertical = m_poly_traits.is_vertical_2_object();
      auto xpoly_min_v = m_poly_traits.construct_min_vertex_2_object();
      auto xpoly_max_v = m_poly_traits.construct_max_vertex_2_object();

      // The first and last points of the subcurves should be equal.
      bool res = equal(xpoly_min_v(cv1), xpoly_min_v(cv2));
      if (! res) return false;
      res = equal(xpoly_max_v(cv1), xpoly_max_v(cv2));
      if (! res) return false;

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
      Comparison_result is_cv1_left_to_right = cmp_endpt(cv1[0]);
      Comparison_result is_cv2_left_to_right = cmp_endpt(cv2[0]);

      while ((i < (n1-1)) || (j < (n2-1))) {
        std::size_t cv1_seg_ind = (SMALLER == is_cv1_left_to_right) ?
          i :  n1 - 1 - i;
        auto point1 = max_vertex(cv1[cv1_seg_ind]);
        std::size_t cv2_seg_ind = (SMALLER == is_cv2_left_to_right) ?
          j :  n2 - 1 - j;
        auto point2 = max_vertex(cv2[cv2_seg_ind]);
        bool res = equal(point1, point2);
        // Easy case - the two points are equal
        if (res) {
          ++i;
          ++j;
        }
        else {
          Comparison_result res_x = cmp_x(point1, point2);
          // Check if the different point is a collinear point situated on
          // the line between its two neighbors.
          if (SMALLER == res_x) {
            Comparison_result res_y_at_x = cmp_y_at_x(point1, cv2[cv2_seg_ind]);
            if (EQUAL == res_y_at_x) ++i;
            else return false;
          }
          else if (LARGER == res_x) {
            Comparison_result res_y_at_x = cmp_y_at_x(point2, cv1[cv1_seg_ind]);
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

  /*! A functor that lexicographically compares the endpoints of a curve. */
  class Compare_endpoints_xy_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! The traits (in case it has state). */
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Compare_endpoints_xy_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Compare the endpoints of an \(x\)-monotone curve lexicographically.
     * (assuming the curve has a designated source and target points).
     * \param cv the curve.
     * \return `SMALLER` if `cv` is oriented left-to-right;
     *         `LARGER` if `cv` is oriented right-to-left.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_endpt = geom_traits->compare_endpoints_xy_2_object();
      return (cmp_endpt(xcv[0]) == SMALLER) ? (SMALLER) : (LARGER);
    }
  };
  //@}

  /// \name Types and functors defined here, required by the
  // ArrangementDirectionalXMonotoneTraits_2 concept.
  //@{

  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(*this); }

  /*! A functor that construct an \f$x\f$-monotone curve with the same endpoints
   * of a given curve, but directed in the opposite direction.
   */
  class Construct_opposite_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor */
    Construct_opposite_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Construct the reversed \f$x\f$-monotone polycurve of the input.
     * Note that the functor constructs the opposites of _all_ subcurves
     * constituting `xcv`.
     * \param xcv the \f$x\f$-monotone polycurve to be reveres
     * \pre xcv contains at least one subcurve
     * \return an \f$x\f$-monotone polycurve with the same graph as the input
     *         `xcv` only with a reverse orientation.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto const_op = geom_traits->construct_opposite_2_object();
      std::vector<X_monotone_subcurve_2> rev_segs(xcv.number_of_subcurves());;
      auto tit = rev_segs.begin();
      for (auto sit = xcv.subcurves_begin(); sit != xcv.subcurves_end(); ++sit)
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

  // The following block defines the nested types Approximate_number_type and
  // Approximate_2 and the member function approximate_2_object() based on the
  // corresponding types and function definitions of the subcurve traits. If
  // the subcurve traits does not provide these definitions, they are defined
  // as dummies. Essentially, the polycurve traits becomes a model of the
  // ArrangementLandmarkTraits concept only if the subcurve traits is a model
  // of this concept.
  //
  // The following implementation is inspired by
  // https://stackoverflow.com/a/11816999/1915421

  template <typename... Ts> using void_t = void;

  template <typename T, typename = void>
  struct has_approximate_2 {
    // Generic implementation
    using Approximate_number_type = void;
    using Approximate_point_2 = void;
    using Approximate_2 = void;
  };

  template <typename T>
  struct has_approximate_2<T, void_t<typename T::Approximate_2>> {
    // Specialization for types holding a nested type T::Approximate_2
    using Approximate_number_type = typename T::Approximate_number_type;
    using Approximate_2 = typename T::Approximate_2;
    using Approximate_point_2 = typename T::Approximate_point_2;
  };

  using Approximate_number_type =
    typename has_approximate_2<Subcurve_traits_2>::Approximate_number_type;
  using Approximate_2 =
    typename has_approximate_2<Subcurve_traits_2>::Approximate_2;
  using Approximate_point_2 =
    typename has_approximate_2<Subcurve_traits_2>::Approximate_point_2;

  /*! Obtain an Approximate_2 functor object. */
  Approximate_2 approximate_2_object_impl(std::false_type) const
  { return subcurve_traits_2()->approximate_2_object(); }

  Approximate_2 approximate_2_object_impl(std::true_type) const { }

  Approximate_2 approximate_2_object() const {
    using Is_void = typename std::is_same<void, Approximate_2>::type;
    return approximate_2_object_impl(Is_void());
  }

  //! A functor that constructs a point.
  class Construct_point_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    /*! Constructor. */
    Construct_point_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

  public:
    /*! Construct a point.
     * Apply perfect forwarding.
     */
    template <typename ... Args>
    Point_2 operator()(Args ... args) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto ctr_point = geom_traits->construct_point_2_object();
      return ctr_point(std::forward<Args>(args)...);
    }
  };

  /*! Obtain a Construct_x_monotone_curve_2 functor object. */
  Construct_point_2 construct_point_2_object() const
  { return Construct_point_2(*this); }

  //! A functor that constructs an \f$x\f$-monotone curve.
  class Construct_x_monotone_curve_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Construct_x_monotone_curve_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Obtain an \f$x\f$-monotone polycurve that consists of one given subcurve.
     * \param seg input subcurve.
     * \pre seg is not degenerated.
     * \return an \f$x\f$-monotone polycurve with one subcurve.
     */
    X_monotone_curve_2 operator()(const X_monotone_subcurve_2& seg) const {
      CGAL_precondition_code
        (
         /* Test that the subcurve is not degenerated. We do this test
          * independently from the subcurve traits in use, as we do not
          * allow a polycurve with degenerated subcurves.
          */
         const auto* geom_traits = m_poly_traits.subcurve_traits_2();
         auto get_min_v = geom_traits->construct_min_vertex_2_object();
         auto get_max_v = geom_traits->construct_max_vertex_2_object();
         auto equal = geom_traits->equal_2_object();

         CGAL_precondition_msg(! equal(get_min_v(seg), get_max_v(seg)),
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

    /*! Construct an \f$x\f$-monotone polycurve, which is well-oriented, from a
     * range of elements.
     * \pre the elements in the range should form a continuous well-oriented
     * \f$x\f$-monotone polycurve.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 operator()(ForwardIterator begin,
                                  ForwardIterator end) const {
      using VT = typename std::iterator_traits<ForwardIterator>::value_type;
      using Is_point = typename std::is_same<VT,Point_2>::type;

      // Dispatch the range to the appropriate implementation.
      return constructor_impl(begin, end, Is_point());
    }

    /*! Construct an \f$x\f$-monotone polycurve from a range of points.
     * The polycurve may be oriented left-to-right or right-to-left
     * depending on the lexicographical order of the points in the input.
     * \pre range contains at least two points.
     * \pre no two consecutive points are the same.
     * \pre the points form an continuous well-oriented \f$x\f$-monotone
     *      polycurve.
     * \post by the construction the returned polycurve is well-oriented.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 constructor_impl(ForwardIterator /* begin */,
                                        ForwardIterator /* end */,
                                        std::true_type) const
    { CGAL_error_msg("Cannot construct a polycurve from a range of points!"); }

    /*! Obtain an \f$x\f$-monotone polycurve from a range of subcurves.
     * \param begin An iterator pointing to the first subcurve in the range.
     * \param end An iterator pointing to the past-the-end subcurve
     * in the range.
     * \pre the range contains at least one subcurve.
     * \pre subcurves correspond to a well-oriented polycurve. That
     *      is, the target of the i-th subcurve is an source of the
     *      (i+1)th subcurve.
     * \pre the sequence of subcurves in the range forms a weak \f$x\f$-monotone
     *      polycurve.
     * \pre the container should support bidirectional iteration.
     * \return a continuous, well-oriented \f$x\f$-monotone polycurve ,which
     *         is directed either left-to-right or right-to-left
     *         depending on the subcurves in the input.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 constructor_impl(ForwardIterator begin,
                                        ForwardIterator end,
                                        std::false_type) const {
      CGAL_precondition_msg
        (
         begin != end,
         "Input range of subcurves has to contain at least one subcurve"
         );

      CGAL_precondition_code
        (
         const auto* geom_traits = m_poly_traits.subcurve_traits_2();
         auto cmp_seg_endpts = geom_traits->compare_endpoints_xy_2_object();
         auto get_min_v = geom_traits->construct_min_vertex_2_object();
         auto get_max_v = geom_traits->construct_max_vertex_2_object();
         auto equal = geom_traits->equal_2_object();

         ForwardIterator curr = begin;
         ForwardIterator next = begin;
         ++next;
         );


      CGAL_precondition_msg
        (
         (next != end) || ! equal(get_max_v(*curr),get_min_v(*curr)),
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
   * entity along the \f$x\f$-axis
   */
  class Parameter_space_in_x_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Parameter_space_in_x_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Obtains the parameter space at the end of a curve along the
     * \f$x\f$-axis. Note that if the curve-end coincides with a pole, then
     * unless the curve coincides with the identification curve, the curve-end
     * is considered to be approaching the boundary, but not on the boundary.
     * If the curve coincides with the identification curve, it is assumed to
     * be smaller than any other object.
     * \param xcv the curve
     * \param ce the curve-end indicator:
     *   `ARR_MIN_END` - the minimal end of `xcv` or
     *   `ARR_MAX_END` - the maximal end of `xcv`
     * \return the parameter space at the ce end of the curve xcv.
     *   `ARR_LEFT_BOUNDARY`  - the curve approaches the identification curve
     *                          from the right at the curve left end.
     *   `ARR_INTERIOR`       - the curve does not approache the identification
     *                          curve.
     *   `ARR_RIGHT_BOUNDARY` - the curve approaches the identification curve
     *                          from the left at the curve right end.
     * \pre xcv does not coincide with the vertical identification curve.
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_endpt = geom_traits->compare_endpoints_xy_2_object();
      Comparison_result direction = cmp_endpt(xcv[0]);
      const X_monotone_subcurve_2& xs =
        (((direction == SMALLER) && (ce == ARR_MIN_END)) ||
         ((direction == LARGER) && (ce == ARR_MAX_END))) ?
        xcv[0] : xcv[xcv.number_of_subcurves()-1];
      return geom_traits->parameter_space_in_x_2_object()(xs, ce);
    }

    /*! Obtains the parameter space at a point along the \f$x\f$-axis.
     * \param p the point.
     * \return the parameter space at `p`.
     * \pre `p` does not lie on the vertical identification curve.
     */
    Arr_parameter_space operator()(const Point_2 p) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->parameter_space_in_x_2_object()(p);
    }
  };

  /*! Obtain a Parameter_space_in_x_2 function object */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const
  { return Parameter_space_in_x_2(*this); }

  /*! A function object that obtains the parameter space of a geometric
   * entity along the \f$y\f$-axis
   */
  class Parameter_space_in_y_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Parameter_space_in_y_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Obtains the parameter space at the end of an curve along the
     * \f$y\f$-axis. Note that if the curve-end coincides with a pole, then
     * unless the curve coincides with the identification curve, the curve-end
     * is considered to be approaching the boundary, but not on the boundary.
     * If the curve coincides with the identification curve, it is assumed to
     * be smaller than any other object.
     * \param xcv the curve
     * \param ce the curve-end indicator:
     *    `ARR_MIN_END` - the minimal end of `xcv` or
     *    `ARR_MAX_END` - the maximal end of `xcv`
     * \return the parameter space at the ce end of the curve xcv.
     *   `ARR_BOTTOM_BOUNDARY`  - the curve approaches the south pole at the
     *                            curve left end.
     *   `ARR_INTERIOR`         - the curve does not approache a contraction
     *                            point.
     *   `ARR_TOP_BOUNDARY`     - the curve approaches the north pole at the
     *                            curve right end.
     * There are no horizontal identification curves!
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xcv,
                                   Arr_curve_end ce) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_endpt = geom_traits->compare_endpoints_xy_2_object();
      Comparison_result direction = cmp_endpt(xcv[0]);
      const X_monotone_subcurve_2& xs =
        (((direction == SMALLER) && (ce == ARR_MIN_END)) ||
         ((direction == LARGER) && (ce == ARR_MAX_END))) ?
        xcv[0] : xcv[xcv.number_of_subcurves()-1];
      return geom_traits->parameter_space_in_y_2_object()(xs, ce);
    }

    /*! Obtains the parameter space at a point along the \f$y\f$-axis.
     * \param p the point.
     * \return the parameter space at `p`.
     * \pre p does not lie on the horizontal identification curve.
     * There are no horizontal identification curves!
     */
    Arr_parameter_space operator()(const Point_2 p) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->parameter_space_in_y_2_object()(p);
    }
  };

  /*! Obtain a Parameter_space_in_y_2 function object */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const
  { return Parameter_space_in_y_2(*this); }

  /*! A functor that compares the \f$x\f$-coordinate of curve-ends and points on
   * the boundary of the parameter space.
   */
  class Compare_x_on_boundary_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Compare_x_on_boundary_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Compare the \f$x\f$-coordinates of a point with the \f$x\f$-coordinate
     * of an \f$x\f$-curve-end on the boundary.
     * \param point the point.
     * \param xcv the \f$x\f$-monotone curve, the endpoint of which is compared.
     * \param ce the \f$x\f$-monotone curve-end indicator:
     *            `ARR_MIN_END` - the minimal end of xcv or
     *            `ARR_MAX_END` - the maximal end of xcv.
     * \return the comparison result:
     *         `SMALLER` - \f$x\f$(`p`) < \f$x\f$(`xcv`, `ce`);
     *         `EQUAL`   - \f$x\f$(`p`) = \f$x\f$(`xcv`, `ce`);
     *         `LARGER`  - \f$x\f$(`p`) > \f$x\f$(`xcv`, `ce`).
     * \pre `p` lies in the interior of the parameter space.
     * \pre the `ce` end of `xcv` lies on the top boundary.
     * \pre `xcv` does not coincide with the vertical identification curve.
     */
    Comparison_result operator()(const Point_2& point,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce) const
    { return operator()(point, xcv, ce, Bottom_or_top_sides_category()); }

    /*! Compare the \f$x\f$-coordinates of 2 curve-ends on the boundary of the
     * parameter space.
     * \param xcv1 the first curve.
     * \param ce1 the first curve-end indicator:
     *            `ARR_MIN_END` - the minimal end of `xcv1` or
     *            `ARR_MAX_END` - the maximal end of `xcv1`.
     * \param xcv2 the second curve.
     * \param ce2 the second  curve-end indicator:
     *            `ARR_MIN_END` - the minimal end of `xcv2` or
     *            `ARR_MAX_END` - the maximal end of `xcv2`.
     * \return the second comparison result:
     *   `SMALLER` - \f$\f$x(`xcv1`, `ce1`) < \f$x\f$(`xcv2`, `ce2`);
     *   `EQUAL`   - \f$x\f$(`xcv1`, `ce1`) = \f$x\f$(`xcv2`, `ce2`);
     *   `LARGER`  - \f$x\f$(`xcv1`, `ce1`) > \f$x\f$(`xcv2`, `ce2`).
     * \pre the `ce1` end of `xcv1` lies on a pole (implying `xcv1` is
     *      vertical).
     * \pre the `ce2` end of `xcv2` lies on a pole (implying `xcv2` is
     *      vertical).
     * \pre `xcv1` does not coincide with the vertical identification curve.
     * \pre `xcv2` does not coincide with the vertical identification curve.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2) const
    { return operator()(xcv1, ce1, xcv2, ce2, Bottom_or_top_sides_category()); }

  private:
    /*! \brief compares the \f$x\f$-coordinates of a point with the
     * \f$x\f$-coordinate of an \f$x\f$-monotone curve-end on the boundary.
     */
    Comparison_result operator()(const Point_2& point,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce,
                                 Arr_boundary_cond_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_endpt = geom_traits->compare_endpoints_xy_2_object();
      Comparison_result direction = cmp_endpt(xcv[0]);
      const X_monotone_subcurve_2& xs =
        (((direction == SMALLER) && (ce == ARR_MIN_END)) ||
         ((direction == LARGER) && (ce == ARR_MAX_END))) ?
        xcv[0] : xcv[xcv.number_of_subcurves()-1];
      return geom_traits->compare_x_on_boundary_2_object()(point, xs, ce);
    }

    /*! \brief compares the \f$x\f$-coordinates of 2 curve-ends on the boundary
     * of the parameter space.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2,
                                 Arr_boundary_cond_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_endpt = geom_traits->compare_endpoints_xy_2_object();
      Comparison_result direction1 = cmp_endpt(xcv1[0]);
      const X_monotone_subcurve_2& xs1 =
        (((direction1 == SMALLER) && (ce1 == ARR_MIN_END)) ||
         ((direction1 == LARGER) && (ce1 == ARR_MAX_END))) ?
        xcv1[0] : xcv1[xcv1.number_of_subcurves()-1];
      Comparison_result direction2 = cmp_endpt(xcv2[0]);
      const X_monotone_subcurve_2& xs2 =
        (((direction2 == SMALLER) && (ce2 == ARR_MIN_END)) ||
         ((direction2 == LARGER) && (ce2 == ARR_MAX_END))) ?
        xcv2[0] : xcv2[xcv2.number_of_subcurves()-1];
      return geom_traits->compare_x_on_boundary_2_object()(xs1, ce1, xs2, ce2);
    }

    size_type get_curve_index(const X_monotone_curve_2& xcv,
                              const Arr_curve_end ce) const
    { return (ce == ARR_MIN_END) ? 0 : xcv.number_of_subcurves() - 1; }

    /*! Given a point \f$p\f$, an x-monotone curve \f$C(t) = (X(t),Y(t))\f$,
     * and an enumerator that specifies either the minimum end or the
     * maximum end of the curve, and thus maps to a parameter value
     * \f$d \in \{0,1\}\f$, compare x_p and limit{t => d} X(t).
     * If the parameter space is unbounded, a precondition ensures that \f$C\f$
     * has a vertical asymptote at its \f$d\f$-end; that is
     * limit{t => d} X(t) is finite.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce,
                                 Arr_has_open_side_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_x_on_boundary = geom_traits->compare_x_on_boundary_2_object();

      size_type index = this->get_curve_index(xcv, ce);
      return cmp_x_on_boundary(p, xcv[index], ce);
    }

    /*! Given two \f$x\f$-monotone curves \f$C_1(t) = (X_1(t),Y_1(t))\f$ and
     * \f$C2_(t) = (X_2(t),Y_2(t))\f$ and two enumerators that specify either
     * the minimum ends or the maximum ends of the curves, and thus map to
     * parameter values \f$d_1 \in \{0,1\}\f$ and \f$d_2 \in \{0,1\}\f$ for
     * \f$C_1\f$ and for \f$C_2\f$, respectively, compare
     * limit{t => d1} X1(t) and limit{t => d2} X2(t).
     * If the parameter space is unbounded, a precondition ensures that
     * \f$C_1\f$ and \f$C_2\f$ have vertical asymptotes at their respective
     * ends; that is, limit{t => d1} X1(t) and limit{t =? d2} X2(t) are finite.
    */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 Arr_curve_end ce1/* for xcv1 */,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce2/*! for xcv2 */,
                                 Arr_has_open_side_tag) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_x_on_boundary = geom_traits->compare_x_on_boundary_2_object();
      size_type index_1 = this->get_curve_index(xcv1, ce1);
      size_type index_2 = this->get_curve_index(xcv2, ce2);
      return cmp_x_on_boundary(xcv1[index_1], ce1, xcv2[index_2], ce2);
    }

    Comparison_result operator()(const X_monotone_curve_2& xcv,
                                 Arr_curve_end ce1/* for xcv */,
                                 const X_monotone_subcurve_2& xseg,
                                 Arr_curve_end ce2/*! for xseg */) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_x_on_boundary = geom_traits->compare_x_on_boundary_2_object();
      size_type index = this->get_curve_index(xcv, ce1);
      return cmp_x_on_boundary(xcv[index], ce1, xseg, ce2);
    }
  };

  /*! Obtain a Compare_x_on_boundary_2 function object. */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const
  { return Compare_x_on_boundary_2(*this); }

  /*! A functor that compares the \f$x\f$-coordinates of curve ends near the
   * boundary of the parameter space.
   */
  class Compare_x_near_boundary_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    Compare_x_near_boundary_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    size_type get_curve_index(const X_monotone_curve_2& xcv,
                              const Arr_curve_end ce) const
    { return (ce == ARR_MIN_END) ? 0 : xcv.number_of_subcurves() - 1; }

    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_x_near_boundary = geom_traits->compare_x_near_boundary_2_object();
      size_type index_1 = this->get_curve_index(xcv1, ce);
      size_type index_2 = this->get_curve_index(xcv2, ce);

      return cmp_x_near_boundary(xcv1[index_1], xcv2[index_2], ce);
    }
  };

  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const
  { return Compare_x_near_boundary_2(*this); }

  /*! A functor that compares the \f$y\f$-coordinate of two given points
   * that lie on the vertical identification curve.
   */
  class Compare_y_on_boundary_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Compare_y_on_boundary_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Compare the \f$y\f$-coordinate of two given points that lie on the
     * vertical identification curve.
     * \param p1 the first point.
     * \param p2 the second point.
     * \return `SMALLER` - `p1` is lexicographically smaller than `p2`;
     *         `EQUAL`   - `p1` and `p2` coincides;
     *         `LARGER`  - `p1` is lexicographically larger than `p2`;
     * \pre `p1` lies on the vertical identification curve.
     * \pre `p2` lies on the vertical identification curve.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->compare_y_on_boundary_2_object()(p1, p2);
    }
  };

  /*! Obtain a Compare_y_on_boundary_2 function object */
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const
  { return Compare_y_on_boundary_2(*this); }

  /*! A functor that compares the \f$y\f$-coordinates of curve-ends near the
   * boundary of the parameter space.
   */
  class Compare_y_near_boundary_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Compare_y_near_boundary_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Compare the \f$y\f$-coordinates of 2 curves at their ends near the
     * boundary of the parameter space.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param ce the curve-end indicator:
     *     `ARR_MIN_END` - the minimal end or
     *     `ARR_MAX_END` - the maximal end
     * \return the second comparison result.
     * \pre the `ce` ends of the curves `xcv1` and `xcv2` lie either on the left
     *      boundary or on the right boundary of the parameter space (implying
     *      that they cannot be vertical).
     * There is no horizontal identification curve!
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1,
                                 const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_endpt = geom_traits->compare_endpoints_xy_2_object();
      Comparison_result direction1 = cmp_endpt(xcv1[0]);
      const X_monotone_subcurve_2& xs1 =
        (((direction1 == SMALLER) && (ce == ARR_MIN_END)) ||
         ((direction1 == LARGER) && (ce == ARR_MAX_END))) ?
        xcv1[0] : xcv1[xcv1.number_of_subcurves()-1];
      Comparison_result direction2 = cmp_endpt(xcv2[0]);
      const X_monotone_subcurve_2& xs2 =
        (((direction2 == SMALLER) && (ce == ARR_MIN_END)) ||
         ((direction2 == LARGER) && (ce == ARR_MAX_END))) ?
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
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Is_on_y_identification_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Determine whether a point lies in the vertical boundary.
     * \param p the point.
     * \return a Boolean indicating whether `p` lies in the vertical boundary.
     */
    bool operator()(const Point_2& p) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->is_on_y_identification_2_object()(p);
    }

    /*! Determine whether an \f$x\f$-monotone curve lies in the vertical
     * boundary.
     * \param xcv the \f$x\f$-monotone curve.
     * \return a Boolean indicating whether `xcv` lies in the vertical boundary.
     */
    bool operator()(const X_monotone_curve_2& xcv) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      for (auto it = xcv.subcurves_begin(); it != xcv.subcurves_end(); ++it)
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
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Is_on_x_identification_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Determine whether a point lies in the vertical boundary.
     * \param p the point.
     * \return a Boolean indicating whether `p` lies in the vertical boundary.
     */
    bool operator()(const Point_2& p) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      return geom_traits->is_on_x_identification_2_object()(p);
    }

    /*! Determine whether an \f$x\f$-monotone curve lies in the vertical
     * boundary.
     * \param `xcv` the \f$x\f$-monotone curve.
     * \return a Boolean indicating whether `xcv` lies in the vertical boundary.
     */
    bool operator()(const X_monotone_curve_2& xcv) const {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      for (auto it = xcv.subcurves_begin(); it != xcv.subcurves_end(); ++it)
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
    size_type operator()(const X_monotone_curve_2& cv) const {
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
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Push_back_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /*! Append a subcurve to an existing \f$x\f$-monotone polycurve at the back.
     */
    void operator()(X_monotone_curve_2& xcv, const X_monotone_subcurve_2& seg)
      const
    { push_back_2_impl<void*>(xcv, seg, All_sides_oblivious_category()); }

  private:
    // Oblivious implementation
    template <typename>
    void push_back_2_impl(X_monotone_curve_2& xcv,
                          const X_monotone_subcurve_2& seg,
                          Arr_all_sides_oblivious_tag) const {
      CGAL_precondition_code
        (
         using size_type = typename X_monotone_curve_2::size_type;
         size_type num_seg = xcv.number_of_subcurves();
         const auto* geom_traits = m_poly_traits.subcurve_traits_2();
         auto cmp_endpts = geom_traits->compare_endpoints_xy_2_object();
         Comparison_result dir = cmp_endpts(seg);
         auto get_max_v = geom_traits->construct_max_vertex_2_object();
         auto get_min_v = geom_traits->construct_min_vertex_2_object();
         auto equal = geom_traits->equal_2_object();
         auto is_vertical = geom_traits->is_vertical_2_object();

         CGAL_precondition_msg((num_seg == 0) ||
                               ((is_vertical(xcv[0]) && is_vertical(seg)) ||
                                (! is_vertical(xcv[0]) && ! is_vertical(seg))),
                               "xcv is vertical and seg is not or vice versa!");

         CGAL_precondition_msg((num_seg == 0) ||
                               (cmp_endpts(xcv[0]) == dir),
                               "xcv and seg do not have the same orientation!");

           CGAL_precondition_msg((num_seg == 0) ||
                                 ! equal(get_min_v(seg), get_max_v(seg)),
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
    template <typename>
    void push_back_2_impl(X_monotone_curve_2& xcv,
                          const X_monotone_subcurve_2& seg,
                          Arr_not_all_sides_oblivious_tag) const {
      CGAL_precondition_code
        (
         using size_type = typename X_monotone_curve_2::size_type;
         size_type num_seg = xcv.number_of_subcurves();
         const auto* geom_traits = m_poly_traits.subcurve_traits_2();
         auto cmp_endpts = geom_traits->compare_endpoints_xy_2_object();
         Comparison_result dir = cmp_endpts(seg);
         auto get_max_v = geom_traits->construct_max_vertex_2_object();
         auto get_min_v = geom_traits->construct_min_vertex_2_object();
         auto equal = geom_traits->equal_2_object();
         auto is_vertical = geom_traits->is_vertical_2_object();
         auto ps_x = geom_traits->parameter_space_in_x_2_object();
         auto ps_y = geom_traits->parameter_space_in_y_2_object();

         CGAL_precondition_msg((num_seg == 0) ||
                               ((is_vertical(xcv[0]) && is_vertical(seg)) ||
                                (! is_vertical(xcv[0]) && ! is_vertical(seg))),
                               "xcv is vertical and seg is not or vice versa!");

         CGAL_precondition_msg((num_seg == 0) ||
                               (cmp_endpts(xcv[0]) == dir),
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
                               "Can not push back any subcurve further.");

         // A subcurve should not be pushed if the polycurve is directed to
         // the left and reaches the boundary.
         CGAL_precondition_msg(((dir != LARGER) ||
                                ((min_x_cv == ARR_INTERIOR) &&
                                 (min_y_cv == ARR_INTERIOR))),
                               "Polycurve reaches the boundary to the left."
                               "Can not push back any subcurve further.");

         // Something like a line should not be pushed if there is already a
         // subcurve present in the polycurve.
         CGAL_precondition_msg((((min_x_seg == ARR_INTERIOR) &&
                                 (min_y_seg == ARR_INTERIOR)) ||
                                ((max_x_seg == ARR_INTERIOR) &&
                                 (max_y_seg == ARR_INTERIOR)) ||
                                (num_seg == 0)),
                               "Subcurve reaching the boundary at both ends "
                               "can not be pushed if there is already one or "
                               "more subcurves present in the polycurve.");

         if ((min_x_seg == ARR_INTERIOR) && (min_y_seg == ARR_INTERIOR) &&
             (max_x_seg == ARR_INTERIOR) && (max_y_seg == ARR_INTERIOR))
         {
           CGAL_precondition_msg((num_seg == 0) ||
                                 ! equal(get_min_v(seg), get_max_v(seg)),
                                 "Seg degenerates to a point!");
         }

         if ((min_x_seg == ARR_INTERIOR) && (min_y_seg == ARR_INTERIOR)) {
           CGAL_precondition_msg((num_seg == 0) ||
                                 (((dir != SMALLER) ||
                                   equal(get_max_v(xcv[num_seg-1]),
                                         get_min_v(seg)))),
                                 "Seg does not connect to the right!");
         }

         if ((max_x_seg == ARR_INTERIOR) && (max_y_seg == ARR_INTERIOR)) {
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
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Push_front_2(const Polycurve_basic_traits_2& traits) :
      m_poly_traits(traits)
    {}

  public:
    /* Append a subcurve `seg` to an existing polycurve `xcv` at the front. */
    void operator()(X_monotone_curve_2& xcv, const X_monotone_subcurve_2& seg)
      const
    { push_front_2_impl<void*>(xcv, seg, All_sides_oblivious_category()); }

  private:
    // Oblivious implementation
    template <typename>
    void push_front_2_impl(X_monotone_curve_2& xcv,
                           const X_monotone_subcurve_2& seg,
                           Arr_all_sides_oblivious_tag) const {
      CGAL_precondition_code
        (
         using size_type = typename X_monotone_curve_2::size_type;
         size_type num_seg = xcv.number_of_subcurves();
         const auto* geom_traits = m_poly_traits.subcurve_traits_2();
         auto cmp_endpts = geom_traits->compare_endpoints_xy_2_object();
         Comparison_result dir = cmp_endpts(seg);
         auto get_max_v = geom_traits->construct_max_vertex_2_object();
         auto get_min_v = geom_traits->construct_min_vertex_2_object();
         auto equal = geom_traits->equal_2_object();
         auto is_vertical = geom_traits->is_vertical_2_object();

         CGAL_precondition_msg((num_seg == 0) ||
                               ((is_vertical(xcv[0]) && is_vertical(seg)) ||
                                (! is_vertical(xcv[0]) && !is_vertical(seg))),
                               "xcv is vertical and seg is not or vice versa!");

         CGAL_precondition_msg((num_seg == 0) ||
                               (cmp_endpts(xcv[0]) == dir),
                               "xcv and seg do not have the same orientation!");


           CGAL_precondition_msg((num_seg == 0) ||
                                 ! equal(get_min_v(seg), get_max_v(seg)),
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
    template <typename>
    void push_front_2_impl(X_monotone_curve_2& xcv,
                           const X_monotone_subcurve_2& seg,
                           Arr_not_all_sides_oblivious_tag) const {
      CGAL_precondition_code
        (
         using size_type = typename X_monotone_curve_2::size_type;
         size_type num_seg = xcv.number_of_subcurves();
         const auto* geom_traits = m_poly_traits.subcurve_traits_2();
         auto cmp_endpts = geom_traits->compare_endpoints_xy_2_object();
         Comparison_result dir = cmp_endpts(seg);
         auto get_max_v = geom_traits->construct_max_vertex_2_object();
         auto get_min_v = geom_traits->construct_min_vertex_2_object();
         auto equal = geom_traits->equal_2_object();
         auto ps_x = geom_traits->parameter_space_in_x_2_object();
         auto ps_y = geom_traits->parameter_space_in_y_2_object();
         auto is_vertical = geom_traits->is_vertical_2_object();

         CGAL_precondition_msg((num_seg == 0) ||
                               ((is_vertical(xcv[0]) && is_vertical(seg)) ||
                                (! is_vertical(xcv[0]) && ! is_vertical(seg))),
                               "xcv is vertical and seg is not or vice versa!");

         CGAL_precondition_msg((num_seg == 0) ||
                               (cmp_endpts(xcv[0]) == dir),
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
                                 ! equal(get_min_v(seg), get_max_v(seg)),
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

  //! A functor that trimps an \f$x\f$-monotone curve.
  class Trim_2 {
  protected:
    using Polycurve_basic_traits_2 =
      Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state).
    const Polycurve_basic_traits_2& m_poly_traits;

    friend class Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

    /*! Constructor. */
    Trim_2(const Polycurve_basic_traits_2& traits) : m_poly_traits(traits) {}

  public:
    /*! \brief returns a trimmed version of the polycurve with `source` and
     * `target` as end points.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv,
                                  const Point_2& source,
                                  const Point_2& target) const
    {
      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto min_vertex = geom_traits->construct_min_vertex_2_object();
      auto max_vertex = geom_traits->construct_max_vertex_2_object();
      auto trim = geom_traits->trim_2_object();

      //check whether src and tgt lies on the polycurve/polycurve.
      CGAL_precondition_code
        (auto cmp_y_at_x_2 = m_poly_traits.compare_y_at_x_2_object());
      CGAL_precondition(cmp_y_at_x_2(source, xcv) == EQUAL);
      CGAL_precondition(cmp_y_at_x_2(target, xcv) == EQUAL);

      /* Check whether the source and the target conform with the
       * direction of the polycurve.
       * since the direction of the poly-line/curve should not be changed.
       * we will interchange the source and the target.
       */

      /* If the curve is oriented from right to left but points are left to
       * right or if the curve is oriented from left to right but points are
       * from right to left, reverse.
       */
      auto [src, trg] =
        (((m_poly_traits.compare_endpoints_xy_2_object()(xcv) == LARGER) &&
          (m_poly_traits.compare_x_2_object()(source, target) == SMALLER)) ||
         ((m_poly_traits.compare_endpoints_xy_2_object()(xcv) == SMALLER) &&
          (m_poly_traits.compare_x_2_object()(source, target) == LARGER))) ?
        std::make_tuple(target, source) : std::make_tuple(source, target);

      // std::cout << "**************the new source: " << source
      //           << "the new target: " << target << std::endl;
      /* Get the source and target subcurve numbers from the polycurve.
       * The trimmed polycurve will have trimmed end subcurves(containing
       * source and target) along with complete
       * subcurves in between them.
       */
      std::size_t src_id = m_poly_traits.locate(xcv, src);
      std::size_t trg_id = m_poly_traits.locate(xcv, trg);
      // std::cout << "source number: " << source_id << "  Target number : "
      //           << target_id << std::endl;
      // std::cout << "target subcurve: " << xcv[target_id] << std::endl;

      std::vector<X_monotone_subcurve_2> trimmed_subcurves;

      Comparison_result orientation =
        m_poly_traits.compare_endpoints_xy_2_object()(xcv);

      auto src_max_vertex = max_vertex(xcv[src_id]);
      auto src_min_vertex = min_vertex(xcv[src_id]);
      auto trg_min_vertex = min_vertex(xcv[trg_id]);
      auto trg_max_vertex = max_vertex(xcv[trg_id]);

      // Push the trimmed version of the source subcurve.
      // if (sorientation == SMALLER && source != src_max_vertex)
      if ((orientation == SMALLER) &&
          ! geom_traits->equal_2_object()(src, src_max_vertex)) {
        if (src_id != trg_id)
          trimmed_subcurves.push_back(trim(xcv[src_id], src, src_max_vertex));
        else trimmed_subcurves.push_back(trim(xcv[src_id], src, trg));
      }
      // else if(orientation == LARGER && source != src_min_vertex)
      else if ((orientation == LARGER) &&
               ! geom_traits->equal_2_object()(src, src_min_vertex)) {
        if (src_id != trg_id)
          trimmed_subcurves.push_back(trim(xcv[src_id], src, src_min_vertex));
        else trimmed_subcurves.push_back(trim(xcv[src_id], src, trg));
      }

      // Push the middle subcurves as they are.
      for (size_t i = src_id+1; i < trg_id; ++i)
        trimmed_subcurves.push_back(xcv[i]);

      // Push the appropriately trimmed target subcurve.
      if (src_id != trg_id) {
        // if (orientation == SMALLER && target != trg_min_vertex)
        if ((orientation == SMALLER) &&
            ! geom_traits->equal_2_object()(trg, trg_min_vertex))
          trimmed_subcurves.push_back(trim(xcv[trg_id], trg_min_vertex, trg));

        // else if (orientation == LARGER && target != trg_max_vertex)
        else if ((orientation == LARGER) &&
                 ! geom_traits->equal_2_object()(trg, trg_max_vertex))
          trimmed_subcurves.push_back(trim(xcv[trg_id], trg_max_vertex, trg));
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

  /*! Obtain the index of the subcurve in the polycurve that contains a point
   * \f$q\f$ in its \f$x\f$-range. The function performs a binary search, so if
   * the point \f$q\f$ is in the \f$x\f$-range of the polycurve with \f$n\f$
   * subcurves, the subcurve containing it can be located in \cgalBigO{log n}
   * operations.
   * \param cv the polycurve curve.
   * \param q the point.
   * \return an index \f$i\f$ such that \f$q\f$ is in the \f$x\f$-range of
   *         `cv[i]`. If \f$q\f$ is not in the \f$x\f$-range of `cv`, returns
   *         `INVALID_INDEX`.
   */
  template <typename Compare>
  std::size_t locate_gen(const X_monotone_curve_2& cv, Compare compare) const {
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
           ((direction == LARGER) && (to < from))) {
      std::size_t mid = (from + to) / 2;
      if (((direction == SMALLER) && (mid > from)) ||
          ((direction == LARGER) && (mid < from))) {
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
    //! The polycurve traits (in case it has state).
    const Subcurve_traits_2& m_subcurve_traits;

    const Point_2& m_point;
    Comparer m_compare;

  public:
    // Constructor
    Compare_points(const Subcurve_traits_2& traits, Comparer compare,
                   const Point_2& p) :
      m_subcurve_traits(traits),
      m_point(p),
      m_compare(compare)
    {}

    // Compare the given curve-end with the stored point.
    Comparison_result operator()(const X_monotone_subcurve_2& xs,
                                 Arr_curve_end ce) {
      auto p = (ce == ARR_MAX_END) ?
        m_subcurve_traits.construct_max_vertex_2_object()(xs) :
        m_subcurve_traits.construct_min_vertex_2_object()(xs);
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
   * \param(in) xcv the given polycurve.
   * \param(in) xs the given curve.
   * \param(in) ce the curve-end indicator.
   */
  std::size_t locate_impl(const X_monotone_curve_2& xcv,
                          const X_monotone_subcurve_2& xs,
                          Arr_curve_end ce,
                          Arr_not_all_sides_oblivious_tag) const {
    const Subcurve_traits_2* geom_traits = subcurve_traits_2();
    if (geom_traits->is_vertical_2_object()(xcv[0])) {
      CGAL_precondition_code
        (
         // Verify that q has the same x-coord as xcv (which is vertical)
         Compare_x_2 cmp_x = compare_x_2_object();
         Comparison_result res = cmp_x(xcv[0], ARR_MIN_END, xs, ce);
         if (res != EQUAL) return INVALID_INDEX;
         );

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
   * \param(in) xcv the given polycurve.
   * \param(in) xs the given curve.
   * \param(in) cd the curve-end indicator.
   */
  std::size_t locate_impl(const X_monotone_curve_2& xcv,
                          const X_monotone_subcurve_2& xs,
                          Arr_curve_end ce,
                          Arr_all_sides_oblivious_tag) const {
    const Subcurve_traits_2* geom_traits = subcurve_traits_2();
    auto p = (ce == ARR_MAX_END) ?
      geom_traits->construct_max_vertex_2_object()(xs) :
      geom_traits->construct_min_vertex_2_object()(xs);
    return locate(xcv, p);
  }

  /*! Locate the index of a curve in a polycurve that contains a point.
   * This implementation is used in the case where at least one side of the
   * parameter space is not oblivious.
   * \param(in) xcv the given polycurve.
   * \param(in) p the query point.
   */
  std::size_t locate_impl(const X_monotone_curve_2& xcv,
                          const Point_2& p,
                          Arr_not_all_sides_oblivious_tag) const {
    const Subcurve_traits_2* geom_traits = subcurve_traits_2();
    if (geom_traits->is_vertical_2_object()(xcv[0])) {
      CGAL_precondition_code
        (
         // Verify that q has the same x-coord as xcv (which is vertical)
         auto cmp_x = compare_x_2_object();
         Comparison_result res = cmp_x(xcv[0], ARR_MIN_END, p);
         if (res != EQUAL) return INVALID_INDEX;
         );

      Compare_point_curve_end<Compare_xy_2> compare(compare_xy_2_object(), p);
      return locate_gen(xcv, compare);
    }

    Compare_point_curve_end<Compare_x_2> compare(compare_x_2_object(), p);
    return locate_gen(xcv, compare);
  }

  /*! Locate the index of a curve in a polycurve that contains a point.
   * This implementation is used in the case where all sides of the parameter
   * space are oblivious.
   * \param(in) xcv the given polycurve.
   * \param(in) p the query point.
   */
  std::size_t locate_impl(const X_monotone_curve_2& xcv, const Point_2& p,
                          Arr_all_sides_oblivious_tag) const
  { return locate(xcv, p); }

  //
  std::size_t locate(const X_monotone_curve_2& xcv, const Point_2& q) const {
    const Subcurve_traits_2* geom_traits = subcurve_traits_2();
    if (geom_traits->is_vertical_2_object()(xcv[0])) {
      CGAL_precondition_code
        (
         // Verify that q has the same x-coord as cv (which is vertical)
         auto min_vertex = geom_traits->construct_min_vertex_2_object();
         auto cmp_x = geom_traits->compare_x_2_object();
         Comparison_result res = cmp_x(min_vertex(xcv[0]), q);
         if (res != EQUAL) return INVALID_INDEX;
         );

      Compare_points<Compare_xy_2> compare(*geom_traits,
                                           compare_xy_2_object(), q);
      return locate_gen(xcv, compare);
    }

    Compare_points<Compare_x_2> compare(*geom_traits, compare_x_2_object(), q);
    return locate_gen(xcv, compare);
  }

  /*! Find the index of the subcurve in the polycurve that is defined to the
   * left (or to the right) of the point `q`.
   * \param cv the polycurve curve.
   * \param q the point.
   * \param to_right `true` if we wish to locate a subcurve to the right of q,
   *                 `false` if we wish to locate a subcurve to its right.
   * \return an index \f$i\f$ such that subcurves[i] is defined to the left (or
   *         to the right) of `q`, or `INVALID_INDEX` if no such subcurve exists.
   */
  std::size_t locate_side(const X_monotone_curve_2& cv,
                          const Point_2& q, const bool& to_right) const {
    // First locate a subcurve subcurves[i] that contains q in its x-range.
    std::size_t i = locate(cv, q);
    if (i == INVALID_INDEX) return INVALID_INDEX;

    auto equal = subcurve_traits_2()->equal_2_object();
    auto cmp_endpts = subcurve_traits_2()->compare_endpoints_xy_2_object();
    auto cmp_x = subcurve_traits_2()->compare_x_2_object();
    auto is_vert = subcurve_traits_2()->is_vertical_2_object();
    auto get_max_v = subcurve_traits_2()->construct_max_vertex_2_object();
    auto get_min_v = subcurve_traits_2()->construct_min_vertex_2_object();

    Comparison_result direction = cmp_endpts(cv[i]);

    if ((! is_vert(cv[0]) && (cmp_x(get_min_v(cv[i]), q) == EQUAL)) ||
        (is_vert(cv[0]) && equal(get_min_v(cv[i]), q))) {
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

    if ((! is_vert(cv[0]) && (cmp_x(get_max_v(cv[i]), q) == EQUAL)) ||
        (is_vert(cv[0]) && equal(get_max_v(cv[i]), q))) {
      // q is the right endpoint of the i'th subcurve:
      if (! to_right) return i;
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
