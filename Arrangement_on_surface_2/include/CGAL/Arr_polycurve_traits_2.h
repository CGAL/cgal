// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University(Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel <efif@post.tau.ac.il>
//            Ron Wein  <wein@post.tau.ac.il>
//            Dror Atariah <dror.atariah@fu-berlin.de>
//            Waqar Khan <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_ARR_POLYCURVE_TRAITS_2_H
#define CGAL_ARR_POLYCURVE_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The traits-class for the general piece-wise (polycurve) type of curves of the
 * arrangement package.
 */

#include <iterator>
#include <type_traits>
#include <variant>

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polycurve_basic_traits_2.h>
#include <CGAL/Arr_geometry_traits/Polycurve_2.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>

namespace CGAL {

template <typename SubcurveTraits_2 = Arr_segment_traits_2<> >
class Arr_polycurve_traits_2 :
    public Arr_polycurve_basic_traits_2<SubcurveTraits_2> {
public:
  using Subcurve_traits_2 = SubcurveTraits_2;

private:
  using Base = Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

public:
  /// \name Types inherited from the polycurve basic traits class.
  //@{
  using Has_left_category = typename Base::Has_left_category;
  using Has_do_intersect_category = typename Base::Has_do_intersect_category;

  using Left_side_category = typename Base::Left_side_category;
  using Bottom_side_category = typename Base::Bottom_side_category;
  using Top_side_category = typename Base::Top_side_category;
  using Right_side_category = typename Base::Right_side_category;

  using All_sides_oblivious_category =
    typename Base::All_sides_oblivious_category;

  using X_monotone_subcurve_2 = typename Base::X_monotone_subcurve_2;
  using Size = typename Base::Size;
  using size_type = typename Base::size_type;

  using Point_2 = typename Base::Point_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  using Compare_x_2 = typename Base::Compare_x_2;
  using Compare_xy_2 = typename Base::Compare_xy_2;
  using Construct_min_vertex_2 = typename Base::Construct_min_vertex_2;
  using Construct_max_vertex_2 = typename Base::Construct_max_vertex_2;
  using Is_vertical_2 = typename Base::Is_vertical_2;
  using Compare_y_at_x_2 = typename Base::Compare_y_at_x_2;
  using Compare_y_at_x_left_2 = typename Base::Compare_y_at_x_left_2;
  using Compare_y_at_x_right_2 = typename Base::Compare_y_at_x_right_2;
  using Equal_2 = typename Base::Equal_2;
  using Compare_endpoints_xy_2 = typename Base::Compare_endpoints_xy_2;
  using Construct_opposite_2 = typename Base::Construct_opposite_2;
  using Approximate_2 = typename Base::Approximate_2;
  using Construct_x_monotone_curve_2 =
    typename Base::Construct_x_monotone_curve_2;
  using Parameter_space_in_x_2 = typename Base::Parameter_space_in_x_2;
  using Parameter_space_in_y_2 = typename Base::Parameter_space_in_y_2;
  using Compare_x_on_boundary_2 = typename Base::Compare_x_on_boundary_2;
  using Compare_x_near_boundary_2 = typename Base::Compare_x_near_boundary_2;
  using Compare_y_on_boundary_2 = typename Base::Compare_y_on_boundary_2;
  using Compare_y_near_boundary_2 = typename Base::Compare_y_near_boundary_2;
  using Is_on_y_identification_2 = typename Base::Is_on_y_identification_2;
  using Is_on_x_identification_2 = typename Base::Is_on_x_identification_2;

  using Trim_2 = typename Base::Trim_2;

  //@}

  /// \name Types and functors inherited from the subcurve geometry traits.
  //@{

  using Has_merge_category = typename Subcurve_traits_2::Has_merge_category;
  using Multiplicity = typename Subcurve_traits_2::Multiplicity;
  using Subcurve_2 = typename Subcurve_traits_2::Curve_2;

  //@}

  // Backward compatibility:
  using Segment_2 = Subcurve_2;

private:
  using Self = Arr_polycurve_traits_2<Subcurve_traits_2>;

public:
  /*! constructs default */
  Arr_polycurve_traits_2() : Base() {}

  /*! constructs with given subcurve traits
   * \param seg_traits an already existing subcurve tarits which is passed will
   *        be used by the class.
   */
  Arr_polycurve_traits_2(const Subcurve_traits_2* geom_traits) :
    Base(geom_traits)
  {}

  /*! A polycurve represents a general continuous piecewise-linear
   * curve, without degenerated subcurves.
   */
  using Curve_2 = internal::Polycurve_2<Subcurve_2, Point_2>;

  /// \name Basic predicate functors(based on the subcurve traits).
  //@{

  /*! \class
   * A functor that obtains the number of points of a polycurve.
   */
  class Number_of_points_2 : public Base::Number_of_points_2 {
  public:
    size_type operator()(const Curve_2& cv) const {
      size_type num_seg = cv.number_of_subcurves();
      return (num_seg == 0) ? 0 : num_seg + 1;
    }
  };

  /*! obtains a `Number_of_points_2` functor object. */
  Number_of_points_2 number_of_points_2_object() const
  { return Number_of_points_2(); }

  ///@}

  /// \name Construction functors(based on the subcurve traits).
  //@{

#ifndef DOXYGEN_RUNNING
  class Push_back_2;
#endif

  //! A functor for subdividing curves into x-monotone curves.
  class Make_x_monotone_2 {
  protected:
    using Polycurve_traits_2 = Arr_polycurve_traits_2<Subcurve_traits_2>;

    //! The traits (in case it has state)
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! constructs. */
    Make_x_monotone_2(const Polycurve_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! subdivides a given curve into x-monotone sub-curves and insert them into
     * a given output iterator.
     *
     * \pre if `cv` is not empty then it must be continuous and well-oriented.
     * \param cv the curve.
     * \param oi an output iterator for the result. Its value type is a variant
     *           that wraps Point_2 or an X_monotone_curve_2 objects.
     * \return the past-the-end iterator.
     */
  private:
    template <typename OutputIterator>
    OutputIterator operator_impl(const Curve_2& cv, OutputIterator oi,
                                 Arr_all_sides_oblivious_tag) const {
      using Make_x_monotone_subresult =
        std::variant<Point_2, X_monotone_subcurve_2>;
      using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;

      // If the polycurve is empty, return.
      if (cv.number_of_subcurves() == 0) return oi;

      auto ctr_x_curve = m_poly_traits.construct_x_monotone_curve_2_object();

      auto make_seg_x_monotone =
        m_poly_traits.subcurve_traits_2()->make_x_monotone_2_object();

      auto cmp_seg_endpts =
        m_poly_traits.subcurve_traits_2()->compare_endpoints_xy_2_object();

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      auto ctr_seg_opposite =
        m_poly_traits.subcurve_traits_2()->construct_opposite_2_object();
#endif

      // Convert the input polycurve to a sequence of variant objects, such
      // that each object wraps an x-monotone subcurve.
      std::vector<Make_x_monotone_subresult> x_seg_objects;
      for (auto its = cv.subcurves_begin(); its != cv.subcurves_end(); ++its)
        make_seg_x_monotone(*its, std::back_inserter(x_seg_objects));
      auto it = x_seg_objects.begin();
      const auto* x_seg_p = std::get_if<X_monotone_subcurve_2>(&(*it));
#if ! defined (CGAL_NO_ASSERTIONS)
      CGAL_assertion(x_seg_p != nullptr);
#endif

      // If the polycurve consists of a single x-monotone subcurve, return.
      if (x_seg_objects.size() == 1) {
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        if (cmp_seg_endpts(x_seg) == LARGER)
          *oi++ = Make_x_monotone_result(ctr_x_curve(ctr_seg_opposite(*x_seg_p)));        else *oi++ = Make_x_monotone_result(ctr_x_curve(*x_seg_p));
#else
        *oi++ = Make_x_monotone_result(ctr_x_curve(*x_seg_p));
#endif
        x_seg_objects.clear();
        return oi;
      }

      CGAL_precondition_code
        (
         // To be used in order to verify continuity and well-orientedness
         // of the input curve cv.
         auto min_seg_v =
           m_poly_traits.subcurve_traits_2()->construct_min_vertex_2_object();
         auto max_seg_v =
           m_poly_traits.subcurve_traits_2()->construct_max_vertex_2_object();
         auto equal = m_poly_traits.subcurve_traits_2()->equal_2_object();
         Point_2 last_target = (cmp_seg_endpts(*x_seg_p) == SMALLER) ?
           max_seg_v(*x_seg_p) : min_seg_v(*x_seg_p);
         Point_2 next_src;
         );

      // The polycurve consists of at least 2 x-monotone subcurves:
      Push_back_2 push_back = m_poly_traits.push_back_2_object();
      auto is_seg_vertical =
        m_poly_traits.subcurve_traits_2()->is_vertical_2_object();

      bool is_start_vertical = is_seg_vertical(*x_seg_p);
      Comparison_result start_dir = cmp_seg_endpts(*x_seg_p);

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      Push_front_2 push_front = m_poly_traits.push_front_2_object();
      X_monotone_curve_2 x_polycurve = (cmp_seg_endpts(x_seg) == LARGER) ?
        ctr_x_curve(ctr_seg_opposite(*x_seg_p)) : ctr_x_curve(*x_seg_p);
#else
      X_monotone_curve_2 x_polycurve = ctr_x_curve(*x_seg_p);
#endif

      for (++it; it != x_seg_objects.end(); ++it) {
        const auto* x_seg_p = std::get_if<X_monotone_subcurve_2>(&(*it));
#if ! defined (CGAL_NO_ASSERTIONS)
        CGAL_assertion(x_seg_p != nullptr);
#endif

        // Test that cv is continuous and well-oriented.
        CGAL_precondition_code
          (
           next_src = (cmp_seg_endpts(*x_seg_p) == SMALLER) ?
             min_seg_v(*x_seg_p) : max_seg_v(*x_seg_p);
           );
        CGAL_precondition_msg
          (
           equal(last_target, next_src),
             "cv must form a continuous and well oriented curve."
           );
        CGAL_precondition_code
          (
           last_target = (cmp_seg_endpts(*x_seg_p) == SMALLER) ?
             max_seg_v(*x_seg_p) : min_seg_v(*x_seg_p);
           );

        if ((cmp_seg_endpts(*x_seg_p) != start_dir) ||
            (is_seg_vertical(*x_seg_p) != is_start_vertical))
        {
            // Construct an x-monotone curve from the sub-range which was found
          *oi++ = Make_x_monotone_result(x_polycurve);
          is_start_vertical = is_seg_vertical(*x_seg_p);
          start_dir = cmp_seg_endpts(*x_seg_p);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
          x_polycurve = (cmp_seg_endpts(*x_seg_p) == LARGER) ?
            ctr_x_curve(ctr_seg_opposite(*x_seg_p)) : ctr_x_curve(*x_seg_p);
#else
          x_polycurve = ctr_x_curve(*x_seg_p);
#endif
        }
        else {
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
          if (cmp_seg_endpts(*x_seg_p) == LARGER)
            push_front(x_polycurve, ctr_seg_opposite(*x_seg_p));
          else
            push_back(x_polycurve, *x_seg_p);
#else
          push_back(x_polycurve, *x_seg_p);
#endif
        }

      } // for loop
      if (x_polycurve.number_of_subcurves() != 0)
        *oi++ = Make_x_monotone_result(x_polycurve);
      x_seg_objects.clear();
      return oi;
    }

    //!
    template <typename OutputIterator>
    OutputIterator operator_impl(const Curve_2& cv, OutputIterator oi,
                                 Arr_not_all_sides_oblivious_tag) const {
      using Make_x_monotone_subresult =
        std::variant<Point_2, X_monotone_subcurve_2>;
      using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;

      // If the polycurve is empty, return.
      if (cv.number_of_subcurves() == 0) return oi;

      auto ctr_x_curve = m_poly_traits.construct_x_monotone_curve_2_object();

      auto make_seg_x_monotone =
        m_poly_traits.subcurve_traits_2()->make_x_monotone_2_object();

      auto cmp_seg_endpts =
        m_poly_traits.subcurve_traits_2()->compare_endpoints_xy_2_object();

      auto ps_x =
        m_poly_traits.subcurve_traits_2()->parameter_space_in_x_2_object();
      auto ps_y =
        m_poly_traits.subcurve_traits_2()->parameter_space_in_y_2_object();

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      typename Subcurve_traits_2::Construct_opposite_2 ctr_seg_opposite =
        m_poly_traits.subcurve_traits_2()->construct_opposite_2_object();
#endif

      // Convert the input polycurve to a sequence of objects, such that
      // each object wraps an x-monotone subcurve.
      std::vector<Make_x_monotone_subresult> x_seg_objects;
      for (auto its = cv.subcurves_begin(); its != cv.subcurves_end(); ++its)
        make_seg_x_monotone(*its, std::back_inserter(x_seg_objects));
      auto it = x_seg_objects.begin();
      const auto* x_seg_p = std::get_if<X_monotone_subcurve_2>(&(*it));
#if ! defined (CGAL_NO_ASSERTIONS)
      CGAL_assertion(x_seg_p != nullptr);
#endif

      // If the polycurve consists of a single x-monotone subcurve, return.
      if (x_seg_objects.size() == 1) {
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        if (cmp_seg_endpts(x_seg) == LARGER)
          *oi++ = Make_x_monotone_result(ctr_x_curve(ctr_seg_opposite(*x_seg_p)));        else *oi++ = Make_x_monotone_result(ctr_x_curve(*x_seg_p));
#else
        *oi++ = Make_x_monotone_result(ctr_x_curve(*x_seg_p));
#endif
        x_seg_objects.clear();
        return oi;
      }

      CGAL_precondition_code
        (
         // To be used in order to verify continuity and well-orientedness
         // of the input curve cv.
         auto min_seg_v =
           m_poly_traits.subcurve_traits_2()->construct_min_vertex_2_object();
         auto max_seg_v =
           m_poly_traits.subcurve_traits_2()->construct_max_vertex_2_object();
         auto equal = m_poly_traits.subcurve_traits_2()->equal_2_object();
         Point_2 last_target = (cmp_seg_endpts(*x_seg_p) == SMALLER) ?
           max_seg_v(*x_seg_p) : min_seg_v(*x_seg_p);
         Point_2 next_src;
         );

      // The polycurve consists of at least 2 x-monotone subcurves:
      Push_back_2 push_back = m_poly_traits.push_back_2_object();
      auto is_seg_vertical =
        m_poly_traits.subcurve_traits_2()->is_vertical_2_object();

      bool is_start_vertical = is_seg_vertical(*x_seg_p);
      Comparison_result start_dir = cmp_seg_endpts(*x_seg_p);

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      Push_front_2 push_front = m_poly_traits.push_front_2_object();
      X_monotone_curve_2 x_polycurve = (cmp_seg_endpts(x_seg) == LARGER) ?
        ctr_x_curve(ctr_seg_opposite(*x_seg_p)) : ctr_x_curve(*x_seg_p);
#else
      X_monotone_curve_2 x_polycurve = ctr_x_curve(*x_seg_p);
#endif

      for (++it; it != x_seg_objects.end(); ++it) {
        const auto* x_seg_p = std::get_if<X_monotone_subcurve_2>(&(*it));
#if ! defined (CGAL_NO_ASSERTIONS)
        CGAL_assertion(x_seg_p != nullptr);
#endif

        // Test that cv is continuous and well-oriented.
        CGAL_precondition_code
          (
           next_src = (cmp_seg_endpts(*x_seg_p) == SMALLER) ?
             min_seg_v(*x_seg_p) : max_seg_v(*x_seg_p);
           );
        CGAL_precondition_msg
          (
           equal(last_target, next_src),
             "cv must form a continuous and well oriented curve."
           );
        CGAL_precondition_code
          (
           last_target = (cmp_seg_endpts(*x_seg_p) == SMALLER) ?
             max_seg_v(*x_seg_p) : min_seg_v(*x_seg_p);
           );

        Arr_curve_end polycurve_target =
          (cmp_seg_endpts(x_polycurve[0]) == SMALLER) ?
          ARR_MAX_END : ARR_MIN_END;
        Arr_curve_end seg_source = (cmp_seg_endpts(*x_seg_p) == SMALLER) ?
          ARR_MIN_END : ARR_MAX_END;
        auto num_segs = x_polycurve.number_of_subcurves();

        if ((cmp_seg_endpts(*x_seg_p) != start_dir) ||
            (is_seg_vertical(*x_seg_p) != is_start_vertical))
        {
            // Construct an x-monotone curve from the sub-range which was found
          *oi++ = Make_x_monotone_result(x_polycurve);
          is_start_vertical = is_seg_vertical(*x_seg_p);
          start_dir = cmp_seg_endpts(*x_seg_p);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
          x_polycurve (cmp_seg_endpts(*x_seg_p) == LARGER) ?
            ctr_x_curve(ctr_seg_opposite(*x_seg_p)) : ctr_x_curve(*x_seg_p);
#else
          x_polycurve = ctr_x_curve(*x_seg_p);
#endif
        }
        else if (ps_x(x_polycurve[num_segs-1], polycurve_target) !=
                 ARR_INTERIOR ||
                 (ps_x(*x_seg_p, seg_source) != ARR_INTERIOR))
        {
          *oi++ = Make_x_monotone_result(x_polycurve);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
          x_polycurve = (cmp_seg_endpts(*x_seg_p) == LARGER) ?
            ctr_seg_opposite(*x_seg_p) : ctr_x_curve(*x_seg_p);
#endif
          x_polycurve = ctr_x_curve(*x_seg_p);
        }

        else {
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
          if (cmp_seg_endpts(*x_seg_p) == LARGER)
            push_front(x_polycurve, ctr_seg_opposite(*x_seg_p));
          else
            push_back(x_polycurve, *x_seg_p);
#else
          push_back(x_polycurve, *x_seg_p);
#endif
        }
      } // for loop
      if (x_polycurve.number_of_subcurves() != 0)
        *oi++ = Make_x_monotone_result(x_polycurve);
      x_seg_objects.clear();
      return oi;
    }

  public:
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    { return operator_impl(cv, oi, All_sides_oblivious_category()); }
  };

  /*! obtains a `Make_x_monotone_2` functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(*this); }

  /* Functor to augment a polycurve by either adding a vertex or a subcurve
   * at the back.
   * TODO: Test all the operator()'s. (Don't forget vertical cases!)
   */
  class Push_back_2 : public Base::Push_back_2 {
  protected:
    using Polycurve_traits_2 = Arr_polycurve_traits_2<Subcurve_traits_2>;

  public:
    /*! constructs. */
    Push_back_2(const Polycurve_traits_2& traits) : Base::Push_back_2(traits) {}

    // Normally, the moment the compiler finds a name, it stops looking. In
    // other words, the compiler first finds the operator() in the current
    // class and stops looking, never finding the one in the base class.
    // Explicitly bring the base operator() into scope, unnecesitating the
    // code below.
    using Base::Push_back_2::operator();

    // /*! appends a subcurve to an existing x-monotone polycurve at the back.
    //  */
    // void operator()(X_monotone_curve_2& xcv,
    //                 const X_monotone_subcurve_2& seg)
    //   const
    // { Base::Push_back_2::operator()(xcv, seg); }

    /* appends a subcurve to an existing polycurve at the back.
     * If the polycurve is empty, the subcurve will be its only subcurve.
     */
    void operator()(Curve_2& cv, const Subcurve_2& seg) const
    { cv.push_back(seg); }
  };

  /*! obtains a `Push_back_2` functor object. */
  Push_back_2 push_back_2_object() const { return Push_back_2(*this); }

  /* Functor to augment a polycurve by either adding a vertex or a subcurve
   * at the front.
   * TODO: Test all the operator()'s. (Don't forget vertical cases!)
   */
  class Push_front_2 : public Base::Push_front_2 {
  protected:
    using Polycurve_traits_2 = Arr_polycurve_traits_2<Subcurve_traits_2>;

  public:
    /*! constructs. */
    Push_front_2(const Polycurve_traits_2& traits) :
      Base::Push_front_2(traits)
    {}

    // Normally, the moment the compiler finds a name, it stops looking. In
    // other words, the compiler first finds the operator() in the current
    // class and stops looking, never finding the one in the base class.
    // Explicitly bring the base operator() into scope, unnecesitating the
    // code below.
    using Base::Push_front_2::operator();

    // /*! appends a subcurve to an existing x-monotone polycurve at the front.
    //  */
    // void operator()(X_monotone_curve_2& xcv,
    //                 const X_monotone_subcurve_2& seg)
    //   const
    // { Base::Push_front_2::operator()(xcv, seg); }

    /* appends a subcurve to an existing polycurve at the front. */
    void operator()(Curve_2& cv, const Subcurve_2& seg) const
    { cv.push_front(seg); }
  };

  /*! obtains a `Push_front_2` functor object. */
  Push_front_2 push_front_2_object() const { return Push_front_2(*this); }

  class Split_2 {
  protected:
    using Polycurve_traits_2 = Arr_polycurve_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state)
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! constructs. */
    Split_2(const Polycurve_traits_2& traits) : m_poly_traits(traits) {}

  public:
    /*! splits a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve(p is its right endpoint).
     * \param c2 Output: The right resulting subcurve(p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& xcv, const Point_2& p,
                    X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      auto min_vertex = geom_traits->construct_min_vertex_2_object();
      auto max_vertex = geom_traits->construct_max_vertex_2_object();
      auto equal = geom_traits->equal_2_object();
      auto cmp_seg_endpts = geom_traits->compare_endpoints_xy_2_object();

      // Make sure the split point is not one of the curve endpoints.
      CGAL_precondition((! equal(m_poly_traits.
                                 construct_min_vertex_2_object()(xcv), p)));
      CGAL_precondition((! equal(m_poly_traits.
                                 construct_max_vertex_2_object()(xcv), p)));

      CGAL_precondition_msg(xcv.number_of_subcurves() > 0,
                            "Cannot split a polycurve of length zero.");

      Comparison_result dir = cmp_seg_endpts(xcv[0]);

      // Locate the subcurve on the polycurve xcv that contains p.
      auto i = m_poly_traits.locate_impl(xcv, p, All_sides_oblivious_category());
      CGAL_precondition(i != Polycurve_traits_2::INVALID_INDEX);

      // Clear the output curves.
      xcv1.clear();
      xcv2.clear();

      // Push all subcurves labeled(0, 1, ... , i-1) into xcv1.
      for (std::size_t j = 0; j < i; ++j) xcv1.push_back(xcv[j]);

      if (dir == SMALLER){
        // Check whether the split point is xcv[i]'s source or target.
        if (equal(max_vertex(xcv[i]), p)) {
          // The entire i'th subcurve belongs to xcv1:
          xcv1.push_back(xcv[i]);
        }
        else if (equal(min_vertex(xcv[i]), p)) {
          // The entire i'th subcurves belongs to xcv2:
          xcv2.push_back(xcv[i]);
        }
        else {
          // The i'th subcurve should be split: The left part(seg1)
          // goes to xcv1, and the right part(seg2) goes to xcv2.
          X_monotone_subcurve_2 seg1, seg2;
          m_poly_traits.subcurve_traits_2()->split_2_object()(xcv[i], p,
                                                              seg1, seg2);

          xcv1.push_back(seg1);
          xcv2.push_back(seg2);
        }
      }
      else {
        if (equal(min_vertex(xcv[i]), p)) {
          xcv1.push_back(xcv[i]);
        }
        else if (equal(max_vertex(xcv[i]), p)) {
          xcv2.push_back(xcv[i]);
        }
        else {
          X_monotone_subcurve_2 seg1, seg2;
          m_poly_traits.subcurve_traits_2()->
            split_2_object()(xcv[i], p, seg1, seg2);

          if (cmp_seg_endpts(seg2) == LARGER){
            xcv1.push_back(seg2);
          }
          else {
            // seg2 has to be reversed
            seg2 = m_poly_traits.subcurve_traits_2()->
              construct_opposite_2_object()(seg2);
            xcv1.push_back(seg2);
          }

          if (cmp_seg_endpts(seg1) == LARGER){
            xcv2.push_back(seg1);
          }
          else {
            // seg2 has to be reversed
            seg1 = m_poly_traits.subcurve_traits_2()->
              construct_opposite_2_object()(seg1);
            xcv1.push_back(seg1);
          }
        }
      }

      // Push all subcurves labeled(i+1, i+2, ... , n-1) into xcv1.
      std::size_t n = xcv.number_of_subcurves();

      for (std::size_t j = i+1; j < n; ++j) xcv2.push_back(xcv[j]);

      if (dir != SMALLER) std::swap(xcv1, xcv2);
    }
  };

  /*! obtains a `Split_2` functor object. */
  Split_2 split_2_object() const { return Split_2(*this); }

  class Intersect_2 {
  protected:
    using Polycurve_traits_2 = Arr_polycurve_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state)
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! constructs. */
    Intersect_2(const Polycurve_traits_2& traits) : m_poly_traits(traits) {}

    /*! finds the intersections of the two given curves and insert them into the
     * given output iterator. As two subcurves may itersect only once, only a
     * single intersection will be contained in the iterator.
     * Note: If the intersection yields an overlap then it will be oriented
     *       from left-to-right.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const {
      // std::cout << "intersect(" << cv1 << ", " << cv2 << ")\n";
      using Intersection_point = std::pair<Point_2, Multiplicity>;
      using Intersection_base_result =
        std::variant<Intersection_point, X_monotone_subcurve_2>;

      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      auto cmp_y_at_x = m_poly_traits.compare_y_at_x_2_object();
      auto equal = geom_traits->equal_2_object();
      auto min_vertex = geom_traits->construct_min_vertex_2_object();
      auto max_vertex = geom_traits->construct_max_vertex_2_object();
      auto intersect = geom_traits->intersect_2_object();
      auto cmp_endpts = geom_traits->compare_endpoints_xy_2_object();
      auto ctr_opposite = geom_traits->construct_opposite_2_object();
      auto cmp_xy = m_poly_traits.compare_xy_2_object();

      Comparison_result dir1 = cmp_endpts(cv1[0]);
      Comparison_result dir2 = cmp_endpts(cv2[0]);
      // std::cout << "dir1: " << dir1 << std::endl;
      // std::cout << "dir2: " << dir2 << std::endl;

      const bool consistent = (dir1 == dir2);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      CGAL_assertion(consistent);
#endif
      std::vector<X_monotone_subcurve_2> ocv; // Used to represent overlaps.
      const bool invert_ocv = ((dir1 == LARGER) && (dir2 == LARGER));

      const std::size_t n1 = cv1.number_of_subcurves();
      const std::size_t n2 = cv2.number_of_subcurves();

      std::size_t i1 = (dir1 == SMALLER) ? 0 : n1-1;
      std::size_t i2 = (dir2 == SMALLER) ? 0 : n2-1;

      Point_2 saved_point;
      X_monotone_curve_2 saved_xcv;
      bool point_saved = false;
      bool xcv_saved = false;

      // Save a point
      auto update_saved_point = [&](const Point_2& p) {
        point_saved = true;
        saved_point = p;
      };

      // Spit out the saved point
      auto spit_saved_point = [&]() {
        *oi++ = std::make_pair(saved_point, 0);
        point_saved = false;
      };

      // Add a subcurve to the saved x-monotone curve
      auto update_saved_xcv = [&](const X_monotone_subcurve_2& sxcv) {
        xcv_saved = true;
        // We maintain the invariant that if the input curves have opposite
        // directions (! consistent), the overalpping curves are directed
        // left=>right. This, however, is not guaranteed by all traits for the
        // subcurves. (It is guaranteed by some traits, but this is
        // insufficient.)  Therefore, we need to enforce it. That is, we make
        // sure the subcurves are also directed left=>right in this case.
        if (! consistent && (cmp_endpts(sxcv) == LARGER)) {
          if (invert_ocv) {
            saved_xcv.push_front(ctr_opposite(sxcv));
            return;
          }
          saved_xcv.push_back(ctr_opposite(sxcv));
          return;
        }
        if (invert_ocv) {
          saved_xcv.push_front(sxcv);
          return;
        }
        saved_xcv.push_back(sxcv);
      };

      // Spit out the saved x-monotone curve
      auto spit_saved_xcv = [&]() {
        *oi++ = saved_xcv;
        saved_xcv.clear();
        xcv_saved = false;
      };

      auto left_res = cmp_xy(cv1[i1], ARR_MIN_END, cv2[i2], ARR_MIN_END);
      if (left_res == SMALLER) {
        // cv1's left endpoint is to the left of cv2's left endpoint.
        // Locate the index i1 of the subcurve in cv1 which contains cv2's
        // left endpoint.
        i1 = m_poly_traits.locate_impl(cv1, cv2[i2], ARR_MIN_END,
                                       All_sides_oblivious_category());
        if (i1 == Polycurve_traits_2::INVALID_INDEX) return oi;
        if (cmp_y_at_x(cv2[i2], ARR_MIN_END, cv1[i1]) == EQUAL) {
          const auto& p = min_vertex(cv2[i2]);
          update_saved_point(p);
          if (equal(max_vertex(cv1[i1]), p)) {
            if (dir1 == SMALLER) {
              ++i1;
              if (i1 == n1) {
                spit_saved_point();
                return oi;
              }
            }
            else {
              if (i1 != 0) --i1;
              else {
                spit_saved_point();
                return oi;
              }
            }
          }
        }
      }
      else if (left_res == LARGER) {
        // cv1's left endpoint is to the right of cv2's left endpoint.
        // Locate the index i2 of the subcurve in cv2 which contains cv1's
        // left endpoint.
        i2 = m_poly_traits.locate_impl(cv2, cv1[i1], ARR_MIN_END,
                                       All_sides_oblivious_category());
        if (i2 == Polycurve_traits_2::INVALID_INDEX) return oi;
        if (cmp_y_at_x(cv1[i1], ARR_MIN_END, cv2[i2]) == EQUAL) {
          const auto& p = min_vertex(cv1[i1]);
          update_saved_point(p);
          if (equal(max_vertex(cv2[i2]), p)) {
            if (dir2 == SMALLER) {
              ++i2;
              if (i2 == n2) {
                spit_saved_point();
                return oi;
              }
            }
            else {
              if (i2 != 0) --i2;
              else {
                spit_saved_point();
                return oi;
              }
            }
          }
        }
      }
      else {
        CGAL_assertion(left_res == EQUAL);
        update_saved_point(min_vertex(cv1[i1]));
      }

      do {
        // std::cout << "i1, i2 = " << i1 <<", " << i2 << std::endl;
        std::vector<Intersection_base_result> xsecs;
        intersect(cv1[i1], cv2[i2], std::back_inserter(xsecs));

        // Iterate over all intersections.
        if (! xsecs.empty()) {
          auto it = xsecs.begin();
          const auto& xsec = *it++;
          if (it == xsecs.end()) {
            // There is exactly one intersection
            const auto* xp_p = std::get_if<Intersection_point>(&xsec);
            if (xp_p) {
              // The intersection is a point
              const auto& xp = xp_p->first;
              auto is_min_end_1 = equal(xp, min_vertex(cv1[i1]));
              auto is_min_end_2 = equal(xp, min_vertex(cv2[i2]));
              if ((is_min_end_1 || is_min_end_2) && (xcv_saved || point_saved)) {
                // It is impossible to have a pending point and a pending
                // x-monotone curve concurrently
                CGAL_assertion((xcv_saved && ! point_saved) ||
                               (! xcv_saved && point_saved));
                if (xcv_saved) spit_saved_xcv();
                else spit_saved_point();
              }
              else {
                auto is_max_end_1 = equal(xp, max_vertex(cv1[i1]));
                auto is_max_end_2 = equal(xp, max_vertex(cv2[i2]));
                if (is_max_end_1 || is_max_end_2) update_saved_point(xp);
                else *oi++ = *xp_p;
              }
            }
            else {
              // The intersection is an x-monotone curve
              auto* sxcv_p = std::get_if<X_monotone_subcurve_2>(&xsec);
              CGAL_assertion(sxcv_p != nullptr);
              const auto& xp_min = min_vertex(*sxcv_p);
              auto is_min_end_1 = equal(xp_min, min_vertex(cv1[i1]));
              auto is_min_end_2 = equal(xp_min, min_vertex(cv2[i2]));
              const auto& xp_max = max_vertex(*sxcv_p);
              auto is_max_end_1 = equal(xp_max, max_vertex(cv1[i1]));
              auto is_max_end_2 = equal(xp_max, max_vertex(cv2[i2]));
              if (is_min_end_1 || is_min_end_2) {
                update_saved_xcv(*sxcv_p);
                if (! is_max_end_1 && ! is_max_end_2) spit_saved_xcv();
              }
              else {
                if (xcv_saved) spit_saved_xcv();
                if (is_max_end_1 || is_max_end_2) update_saved_xcv(*sxcv_p);
                else *oi++ = *sxcv_p;
              }
              point_saved = false;
            }
          }
          else {
            // There is more than one intersection
            // Handle the first intersection
            const auto* xp_p = std::get_if<Intersection_point>(&xsec);
            if (xp_p) {
              const auto& xp = xp_p->first;
              auto is_min_end_1 = equal(xp, min_vertex(cv1[i1]));
              auto is_min_end_2 = equal(xp, min_vertex(cv2[i2]));
              if ((is_min_end_1 || is_min_end_2) && (xcv_saved || point_saved)) {
                // It is impossible to have a pending point and a pending
                // x-monotone curve concurrently
                CGAL_assertion((xcv_saved && ! point_saved) ||
                               (! xcv_saved && point_saved));
                if (xcv_saved) spit_saved_xcv();
                else spit_saved_point();
              }
              else *oi++ = *xp_p;
            }
            else {
              const auto* sxcv_p = std::get_if<X_monotone_subcurve_2>(&xsec);
              CGAL_assertion(sxcv_p);
              const auto& xp = min_vertex(*sxcv_p);
              auto is_min_end_1 = equal(xp, min_vertex(cv1[i1]));
              auto is_min_end_2 = equal(xp, min_vertex(cv2[i2]));
              if (is_min_end_1 || is_min_end_2) {
                update_saved_xcv(*sxcv_p);
                spit_saved_xcv();
              }
              else {
                if (xcv_saved) spit_saved_xcv();
                else *oi++ = *sxcv_p;
              }
              point_saved = false;
            }

            // Handle all but the last intersection
            for (; it != std::prev(xsecs.end()); ++it) {
              const auto& xsec = *it;
              xp_p = std::get_if<Intersection_point>(&xsec);
              if (xp_p) *oi++ = *xp_p;
              else {
                const auto* sxcv_p = std::get_if<X_monotone_subcurve_2>(&xsec);
                CGAL_assertion(sxcv_p);
                *oi++ = *sxcv_p;
              }
            }

            // Handle the last intersection
            const auto& xsec = *it;
            xp_p = std::get_if<Intersection_point>(&xsec);
            if (xp_p) {
              const auto& xp = xp_p->first;
              auto is_max_end_1 = equal(xp, max_vertex(cv1[i1]));
              auto is_max_end_2 = equal(xp, max_vertex(cv2[i2]));
              if (is_max_end_1 || is_max_end_2) update_saved_point(xp);
              else *oi++ = *xp_p;
            }
            else {
              const auto* sxcv_p = std::get_if<X_monotone_subcurve_2>(&xsec);
              CGAL_assertion(sxcv_p);
              const auto& xp = max_vertex(*sxcv_p);
              auto is_max_end_1 = equal(xp, max_vertex(cv1[i1]));
              auto is_max_end_2 = equal(xp, max_vertex(cv2[i2]));
              if (is_max_end_1 || is_max_end_2) update_saved_xcv(*sxcv_p);
              else *oi++ = *sxcv_p;
            }
          }
        }

        // Advance the indices
        auto right_res = cmp_xy(cv1[i1], ARR_MAX_END, cv2[i2], ARR_MAX_END);
        if (right_res != LARGER) {
          if (dir1 == SMALLER) {
            ++i1;
            if (i1 == n1) break;
          }
          else {
            if (i1 != 0) --i1;
            else break;
          }
        }
        if (right_res != SMALLER) {
          if (dir2 == SMALLER) {
            ++i2;
            if (i2 == n2) break;
          }
          else {
            if (i2 != 0) --i2;
            else break;
          }
        }
      } while (true);

      if (xcv_saved) spit_saved_xcv();
      else if (point_saved) spit_saved_point();

      return oi;
    }

  private:
    template <typename OutputIterator>
    inline OutputIterator output_ocv(std::vector<X_monotone_subcurve_2>& ocv,
                                     bool invert_ocv, OutputIterator oi) const {
      X_monotone_curve_2 curve;
      if (invert_ocv) std::reverse(ocv.begin(), ocv.end());
      for (X_monotone_subcurve_2& sc : ocv) curve.push_back(sc);
      *oi++ = curve;
      ocv.clear();
      return oi;
    }
  };

  /*! obtains an `Intersect_2` functor object. */
  Intersect_2 intersect_2_object() const { return Intersect_2(*this); }

  class Are_mergeable_2 {
  protected:
    using Polycurve_traits_2 = Arr_polycurve_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state)
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! constructs. */
    Are_mergeable_2(const Polycurve_traits_2& traits) : m_poly_traits(traits) {}

    /*! checks whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return(true) if the two curves are mergeable, that is, they share a
     * common endpoint and the same orientation;(false) otherwise.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const {
      const Subcurve_traits_2* geom_traits = m_poly_traits.subcurve_traits_2();
      auto min_vertex = m_poly_traits.construct_min_vertex_2_object();
      auto max_vertex = m_poly_traits.construct_max_vertex_2_object();
      auto equal = geom_traits->equal_2_object();
      auto is_seg_vertical = geom_traits->is_vertical_2_object();

      auto dir1 = m_poly_traits.compare_endpoints_xy_2_object()(cv1);
      auto dir2 = m_poly_traits.compare_endpoints_xy_2_object()(cv2);

      if (dir1 != dir2) return false;

      bool ver1 = is_seg_vertical(cv1[0]);
      bool ver2 = is_seg_vertical(cv2[0]);

      return (((// Both are directed from left-to-right
                (dir1 == SMALLER) &&
                ((equal(max_vertex(cv1),min_vertex(cv2))) ||
                 (equal(max_vertex(cv2),min_vertex(cv1))))) ||
               (// Both are directed from right-to-left
                (dir1 == LARGER) &&
                ((equal(min_vertex(cv1),max_vertex(cv2))) ||
                 (equal(max_vertex(cv1),min_vertex(cv2))))
                )) &&
              (// Either both should be vertical or both should
               // be NOT vertical.
               (ver1 && ver2) || (!ver1 && !ver2)));
    }
  };

  /*! obtains an `Are_mergeable_2` functor object. */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(*this); }

  /*! \class Merge_2
   * A functor that merges two x-monotone curves into one.
   */
  /* Roadmap: Allow merging of overlapping polycurves. This means also
   *          changing the subcurve traits class.
   */
  class Merge_2 {
  protected:
    using Geometry_traits = Arr_polycurve_traits_2<Subcurve_traits_2>;

    //! The traits (in case it has state)
    const Geometry_traits& m_poly_traits;

  public:
    /*! constructs
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Geometry_traits& traits) : m_poly_traits(traits) {}

    /*! merges two given x-monotone curves into a single curve(segment).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable.
     */
    void operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const {
      CGAL_precondition(m_poly_traits.are_mergeable_2_object()(cv1, cv2));

      auto min_vertex = m_poly_traits.construct_min_vertex_2_object();
      auto max_vertex = m_poly_traits.construct_max_vertex_2_object();
      auto cmp_endpts = m_poly_traits.compare_endpoints_xy_2_object();
      Equal_2 equal = m_poly_traits.equal_2_object();

      c.clear();
      if (// Either both are left-to-right and cv2 is to the right of cv1
          ((cmp_endpts(cv1)==SMALLER) &&
           (equal(max_vertex(cv1), min_vertex(cv2)))) ||
          // or both are right-to-left and cv2 is to the left of cv1
          ((cmp_endpts(cv1)==LARGER) &&
           (equal(min_vertex(cv1), max_vertex(cv2)))))
      {
        const std::size_t n1 = cv1.number_of_subcurves();
        const std::size_t n2 = cv2.number_of_subcurves();
        std::size_t i;

        // cv2 extends cv1 to the right:
        for (i = 0; i < n1 - 1; ++i) c.push_back(cv1[i]);

        // Try to merge the to contiguous line subcurves:
        if (m_poly_traits.subcurve_traits_2()->
            are_mergeable_2_object()(cv1[n1 - 1], cv2[0]))
        {
          X_monotone_subcurve_2 seg;
          m_poly_traits.subcurve_traits_2()->
            merge_2_object()(cv1[n1 - 1], cv2[0], seg);
          c.push_back(seg);
        }
        else {
          c.push_back(cv1[n1 - 1]);
          c.push_back(cv2[0]);
        }

        for (i = 1; i < n2; ++i) c.push_back(cv2[i]);
      }
      else
        return this->operator()(cv2,cv1,c);
    }
  };

  /*! obtains a `Merge_2` functor object. */
  Merge_2 merge_2_object() const { return Merge_2(*this); }
  ///@}

  /*! \class
   * A functor that constructs a (general) polycurve.
   */
  class Construct_curve_2 {
  protected:
    using Polycurve_traits_2 = Arr_polycurve_traits_2<Subcurve_traits_2>;

    //! The polycurve traits (in case it has state)
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! constructs. */
    Construct_curve_2(const Polycurve_traits_2& traits) :
      m_poly_traits(traits)
    {}

    /*! obtains a polycurve that consists of one given subcurve. */
    Curve_2 operator()(const Subcurve_2& seg) const { return Curve_2(seg); }

    /* constructs a well-oriented polycurve from a range of either
     * `SubcurveTraits::Point_2` or `SubcurveTraits::Subcurve_2`.
     */
    template <typename ForwardIterator>
    Curve_2 operator()(ForwardIterator begin, ForwardIterator end) const {
      using VT = typename std::iterator_traits<ForwardIterator>::value_type;
      using Is_point = typename std::is_same<VT, Point_2>::type;
      // Dispatch the range to the appropriate implementation.
      return constructor_impl(begin, end, Is_point());
    }

    /*! constructs a polycurve from a range of points.
     * \pre The range contains at least two points
     * \pre Consecutive points are disjoint.
     * \return Well-oriented polycurve connecting the given
     *         points. The order of the vertices is determined by
     *         their order in the range.  Furthermore, the
     *         orientation of the polycurve is induced by their
     *         order.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl(ForwardIterator /* begin */,
                             ForwardIterator /* end */,
                             std::true_type) const
    {  CGAL_error_msg("Cannot construct a polycurve from a range of points!"); }

    /*! Construction implementation from a range of subcurves.
     *  Note that the subcurves in the range are NOT necessarily x-monotone,
     *  thus it is impossible to test (even in precondition) whether the input
     *  forms a continuous and well oriented polycurve.
     *  \pre Range should contain at least one subcurve.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end,
                             std::false_type) const {
      // Range has to contain at least one subcurve
      CGAL_precondition(begin != end);
      return Curve_2(begin, end);
    }
  };

  /*! obtains a `Construct_curve_2` functor object. */
  Construct_curve_2 construct_curve_2_object() const
  { return Construct_curve_2(*this); }
};

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
