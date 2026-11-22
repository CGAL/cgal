// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University(Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Efi Fogel <efif@post.tau.ac.il>
//             Ron Wein  <wein@post.tau.ac.il>
//             Dror Atariah <dror.atariah@fu-berlin.de>
//             Waqar Khan <wkhan@mpi-inf.mpg.de>

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
#include <utility>
#include <variant>
#include <vector>

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polycurve_basic_traits_2.h>
#include <CGAL/Arr_geometry_traits/Polycurve_2.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>

namespace CGAL {

template <typename SubcurveTraits_2 = Arr_segment_traits_2<>>
class Arr_polycurve_traits_2 : public Arr_polycurve_basic_traits_2<SubcurveTraits_2> {
public:
  using Subcurve_traits_2 = SubcurveTraits_2;

private:
  using Base = Arr_polycurve_basic_traits_2<Subcurve_traits_2>;

public:
  /// \name Types inherited from the polycurve basic traits class.
  //@{
  using Has_left_category = typename Base::Has_left_category;

  using Left_side_category = typename Base::Left_side_category;
  using Bottom_side_category = typename Base::Bottom_side_category;
  using Top_side_category = typename Base::Top_side_category;
  using Right_side_category = typename Base::Right_side_category;

  using All_sides_oblivious_category = typename Base::All_sides_oblivious_category;

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
  using Construct_x_monotone_curve_2 = typename Base::Construct_x_monotone_curve_2;
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
   * \param seg_traits an already existing subcurve traits which is passed will
   *        be used by the class.
   */
  Arr_polycurve_traits_2(const Subcurve_traits_2* geom_traits) : Base(geom_traits) {}

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
  Number_of_points_2 number_of_points_2_object() const { return Number_of_points_2(); }

  ///@}

  /// \name Construction functors(based on the subcurve traits).
  //@{

  //! A functor for subdividing curves into x-monotone curves.
  class Make_x_monotone_2 {
  protected:
    using Polycurve_traits_2 = Arr_polycurve_traits_2<Subcurve_traits_2>;

    //! The traits (in case it has state)
    const Polycurve_traits_2& m_poly_traits;

  public:
    /*! constructs. */
    Make_x_monotone_2(const Polycurve_traits_2& traits) : m_poly_traits(traits) {}

    /*! subdivides a given curve into x-monotone sub-curves and insert them into
     * a given output iterator.
     *
     * \pre if `cv` is not empty then it must be continuous and well-oriented.
     * \param cv the curve.
     * \param oi an output iterator for the result. Its value type is a variant
     *           that wraps Point_2 or an X_monotone_curve_2 objects.
     * \return the past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    { return operator_impl(cv, oi, All_sides_oblivious_category()); }

  private:
    template <typename OutputIterator, typename SidesCategory>
    OutputIterator operator_impl(const Curve_2& cv, OutputIterator oi, SidesCategory) const {
      using Make_x_monotone_subresult = std::variant<Point_2, X_monotone_subcurve_2>;
      using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;
      constexpr bool oblivious = std::is_same_v<SidesCategory, Arr_all_sides_oblivious_tag>;

      const auto* geom_traits = m_poly_traits.subcurve_traits_2();
      auto make_seg_x_monotone = geom_traits->make_x_monotone_2_object();
      auto cmp_seg_endpts = geom_traits->compare_endpoints_xy_2_object();
      auto is_seg_vertical = geom_traits->is_vertical_2_object();

      CGAL_precondition_code
        (
         // To be used in order to verify that the subcurves are not
         // degenerate, and that cv is continuous and well-oriented.
         auto min_seg_v = geom_traits->construct_min_vertex_2_object();
         auto max_seg_v = geom_traits->construct_max_vertex_2_object();
         auto equal = geom_traits->equal_2_object();
         Point_2 last_target;
         bool is_first = true;
         );

      // The x-monotone polycurve under construction. Its subcurves have the
      // same orientation, and they are either all vertical or all non-vertical.
      X_monotone_curve_2 x_polycurve;
      Comparison_result dir = EQUAL;    // the orientation of its subcurves
      bool is_vertical = false;         // are its subcurves vertical?

      // Output x_polycurve (by moving it) and clear it.
      auto output = [&]() {
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
        if (dir == LARGER) x_polycurve = m_poly_traits.construct_opposite_2_object()(x_polycurve);
#endif
        *oi++ = Make_x_monotone_result(std::move(x_polycurve));
        x_polycurve.clear();
      };

      // Subdivide every subcurve of cv into x-monotone subcurves, and append
      // each one, by moving it, to x_polycurve, unless it starts a new
      // x-monotone polycurve; in this case, output x_polycurve first. The
      // subcurves are processed one at a time, so that only the (few)
      // x-monotone subcurves of a single subcurve are stored at a time.
      std::vector<Make_x_monotone_subresult> x_seg_objects;
      for (auto its = cv.subcurves_begin(); its != cv.subcurves_end(); ++its) {
        x_seg_objects.clear();
        make_seg_x_monotone(*its, std::back_inserter(x_seg_objects));
        for (auto& obj : x_seg_objects) {
          auto* x_seg_p = std::get_if<X_monotone_subcurve_2>(&obj);
          CGAL_assertion(x_seg_p != nullptr);
          const Comparison_result seg_dir = cmp_seg_endpts(*x_seg_p);
          const bool is_seg_vert = is_seg_vertical(*x_seg_p);

          // Test that the subcurve is not degenerate, and that cv is
          // continuous and well-oriented.
          CGAL_precondition_msg(! equal(min_seg_v(*x_seg_p), max_seg_v(*x_seg_p)),
                                "Cannot construct a degenerated subcurve");
          CGAL_precondition_msg(is_first ||
                                equal(last_target, (seg_dir == SMALLER) ? min_seg_v(*x_seg_p) : max_seg_v(*x_seg_p)),
                                "cv must form a continuous and well oriented curve.");
          CGAL_precondition_code
            (
             last_target = (seg_dir == SMALLER) ? max_seg_v(*x_seg_p) : min_seg_v(*x_seg_p);
             is_first = false;
             );

          bool is_new = (x_polycurve.number_of_subcurves() == 0) || (seg_dir != dir) || (is_seg_vert != is_vertical);
          if constexpr (! oblivious) {
            // An x-monotone polycurve does not pass through a point on the
            // boundary; it may only start or end there.
            if (! is_new) {
              auto ps_x = geom_traits->parameter_space_in_x_2_object();
              const auto n = x_polycurve.number_of_subcurves();
              const Arr_curve_end target = (dir == SMALLER) ? ARR_MAX_END : ARR_MIN_END;
              const Arr_curve_end source = (dir == SMALLER) ? ARR_MIN_END : ARR_MAX_END;
              is_new = (ps_x(x_polycurve[n-1], target) != ARR_INTERIOR) || (ps_x(*x_seg_p, source) != ARR_INTERIOR);
            }
          }

          if (is_new) {
            if (x_polycurve.number_of_subcurves() != 0) output();
            dir = seg_dir;
            is_vertical = is_seg_vert;
          }
          x_polycurve.push_back(std::move(*x_seg_p));
        }
      }
      if (x_polycurve.number_of_subcurves() != 0) output();
      return oi;
    }
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
    // Explicitly bring the base operator() into scope, unnecessitating the
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
    void operator()(Curve_2& cv, const Subcurve_2& seg) const { cv.push_back(seg); }
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
    Push_front_2(const Polycurve_traits_2& traits) : Base::Push_front_2(traits) {}

    // Normally, the moment the compiler finds a name, it stops looking. In
    // other words, the compiler first finds the operator() in the current
    // class and stops looking, never finding the one in the base class.
    // Explicitly bring the base operator() into scope, unnecessitating the
    // code below.
    using Base::Push_front_2::operator();

    // /*! appends a subcurve to an existing x-monotone polycurve at the front.
    //  */
    // void operator()(X_monotone_curve_2& xcv,
    //                 const X_monotone_subcurve_2& seg)
    //   const
    // { Base::Push_front_2::operator()(xcv, seg); }

    /* appends a subcurve to an existing polycurve at the front. */
    void operator()(Curve_2& cv, const Subcurve_2& seg) const { cv.push_front(seg); }
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
      auto max_vertex = geom_traits->construct_max_vertex_2_object();
      auto cmp_xy = geom_traits->compare_xy_2_object();
      auto cmp_seg_endpts = geom_traits->compare_endpoints_xy_2_object();
      auto split = geom_traits->split_2_object();
      auto opposite = geom_traits->construct_opposite_2_object();

      // Make sure the split point is not one of the curve endpoints.
      CGAL_precondition_code(auto equal = geom_traits->equal_2_object());
      CGAL_precondition((! equal(m_poly_traits.construct_min_vertex_2_object()(xcv), p)));
      CGAL_precondition((! equal(m_poly_traits.construct_max_vertex_2_object()(xcv), p)));

      CGAL_precondition_msg(xcv.number_of_subcurves() > 0, "Cannot split a polycurve of length zero.");

      const std::size_t n = xcv.number_of_subcurves();
      const Comparison_result dir = cmp_seg_endpts(xcv[0]);

      // Obtains the index of the subcurve whose rank, in left-to-right
      // order, is k.
      auto index = [&](std::size_t k) { return (dir == SMALLER) ? k : n - 1 - k; };

      // Locate the subcurve that contains p by scanning the subcurves from
      // left to right; obtain its rank k and the xy-order of p and its right
      // endpoint. The right endpoint of the rightmost subcurve is not
      // examined; by the precondition, it is not p. The scan replaces a
      // binary search followed by equality tests: the split takes time linear
      // in the number of subcurves anyway, as they are copied, and the split
      // point is typically near the left end of the polycurve (e.g., the
      // surface-sweep splits off the part to the left of the current event),
      // where the scan terminates early. Measured on the polyline benchmarks
      // (aggregated and incremental insertion), k <= 1 in 89%-100% of the
      // splits, unless the polycurves are very long (about 50 subcurves). A
      // binary search, which costs about log(n) comparisons wherever p lies,
      // increased the instruction count by up to 2%, and gained only 0.2% on
      // very long polycurves. A galloping search saves comparisons only on
      // very long polycurves.
      std::size_t k = 0;
      Comparison_result res = SMALLER;
      for (; k + 1 < n; ++k) if ((res = cmp_xy(p, max_vertex(xcv[index(k)]))) != LARGER) break;
      if (k + 1 == n) res = SMALLER;
      const std::size_t i = index(k);

      // Clear the output curves.
      xcv1.clear();
      xcv2.clear();

      // The subcurves xcv[0], ..., xcv[i-1] form the beginning of the
      // resulting polycurve that contains xcv[0], namely `head`; the
      // subcurves xcv[i+1], ..., xcv[n-1] form the end of the other resulting
      // polycurve, namely `tail`. If xcv is directed left-to-right, `head`
      // is the left polycurve (xcv1); otherwise it is the right one (xcv2).
      X_monotone_curve_2& head = (dir == SMALLER) ? xcv1 : xcv2;
      X_monotone_curve_2& tail = (dir == SMALLER) ? xcv2 : xcv1;

      // Reserve storage for the output curves; the reservation is idle when
      // the capacity suffices, and as the curves are empty, it moves nothing
      // otherwise.
      head.reserve(i + 1);
      tail.reserve(n - i);

      for (std::size_t j = 0; j < i; ++j) head.push_back(xcv[j]);
      if (res == EQUAL) {
        // p is the right endpoint of xcv[i], which therefore belongs, as a
        // whole, to the left polycurve (xcv1).
        ((dir == SMALLER) ? head : tail).push_back(xcv[i]);
      }
      else {
        // p lies in the interior of xcv[i]. Split xcv[i] and append the two
        // parts, oriented like xcv, to `head` and `tail`.
        X_monotone_subcurve_2 seg1, seg2;   // the left and right parts
        split(xcv[i], p, seg1, seg2);
        if (cmp_seg_endpts(seg1) != dir) seg1 = opposite(seg1);
        if (cmp_seg_endpts(seg2) != dir) seg2 = opposite(seg2);
        head.push_back(std::move((dir == SMALLER) ? seg1 : seg2));
        tail.push_back(std::move((dir == SMALLER) ? seg2 : seg1));
      }
      for (std::size_t j = i+1; j < n; ++j) tail.push_back(xcv[j]);
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
     * given output iterator. As two subcurves may intersect only once, only a
     * single intersection will be contained in the iterator.
     * Note: If the intersection yields an overlap then it will be oriented
     *       from left-to-right.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1, const X_monotone_curve_2& cv2, OutputIterator oi) const {
      using Intersection_point = std::pair<Point_2, Multiplicity>;
      using Intersection_base_result = std::variant<Intersection_point, X_monotone_subcurve_2>;

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

      const bool consistent = (dir1 == dir2);
#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      CGAL_assertion(consistent);
#endif
      const bool invert_ocv = ((dir1 == LARGER) && (dir2 == LARGER));

      const std::size_t n1 = cv1.number_of_subcurves();
      const std::size_t n2 = cv2.number_of_subcurves();

      Point_2 saved_point;
      // The subcurves of the saved x-monotone curve in the order of the
      // traversal, that is, from left to right.
      std::vector<X_monotone_subcurve_2> saved_subcurves;
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

      // Add a subcurve (by moving it) to the saved x-monotone curve
      auto update_saved_xcv = [&](X_monotone_subcurve_2&& sxcv) {
        xcv_saved = true;
        // We maintain the invariant that if the input curves have opposite
        // directions (! consistent), the overlapping curves are directed
        // left=>right. This, however, is not guaranteed by all traits for the
        // subcurves. (It is guaranteed by some traits, but this is
        // insufficient.)  Therefore, we need to enforce it. That is, we make
        // sure the subcurves are also directed left=>right in this case.
        if (! consistent && (cmp_endpts(sxcv) == LARGER)) saved_subcurves.push_back(ctr_opposite(sxcv));
        else saved_subcurves.push_back(std::move(sxcv));
      };

      // Spit out the saved x-monotone curve. If both input curves are
      // directed right=>left (invert_ocv), so is the saved curve; thus, its
      // subcurves, saved from left to right, are moved in reverse order.
      auto spit_saved_xcv = [&]() {
        if (invert_ocv)
          *oi++ = X_monotone_curve_2(std::make_move_iterator(saved_subcurves.rbegin()),
                                     std::make_move_iterator(saved_subcurves.rend()));
        else
          *oi++ = X_monotone_curve_2(std::make_move_iterator(saved_subcurves.begin()),
                                     std::make_move_iterator(saved_subcurves.end()));
        saved_subcurves.clear();
        xcv_saved = false;
      };

      // Spit out an overlapping subcurve (as an x-monotone curve).
      auto spit_xcv = [&](X_monotone_subcurve_2&& sxcv) {
        CGAL_assertion(! xcv_saved);
        update_saved_xcv(std::move(sxcv));
        spit_saved_xcv();
      };

      // Determines whether an end of the range common to the subcurves xs1
      // and xs2 (whose ranges overlap) lies on both subcurves; that is,
      // whether it is an endpoint of both, or the endpoint of one lies on
      // the other.
      auto is_on_both = [&](const X_monotone_subcurve_2& xs1, const X_monotone_subcurve_2& xs2,
                            Arr_curve_end ce, Comparison_result res) {
        if (res == EQUAL) return true;
        // The end of the common range is the end of the subcurve that
        // starts last (ARR_MIN_END) or ends first (ARR_MAX_END).
        const bool first = (ce == ARR_MIN_END) ? (res == LARGER) : (res == SMALLER);
        return first ? (cmp_y_at_x(xs1, ce, xs2) == EQUAL) : (cmp_y_at_x(xs2, ce, xs1) == EQUAL);
      };

      // Fast path: each curve consists of a single subcurve; thus, there are
      // no joints, at which intersections reported for adjacent subcurves
      // must be merged. An intersection point at an endpoint of a subcurve
      // is reported with multiplicity 0 (undefined), as in the general case
      // below; such a point must be the first intersection, at the left end
      // of the common range, or the last one, at its right end.
      if ((n1 == 1) && (n2 == 1)) {
        const auto& xs1 = cv1[0];
        const auto& xs2 = cv2[0];
        std::vector<Intersection_base_result> xsecs;
        intersect(xs1, xs2, std::back_inserter(xsecs));
        for (auto it = xsecs.begin(); it != xsecs.end(); ++it) {
          if (const auto* xp_p = std::get_if<Intersection_point>(&*it)) {

            auto is_at_end = [&](Arr_curve_end ce)
            { return is_on_both(xs1, xs2, ce, cmp_xy(xs1, ce, xs2, ce)); };

            if (((it == xsecs.begin()) && is_at_end(ARR_MIN_END)) ||
                ((std::next(it) == xsecs.end()) && is_at_end(ARR_MAX_END)))
              *oi++ = std::make_pair(xp_p->first, 0);
            else *oi++ = *xp_p;
          }
          else spit_xcv(std::move(std::get<X_monotone_subcurve_2>(*it)));
        }
        return oi;
      }

      std::size_t i1 = (dir1 == SMALLER) ? 0 : n1-1;
      std::size_t i2 = (dir2 == SMALLER) ? 0 : n2-1;

      auto left_res = cmp_xy(cv1[i1], ARR_MIN_END, cv2[i2], ARR_MIN_END);
      if (left_res == SMALLER) {
        // cv1's left endpoint is to the left of cv2's left endpoint.
        // Locate the index i1 of the subcurve in cv1 which contains cv2's
        // left endpoint.
        i1 = m_poly_traits.locate_impl(cv1, cv2[i2], ARR_MIN_END, All_sides_oblivious_category());
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

      // Traverse the pairs of subcurves cv1[i1] and cv2[i2] whose x-ranges
      // overlap, from left to right. Let [L, R] denote the range common to
      // the subcurves of a pair, where L < R. An intersection of the pair at
      // L or at R is also an intersection of the adjacent pair (unless L or
      // R is an endpoint of cv1 or cv2); it is kept pending, so that it is
      // reported only once, and merged with an adjacent overlap. The
      // following invariant holds: an intersection is pending iff L lies on
      // both subcurves, that is, iff the first intersection of the pair is
      // at L. (Initially, this has been established above; afterwards, L of
      // a pair is R of the preceding pair.) Similarly, the last intersection
      // of the pair is at R iff R lies on both subcurves. Thus, no equality
      // test between a (constructed) intersection point and an endpoint of a
      // subcurve is necessary.
      std::vector<Intersection_base_result> xsecs;
      do {
        auto right_res = cmp_xy(cv1[i1], ARR_MAX_END, cv2[i2], ARR_MAX_END);

        // Determine whether R lies on both subcurves.
        auto is_at_right = [&]() { return is_on_both(cv1[i1], cv2[i2], ARR_MAX_END, right_res); };

        xsecs.clear();
        intersect(cv1[i1], cv2[i2], std::back_inserter(xsecs));

        if (! xsecs.empty()) {
          // It is impossible to have a pending point and a pending
          // x-monotone curve concurrently
          CGAL_assertion(! (xcv_saved && point_saved));
          const bool pending = xcv_saved || point_saved;

          auto it = xsecs.begin();
          auto& xsec = *it++;
          if (it == xsecs.end()) {
            // There is exactly one intersection
            if (const auto* xp_p = std::get_if<Intersection_point>(&xsec)) {
              // The intersection is a point
              if (pending) {
                // The point is L, which is pending.
                if (xcv_saved) spit_saved_xcv();
                else spit_saved_point();
              }
              else if (is_at_right()) update_saved_point(xp_p->first);
              else *oi++ = *xp_p;
            }
            else {
              // The intersection is an x-monotone curve
              auto& sxcv = std::get<X_monotone_subcurve_2>(xsec);
              point_saved = false;
              const bool at_right = is_at_right();
              if (at_right || pending) update_saved_xcv(std::move(sxcv));
              else spit_xcv(std::move(sxcv));
              if (pending && ! at_right) spit_saved_xcv();
            }
          }
          else {
            // There is more than one intersection
            // Handle the first intersection
            if (std::get_if<Intersection_point>(&xsec)) {
              if (pending) {
                // The point is L, which is pending.
                if (xcv_saved) spit_saved_xcv();
                else spit_saved_point();
              }
              else *oi++ = std::get<Intersection_point>(xsec);
            }
            else {
              auto& sxcv = std::get<X_monotone_subcurve_2>(xsec);
              point_saved = false;
              if (pending) {
                update_saved_xcv(std::move(sxcv));
                spit_saved_xcv();
              }
              else spit_xcv(std::move(sxcv));
            }

            // Handle all but the last intersection
            for (; it != std::prev(xsecs.end()); ++it) {
              if (const auto* xp_p = std::get_if<Intersection_point>(&*it)) *oi++ = *xp_p;
              else spit_xcv(std::move(std::get<X_monotone_subcurve_2>(*it)));
            }

            // Handle the last intersection
            const bool at_right = is_at_right();
            if (const auto* xp_p = std::get_if<Intersection_point>(&*it)) {
              if (at_right) update_saved_point(xp_p->first);
              else *oi++ = *xp_p;
            }
            else {
              auto& sxcv = std::get<X_monotone_subcurve_2>(*it);
              if (at_right) update_saved_xcv(std::move(sxcv));
              else spit_xcv(std::move(sxcv));
            }
          }
        }

        // Advance the indices
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
    bool operator()(const X_monotone_curve_2& cv1, const X_monotone_curve_2& cv2) const {
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
                ((equal(max_vertex(cv1),min_vertex(cv2))) || (equal(max_vertex(cv2),min_vertex(cv1))))) ||
               (// Both are directed from right-to-left
                (dir1 == LARGER) &&
                ((equal(min_vertex(cv1),max_vertex(cv2))) || (equal(max_vertex(cv1),min_vertex(cv2)))))) &&
              (// Either both should be vertical or both should
               // be NOT vertical.
               (ver1 && ver2) || (!ver1 && ! ver2)));
    }
  };

  /*! obtains an `Are_mergeable_2` functor object. */
  Are_mergeable_2 are_mergeable_2_object() const { return Are_mergeable_2(*this); }

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
    void operator()(const X_monotone_curve_2& cv1, const X_monotone_curve_2& cv2, X_monotone_curve_2& c) const {
      CGAL_precondition(m_poly_traits.are_mergeable_2_object()(cv1, cv2));

      auto min_vertex = m_poly_traits.construct_min_vertex_2_object();
      auto max_vertex = m_poly_traits.construct_max_vertex_2_object();
      auto cmp_endpts = m_poly_traits.compare_endpoints_xy_2_object();
      Equal_2 equal = m_poly_traits.equal_2_object();

      c.clear();
      if (// Either both are left-to-right and cv2 is to the right of cv1
          ((cmp_endpts(cv1)==SMALLER) && (equal(max_vertex(cv1), min_vertex(cv2)))) ||
          // or both are right-to-left and cv2 is to the left of cv1
          ((cmp_endpts(cv1)==LARGER) && (equal(min_vertex(cv1), max_vertex(cv2))))) {
        const std::size_t n1 = cv1.number_of_subcurves();
        const std::size_t n2 = cv2.number_of_subcurves();
        std::size_t i;
        c.reserve(n1 + n2);

        // cv2 extends cv1 to the right:
        for (i = 0; i < n1 - 1; ++i) c.push_back(cv1[i]);

        // Try to merge the to contiguous line subcurves:
        if (m_poly_traits.subcurve_traits_2()->are_mergeable_2_object()(cv1[n1 - 1], cv2[0])) {
          X_monotone_subcurve_2 seg;
          m_poly_traits.subcurve_traits_2()->
            merge_2_object()(cv1[n1 - 1], cv2[0], seg);
          c.push_back(std::move(seg));
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
    Construct_curve_2(const Polycurve_traits_2& traits) : m_poly_traits(traits) {}

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
    Curve_2 constructor_impl(ForwardIterator /* begin */, ForwardIterator /* end */, std::true_type) const
    {  CGAL_error_msg("Cannot construct a polycurve from a range of points!"); }

    /*! Construction implementation from a range of subcurves.
     *  Note that the subcurves in the range are NOT necessarily x-monotone,
     *  thus it is impossible to test (even in precondition) whether the input
     *  forms a continuous and well oriented polycurve.
     *  \pre Range should contain at least one subcurve.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end, std::false_type) const {
      // Range has to contain at least one subcurve
      CGAL_precondition(begin != end);
      return Curve_2(begin, end);
    }
  };

  /*! obtains a `Construct_curve_2` functor object. */
  Construct_curve_2 construct_curve_2_object() const { return Construct_curve_2(*this); }
};

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
