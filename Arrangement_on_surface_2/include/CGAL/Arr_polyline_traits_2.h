// Copyright (c) 2000-2014
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Efi Fogel       <efifogel@gmail.com>
//             Simon Giraudot  <simon.giraudot@geometryfactory.com>
//             Waqar Khan      <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_ARR_POLYLINE_TRAITS_2_H
#define CGAL_ARR_POLYLINE_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The traits-class for the linear piecewise (polyline) type of curves of the
 * arrangement package.
 */

#include <iterator>
#include <type_traits>
#include <vector>

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polycurve_traits_2.h>

namespace CGAL {

// If no template instantiation is provided, it will be instantiated
// with Arr_segment_traits.
template <typename SegmentTraits_2 = Arr_segment_traits_2<> >
class Arr_polyline_traits_2 : public Arr_polycurve_traits_2<SegmentTraits_2> {
public:
  using Segment_traits_2 = SegmentTraits_2;

  // For completeness
  using Segment_2 = typename Segment_traits_2::Curve_2;
  using X_monotone_segment_2 = typename Segment_traits_2::X_monotone_curve_2;

private:
  using Self = Arr_polyline_traits_2<Segment_traits_2>;
  using Base = Arr_polycurve_traits_2<Segment_traits_2>;

public:
  /// \name Types inherited from the polycurve traits class.
  //@{

  // Tag definitions:
  using Has_left_category = typename Base::Has_left_category;
  using Has_merge_category = typename Base::Has_merge_category;

  using Left_side_category = typename Base::Left_side_category;
  using Bottom_side_category = typename Base::Bottom_side_category;
  using Top_side_category = typename Base::Top_side_category;
  using Right_side_category = typename Base::Right_side_category;

  using All_sides_oblivious_category = typename Base::All_sides_oblivious_category;

  using X_monotone_subcurve_2 = typename Base::X_monotone_subcurve_2;
  using Subcurve_2 = typename Base::Subcurve_2;

  using Point_2 = typename Base::Point_2;
  using Curve_2 = typename Base::Curve_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  using Multiplicity = typename Base::Multiplicity;

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
  using Parameter_space_in_x_2 = typename Base::Parameter_space_in_x_2;
  using Parameter_space_in_y_2 = typename Base::Parameter_space_in_y_2;
  using Compare_x_on_boundary_2 = typename Base::Compare_x_on_boundary_2;
  using Compare_x_near_boundary_2 = typename Base::Compare_x_near_boundary_2;
  using Compare_y_on_boundary_2 = typename Base::Compare_y_on_boundary_2;
  using Compare_y_near_boundary_2 = typename Base::Compare_y_near_boundary_2;
  using Is_on_y_identification_2 = typename Base::Is_on_y_identification_2;
  using Is_on_x_identification_2 = typename Base::Is_on_x_identification_2;

  using Split_2 = typename Base::Split_2;
  using Intersect_2 = typename Base::Intersect_2;
  using Are_mergeable_2 = typename Base::Are_mergeable_2;
  using Merge_2 = typename Base::Merge_2;

  using Number_of_points_2 = typename Base::Number_of_points_2;
  using Trim_2 = typename Base::Trim_2;

  //@}

public:
  /*! constructs default. */
  Arr_polyline_traits_2() : Base() {}

  /*! constructs from a segment traits
   * \param geom_traits an already existing segment traits which is passed will
   *        be used by the class.
   */
  Arr_polyline_traits_2(const Segment_traits_2* geom_traits) : Base(geom_traits) {}

  /// \name Types and functors defined here.

  /* Functor to augment a polyline by either adding a vertex or a segment
   * at the back.
   */
  class Push_back_2 : public Base::Push_back_2 {
  protected:
    using Polyline_traits_2 = Arr_polyline_traits_2<Segment_traits_2>;

  public:
    /*! constructs. */
    Push_back_2(const Polyline_traits_2& traits) : Base::Push_back_2(traits) {}

    // Normally, the moment the compiler finds a name, it stops looking. In
    // other words, the compiler first finds the operator() in the current
    // class and stops looking, never finding the one in the base class.
    // Explicitly bring the base operator() into scope, unnecessitating the
    // code below.
    using Base::Push_back_2::operator();

    // /*! appends a segment `seg` to an existing polyline `cv` at the back. */
    // void operator()(Curve_2& cv, const Segment_2& seg) const
    // { Base::Push_back_2::operator() (cv, seg); }

    // /* appends a segment `seg` to an existing polyline `xcv` at the back. */
    // void operator()(X_monotone_curve_2& xcv,
    //                 const X_monotone_subcurve_2& seg) const
    // { Base::Push_back_2::operator()(xcv, seg); }

    /* appends a point `p` to an existing polyline `cv` at the back. */
    void operator()(Curve_2& cv, const Point_2& p) const {
      typedef typename Curve_2::size_type size_type;
      size_type num_seg = cv.number_of_subcurves();
      CGAL_precondition(num_seg > 0);
      int last_seg = num_seg-1;

      const Segment_traits_2* seg_traits = this->m_poly_traits.subcurve_traits_2();
      auto cmp_seg_endpts = seg_traits->compare_endpoints_xy_2_object();
      auto ctr = seg_traits->construct_curve_2_object();

      /* Since it is desired to maintain `cv` well-oriented, we have
       * to append the segment [cv[last_seg].target(),p]. The
       * following test determines which end of the last segment is
       * the target, i.e. the actual end of `cv`.
       */
      if (cmp_seg_endpts(cv[last_seg]) == SMALLER) {
        auto get_max_v = seg_traits->construct_max_vertex_2_object();
        cv.push_back(ctr(get_max_v(cv[last_seg]), p));
      }
      else {
        auto get_min_v = seg_traits->construct_min_vertex_2_object();
        cv.push_back(ctr(get_min_v(cv[last_seg]), p));
      }
    }

    /* appends a point `p` to an existing polyline `xcv` at the back. */
    void operator()(X_monotone_curve_2& xcv, const Point_2& p) const {
      typedef typename X_monotone_curve_2::size_type size_type;
      size_type num_seg = xcv.number_of_subcurves();
      CGAL_precondition(num_seg > 0);

      const Segment_traits_2* seg_traits = this->m_poly_traits.subcurve_traits_2();
      CGAL_precondition_code
        (
         auto comp_x = seg_traits->compare_x_2_object();
         auto comp_xy = seg_traits->compare_xy_2_object();
         auto is_vertical = this->m_poly_traits.is_vertical_2_object();
         );

      if (seg_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) {
        // xcv is oriented left-to-right
        auto get_max_v = seg_traits->construct_max_vertex_2_object();
        CGAL_precondition
          (
           (! is_vertical(xcv) && (comp_x(get_max_v(xcv[num_seg-1]), p) == SMALLER)) ||
           (is_vertical(xcv) &&
            (comp_x(get_max_v(xcv[num_seg-1]), p) == EQUAL) &&
            (comp_xy(get_max_v(xcv[num_seg-1]), p) == SMALLER))
           );
        xcv.push_back(X_monotone_subcurve_2(get_max_v(xcv[num_seg-1]), p));
      }
      else {
        // xcv is oriented right-to-left
        typename Segment_traits_2::Construct_min_vertex_2 get_min_v = seg_traits->construct_min_vertex_2_object();
        CGAL_precondition
          (
           (! is_vertical(xcv) && (comp_x(get_min_v(xcv[num_seg-1]), p) == LARGER)) ||
           (is_vertical(xcv) &&
            (comp_x(get_min_v(xcv[num_seg-1]), p) == EQUAL) &&
            (comp_xy(get_min_v(xcv[num_seg-1]), p) == LARGER))
           );
        xcv.push_back(X_monotone_subcurve_2(get_min_v(xcv[num_seg-1]), p));
      }
    }
  };

  /*! obtains a Push_back_2 functor object. */
  Push_back_2 push_back_2_object() const { return Push_back_2(*this); }

  /* Functor to augment a polyline by either adding a vertex or a segment
   * at the front.
   * TODO: Test all the operator()'s. (Don't forget vertical cases!)
   */
  class Push_front_2 : public Base::Push_front_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;

  public:
    /*! constructs. */
    Push_front_2(const Polyline_traits_2& traits) : Base::Push_front_2(traits) {}

    // Normally, the moment the compiler finds a name, it stops looking. In
    // other words, the compiler first finds the operator() in the current
    // class and stops looking, never finding the one in the base class.
    // Explicitly bring the base operator() into scope, unnecessitating the
    // code below.
    using Base::Push_front_2::operator();

    // /*! appends a segment `seg` to an existing polyline `cv` at the front. */
    // void operator()(Curve_2& cv, const Subcurve_2& seg) const
    // { Base::Push_front_2::operator()(cv, seg); }

    // /*! appends a segment `seg` to an existing polyline `xcv` at the front. */
    // void operator()(X_monotone_curve_2& xcv,
    //                 const X_monotone_subcurve_2& seg) const
    // { Base::Push_front_2::operator()(xcv, seg); }

    /* Append a point `p` to an existing polyline `cv` at the front. */
    void operator()(Curve_2& cv, const Point_2& p) const {
      CGAL_precondition_code
        (
         typedef typename Curve_2::size_type size_type;
         size_type num_seg = cv.number_of_subcurves();
         );
      CGAL_precondition(num_seg > 0);

      const Segment_traits_2* geom_traits = this->m_poly_traits.subcurve_traits_2();
      auto cmp_seg_endpts = geom_traits->compare_endpoints_xy_2_object();
      auto ctr = geom_traits->construct_curve_2_object();

      if (cmp_seg_endpts(cv[0]) == SMALLER) {
        auto get_min_v = geom_traits->construct_min_vertex_2_object();
        cv.push_front(ctr(p, get_min_v(cv[0])));
      }
      else {
        auto get_max_v = geom_traits->construct_max_vertex_2_object();
        cv.push_front(ctr(p, get_max_v(cv[0])));
      }
    }

    /*! appends a point `p` to an existing polyline `xcv` at the front. */
    void operator()(X_monotone_curve_2& xcv, const Point_2& p) const {
      const auto* geom_traits = this->m_poly_traits.subcurve_traits_2();
      CGAL_precondition_code
        (
         typedef typename X_monotone_curve_2::size_type size_type;
         size_type num_seg = xcv.number_of_subcurves();
         auto comp_x = geom_traits->compare_x_2_object();
         auto comp_xy = geom_traits->compare_xy_2_object();
         auto is_vertical = this->m_poly_traits.is_vertical_2_object();
         );
      CGAL_precondition(num_seg > 0);

      if (geom_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) {
        // xcv is oriented left-to-right
        auto get_max_v = geom_traits->construct_max_vertex_2_object();
        CGAL_precondition
          (
           (! is_vertical(xcv) && (comp_x(get_max_v(xcv[0]), p) == LARGER)) ||
           (is_vertical(xcv) &&
            (comp_x(get_max_v(xcv[0]), p) == EQUAL) &&
            (comp_xy(get_max_v(xcv[0]), p) == LARGER))
           );
        xcv.push_front(X_monotone_subcurve_2(p, get_max_v(xcv[0])));
      }
      else {
        // xcv is oriented right-to-left
        auto get_min_v = geom_traits->construct_min_vertex_2_object();
        CGAL_precondition
          (
           (! is_vertical(xcv) && (comp_x(get_min_v(xcv[0]), p) == SMALLER)) ||
           (is_vertical(xcv) &&
            (comp_x(get_min_v(xcv[0]), p) == EQUAL) &&
            (comp_xy(get_min_v(xcv[0]), p) == SMALLER))
           );
        xcv.push_front(X_monotone_subcurve_2(p, get_min_v(xcv[0])));
      }
    }
  };

  /*! obtains a Push_front_2 functor object. */
  Push_front_2 push_front_2_object() const { return Push_front_2(*this); }

  /*! constructs a general curve from :
   * 1. two points,
   * 2. a sub curve,
   * 3. a range of points, or
   * 4. a range of subcurves.
   */
  class Construct_curve_2 : public Base::Construct_curve_2 {
  protected:
    using Polyline_traits_2 = Arr_polyline_traits_2<Segment_traits_2>;

  public:
    /*! constructs.
     */
    Construct_curve_2(const Polyline_traits_2& traits) : Base::Construct_curve_2(traits) {}

    /* Obtain an polyline connecting two given endpoints.
     */
    Curve_2 operator()(const Point_2& p, const Point_2& q) const {
      const Segment_traits_2* seg_traits = this->m_poly_traits.subcurve_traits_2();
      return seg_traits->construct_curve_2_object()(p, q);
    }

    /* Obtain a polyline consists of one given segment.
     */
    Curve_2 operator()(const Subcurve_2& seg) const { return Base::Construct_curve_2::operator()(seg); }

    /* Construct a well-oriented polyline from a range of either
     * `Segment_traits::Point_2` or `Segment_traits::Subcurve_2`.
     */
    template <typename ForwardIterator>
    Curve_2 operator()(ForwardIterator begin, ForwardIterator end) const {
      typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
      typedef typename std::is_same<VT, Point_2>::type Is_point;
      // Dispatch the range to the appropriate implementation.
      return constructor_impl(begin, end, Is_point());
    }

    /* Construct a well-oriented polyline from a range of either
     * `Segment_traits::Point_2` or `Segment_traits::Subcurve_2`.
     */
    template <typename Range>
    Curve_2 operator()(const Range& range) const { return operator()(range.begin(), range.end()); }

    /*! Construction implementation from a range of segments.
     *  Note that the segments in the range are NOT necessarily x-monotone,
     *  thus it is impossible to test (even in precondition) whether the input
     *  forms a continuous and well oriented polyline.
     *  \pre Range should contain at least one segment.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end, std::false_type) const
    { return Base::Construct_curve_2::operator()(begin, end); }

    /*! constructs a polyline from a range of points.
     * \pre The range contains at least two points
     * \pre Consecutive points are disjoint.
     * \return Well-oriented polyline connecting the given points. The order
     *         of the vertices is determined by their order in the range.
     *         Furthermore, the orientation of the polyline is induced by
     *         their order.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end, std::true_type) const {
      // The range must not contain a single point.
      CGAL_precondition_msg(std::distance(begin, end) != 1,
                            "Range of points must not contain a single point");
      // The range must contain at least one segment.
      CGAL_precondition(begin != end);

      CGAL_precondition_code(auto equal = this->m_poly_traits.equal_2_object());
      auto ctr_seg = this->m_poly_traits.subcurve_traits_2()->construct_curve_2_object();

      // Construct the segments directly in the polyline, which is reserved
      // in advance.
      Curve_2 cv;
      if (begin == end) return cv;
      using Category = typename std::iterator_traits<ForwardIterator>::iterator_category;
      if constexpr (std::is_base_of_v<std::forward_iterator_tag, Category>)
        cv.reserve(std::distance(begin, end) - 1);
      for (auto curr = begin, next = std::next(begin); next != end; curr = next++) {
        CGAL_precondition_msg(! equal(*curr, *next), "Cannot construct a degenerated segment");
        cv.push_back(ctr_seg(*curr, *next));
      }
      return cv;
    }
  };

  /*! obtains a Construct_curve_2 functor object. */
  Construct_curve_2 construct_curve_2_object() const { return Construct_curve_2(*this); }

  /*! constructs an x-monotone curve. */
  class Construct_x_monotone_curve_2 :
    public Base::Construct_x_monotone_curve_2 {
  protected:
    using Polyline_traits_2 = Arr_polyline_traits_2<Segment_traits_2>;

  public:
    /*! constructs.
     */
    Construct_x_monotone_curve_2(const Polyline_traits_2& traits) :
      Base::Construct_x_monotone_curve_2(traits)
    {}

    /*! obtains an x-monotone polyline connecting two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q must not be the same.
     * \return A segment connecting p and q.
     */
    X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q) const {
      CGAL_precondition_msg
        (!this->m_poly_traits.subcurve_traits_2()->equal_2_object()(p,q),
         "Cannot construct a degenerated segment as a polyline");
      X_monotone_subcurve_2 seg = this->m_poly_traits.subcurve_traits_2()->construct_x_monotone_curve_2_object()(p, q);

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      if (this->m_poly_traits.subcurve_traits_2()->compare_xy_2_object()(p,q) == LARGER)
        seg = this->m_poly_traits.subcurve_traits_2()->construct_opposite_2_object()(seg);
#endif

      return X_monotone_curve_2(seg);
    }

    /*! obtains an x-monotone polyline that consists of one given segment.
     * \param seg input segment.
     * \pre seg is not degenerated.
     * \return An x-monotone polyline with one segment.
     */
    X_monotone_curve_2 operator()(const X_monotone_subcurve_2& seg) const
    { return Base::Construct_x_monotone_curve_2::operator()(seg); }

    /*! constructs an x-monotone polyline from a range of elements.
     * \pre Range should from a continuous well-oriented x-monotone polyline.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 operator()(ForwardIterator begin, ForwardIterator end) const {
      typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
      typedef typename std::is_same<VT, Point_2>::type Is_point;
      // Dispatch the range to the appropriate implementation.
      return constructor_impl(begin, end, Is_point());
    }

    /*! Construction implementation from a range of segments.
     * \param begin An iterator pointing to the first segment in the range.
     * \param end An iterator pointing to the past-the-end segment
     * in the range.
     * \pre The range contains at least one segment.
     * \pre Segments correspond to a well-oriented polyline. That
     *      is, the target of the i-th segment is an source of the
     *      (i+1)th segment.
     * \pre The sequence of segments in the range forms a weak x-monotone
     *      polyline.
     * \pre The container should support bidirectional iteration.
     * \return A continuous, well-oriented x-monotone polyline which
     *         is directed either left-to-right or right-to-left
     *         depending on the segments in the input.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end, std::false_type) const
    { return Base::Construct_x_monotone_curve_2::operator()(begin, end); }

    /*! Construction of an x-monotone polyline from a range of points.
     * The polyline may be oriented left-to-right or right-to-left
     * depending on the lexicographical order of the points in the
     * input.
     * \pre Range contains at least two points.
     * \pre No two consecutive points are the same.
     * \pre The points form an continuous well-oriented x-monotone polyline.
     * \post By the construction the returned polyline is well-oriented.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end, std::true_type) const {
      // The range must not contain a single point.
      CGAL_precondition_msg(std::distance(begin, end) != 1,
                            "Range of points must not contain a single point");

      // Container of the segments to be constructed from the range of points
      std::vector<X_monotone_subcurve_2> segs;
      using Category = typename std::iterator_traits<ForwardIterator>::iterator_category;
      if constexpr (std::is_base_of_v<std::forward_iterator_tag, Category>)
        segs.reserve(std::distance(begin, end) - 1);
      // Make sure the range of points contains at least two points.
      ForwardIterator next = begin;
      ForwardIterator curr = next++;

      CGAL_precondition_code
        (
         const auto* geom_traits = this->m_poly_traits.subcurve_traits_2();
         // Initialize two comparison functors
         auto compare_x = geom_traits->compare_x_2_object();
         auto compare_xy = geom_traits->compare_xy_2_object();
         // Make sure there is no changed of directions.
         // Saves the comp_x between the first two points
         const Comparison_result cmp_x_res = compare_x(*curr, *next);
         // Save the comp_xy between the first two points
         const Comparison_result cmp_xy_res = compare_xy(*curr, *next);
         );

      // Assure that the first two points are not the same.
      // Note that this also assures that non of the consecutive
      // points are equal in the whole range.
      CGAL_precondition(cmp_xy_res != EQUAL);

      while (next != end) {
        CGAL_precondition(compare_xy(*curr, *next) == cmp_xy_res);
        CGAL_precondition(compare_x(*curr, *next) == cmp_x_res);

        segs.emplace_back(*curr, *next);
        curr = next++;
      }

      // Move the segments into the polyline; the construction from a range of
      // segments also handles CGAL_ALWAYS_LEFT_TO_RIGHT.
      return operator()(std::make_move_iterator(segs.begin()), std::make_move_iterator(segs.end()));
    }
  };

  /*! obtains a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(*this); }

  //
  class Approximate_2 : public Base::Approximate_2 {
  protected:
    using Traits = Arr_polyline_traits_2<Segment_traits_2>;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! constructs
     * \param traits the traits.
     */
    Approximate_2(const Traits& traits) :
      Base::Approximate_2(*(traits.subcurve_traits_2())),
      m_traits(traits)
    {}

    friend class Arr_polyline_traits_2<Segment_traits_2>;

  public:
    using Approximate_number_type = typename Base::Approximate_number_type;
    using Approximate_point_2 = typename Base::Approximate_point_2;
    using Approximate_kernel = typename Base::Approximate_kernel;

    /*! obtains an approximation of a point coordinate.
     * \param p the exact point.
     * \param i the coordinate index (either 0 or 1).
     * \pre `i` is either 0 or 1.
     * \return An approximation of `p`'s \f$x\f$-coordinate (if `i` == 0), or an
     *         approximation of `p`'s \f$y\f$-coordinate (if `i` == 1).
     */
    Approximate_number_type operator()(const Point_2& p, int i) const { return Base::Approximate_2::operator()(p, i); }

    Approximate_point_2 operator()(const Point_2& p) const { return Base::Approximate_2::operator()(p); }

    /*! obtains an approximation of an \f$x\f$-monotone curve. */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv, double /* error */,
                              OutputIterator oi, bool l2r = true) const {
      if (l2r) {
        for (auto it = xcv.points_begin(); it != xcv.points_end(); ++it) {
          const auto& p = *it;
          auto x = CGAL::to_double(p.x());
          auto y = CGAL::to_double(p.y());
          *oi++ = Approximate_point_2(x, y);
        }
        return oi;
      }

      for (auto it = xcv.points_rbegin(); it != xcv.points_rend(); ++it) {
        const auto& p = *it;
        auto x = CGAL::to_double(p.x());
        auto y = CGAL::to_double(p.y());
        *oi++ = Approximate_point_2(x, y);
      }
      return oi;
    }
  };

  /*! obtains an Approximate_2 functor object. */
  Approximate_2 approximate_2_object() const { return Approximate_2(*this); }

  /*! obtains the segment traits.
   * \return the segment traits.
   */
  CGAL_DEPRECATED const Segment_traits_2* segment_traits_2() const { return this->subcurve_traits_2(); }
};

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
