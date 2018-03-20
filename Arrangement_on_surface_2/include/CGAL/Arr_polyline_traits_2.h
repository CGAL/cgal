// Copyright (c) 2000-2014
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// Author(s)     : Waqar Khan <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_ARR_POLYLINE_TRAITS_2_H
#define CGAL_ARR_POLYLINE_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The traits-class for the linear piece-wiese(polyline) type of curves of the
 * arrangement package.
 */

#include <list>
#include <iterator>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polycurve_traits_2.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>

namespace CGAL {

// If no template instantiation is provided, it will be instantiated
// with Arr_segment_traits.
template <typename SegmentTraits_2 = Arr_segment_traits_2<> >
class Arr_polyline_traits_2 : public Arr_polycurve_traits_2<SegmentTraits_2> {
public:
  typedef SegmentTraits_2                             Segment_traits_2;

  // For completeness
  typedef typename Segment_traits_2::Curve_2            Segment_2;
  typedef typename Segment_traits_2::X_monotone_curve_2 X_monotone_segment_2;

private:
  typedef Arr_polyline_traits_2<Segment_traits_2>     Self;
  typedef Arr_polycurve_traits_2<Segment_traits_2>    Base;

public:
  /// \name Types inherited from the polycurve traits class.
  //@{

  // Tag definitions:
  typedef typename Base::Has_left_category            Has_left_category;
  typedef typename Base::Has_merge_category           Has_merge_category;
  typedef typename Base::Has_do_intersect_category
    Has_do_intersect_category;

  typedef typename Base::Left_side_category           Left_side_category;
  typedef typename Base::Bottom_side_category         Bottom_side_category;
  typedef typename Base::Top_side_category            Top_side_category;
  typedef typename Base::Right_side_category          Right_side_category;

  typedef typename Base::Are_all_sides_oblivious_tag
    Are_all_sides_oblivious_tag;

  typedef typename Base::X_monotone_subcurve_2        X_monotone_subcurve_2;
  typedef typename Base::Subcurve_2                   Subcurve_2;

  typedef typename Base::Point_2                      Point_2;
  typedef typename Base::Curve_2                      Curve_2;
  typedef typename Base::X_monotone_curve_2           X_monotone_curve_2;

  typedef typename Base::Multiplicity                 Multiplicity;

  typedef typename Base::Compare_x_2                  Compare_x_2;
  typedef typename Base::Compare_xy_2                 Compare_xy_2;
  typedef typename Base::Construct_min_vertex_2       Construct_min_vertex_2;
  typedef typename Base::Construct_max_vertex_2       Construct_max_vertex_2;
  typedef typename Base::Is_vertical_2                Is_vertical_2;
  typedef typename Base::Compare_y_at_x_2             Compare_y_at_x_2;
  typedef typename Base::Compare_y_at_x_left_2        Compare_y_at_x_left_2;
  typedef typename Base::Compare_y_at_x_right_2       Compare_y_at_x_right_2;
  typedef typename Base::Equal_2                      Equal_2;
  typedef typename Base::Compare_endpoints_xy_2       Compare_endpoints_xy_2;
  typedef typename Base::Construct_opposite_2         Construct_opposite_2;
  typedef typename Base::Approximate_2                Approximate_2;
  typedef typename Base::Parameter_space_in_x_2       Parameter_space_in_x_2;
  typedef typename Base::Parameter_space_in_y_2       Parameter_space_in_y_2;
  typedef typename Base::Compare_x_on_boundary_2      Compare_x_on_boundary_2;
  typedef typename Base::Compare_x_at_limit_2         Compare_x_at_limit_2;
  typedef typename Base::Compare_x_near_limit_2       Compare_x_near_limit_2;
  typedef typename Base::Compare_y_on_boundary_2      Compare_y_on_boundary_2;
  typedef typename Base::Compare_y_near_boundary_2    Compare_y_near_boundary_2;
  typedef typename Base::Is_on_y_identification_2     Is_on_y_identification_2;
  typedef typename Base::Is_on_x_identification_2     Is_on_x_identification_2;

  typedef typename Base::Split_2                      Split_2;
  typedef typename Base::Intersect_2                  Intersect_2;
  typedef typename Base::Are_mergeable_2              Are_mergeable_2;
  typedef typename Base::Merge_2                      Merge_2;

  typedef typename Base::Number_of_points_2           Number_of_points_2;
  typedef typename Base::Trim_2                       Trim_2;

  //@}

public:
  /*! Construct default. */
  Arr_polyline_traits_2() : Base() {}

  /*! Construct from a segment traits
   * \param geom_traits an already existing segment tarits which is passed will
   *        be used by the class.
   */
  Arr_polyline_traits_2(const Segment_traits_2* geom_traits) :
    Base(geom_traits)
  {}

  /// \name Types and functors defined here.
  //@{

  /* Functor to augment a polyline by either adding a vertex or a segment
   * at the back.
   */
  class Push_back_2 : public Base::Push_back_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;

  public:
    /*! Constructor. */
    Push_back_2(const Polyline_traits_2& traits) :
      Base::Push_back_2(traits)
    {}

    // Normally, the moment the compiler finds a name, it stops looking. In
    // other words, the compiler first finds the operator() in the current
    // class and stops looking, never finding the one in the base class.
    // Explicitly bring the base operator() into scope, unnecesitating the
    // code below.
    using Base::Push_back_2::operator();

    // /*! Append a segment `seg` to an existing polyline `cv` at the back. */
    // void operator()(Curve_2& cv, const Segment_2& seg) const
    // { Base::Push_back_2::operator() (cv, seg); }

    // /* Append a segment `seg` to an existing polyline `xcv` at the back. */
    // void operator()(X_monotone_curve_2& xcv,
    //                 const X_monotone_subcurve_2& seg) const
    // { Base::Push_back_2::operator()(xcv, seg); }

    /* Append a point `p` to an existing polyline `cv` at the back. */
    void operator()(Curve_2& cv, const Point_2& p) const
    {
      typedef typename Curve_2::size_type size_type;
      size_type num_seg = cv.number_of_subcurves();
      CGAL_precondition(num_seg > 0);
      int last_seg = num_seg-1;

      const Segment_traits_2* seg_traits =
        this->m_poly_traits.subcurve_traits_2();
      typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
        seg_traits->compare_endpoints_xy_2_object();

      /* Since it is desired to maintain `cv` well-oriented, we have
       * to append the segment [cv[last_seg].target(),p]. The
       * following test determines which end of the last segment is
       * the target, i.e. the actual end of `cv`.
       */
      if (cmp_seg_endpts(cv[last_seg]) == SMALLER) {
        typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
          seg_traits->construct_max_vertex_2_object();
        cv.push_back(Subcurve_2(get_max_v(cv[last_seg]), p));
      }
      else {
        typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
          seg_traits->construct_min_vertex_2_object();
        cv.push_back(Subcurve_2(get_min_v(cv[last_seg]), p));
      }
    }

    /* Append a point `p` to an existing polyline `xcv` at the back. */
    void operator()(X_monotone_curve_2& xcv, const Point_2& p) const
    {
      typedef typename X_monotone_curve_2::size_type size_type;
      size_type num_seg = xcv.number_of_subcurves();
      CGAL_precondition(num_seg > 0);

      const Segment_traits_2* seg_traits =
        this->m_poly_traits.subcurves_traits_2();
      CGAL_precondition_code
        (
         typename Segment_traits_2::Compare_x_2 comp_x =
           seg_traits->compare_x_2_object();
         typename Segment_traits_2::Compare_xy_2 comp_xy =
           seg_traits->compare_xy_2_object();
         typename Base::Is_vertical_2 is_vertical =
           this->m_poly_traits.is_vertical_2_object();
         );

      if (seg_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) {
        // xcv is oriented left-to-right
        typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
          seg_traits->construct_max_vertex_2_object();
        CGAL_precondition
          (
           (!is_vertical(xcv) &&
            (comp_x(get_max_v(xcv[num_seg-1]), p) == SMALLER)) ||
           (is_vertical(xcv) &&
            (comp_x(get_max_v(xcv[num_seg-1]), p) == EQUAL) &&
            (comp_xy(get_max_v(xcv[num_seg-1]), p) == SMALLER))
           );
        xcv.push_back(X_monotone_subcurve_2(get_max_v(xcv[num_seg-1]), p));
      }
      else {
        // xcv is oriented right-to-left
        typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
          seg_traits->construct_min_vertex_2_object();
        CGAL_precondition
          (
           (!is_vertical(xcv) &&
            (comp_x(get_min_v(xcv[num_seg-1]), p) == LARGER)) ||
           (is_vertical(xcv) &&
            (comp_x(get_min_v(xcv[num_seg-1]), p) == EQUAL) &&
            (comp_xy(get_min_v(xcv[num_seg-1]), p) == LARGER))
           );
        xcv.push_back(X_monotone_subcurve_2(get_min_v(xcv[num_seg-1]), p));
      }
    }
  };

  /*! Obtain a Push_back_2 functor object. */
  Push_back_2 push_back_2_object() const { return Push_back_2(*this); }

  /* Functor to augment a polyline by either adding a vertex or a segment
   * at the front.
   * TODO: Test all the operator()'s. (Don't forget vertical cases!)
   */
  class Push_front_2 : public Base::Push_front_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;

  public:
    /*! Constructor. */
    Push_front_2(const Polyline_traits_2& traits) :
      Base::Push_front_2(traits)
    {}

    // Normally, the moment the compiler finds a name, it stops looking. In
    // other words, the compiler first finds the operator() in the current
    // class and stops looking, never finding the one in the base class.
    // Explicitly bring the base operator() into scope, unnecesitating the
    // code below.
    using Base::Push_front_2::operator();

    // /*! Append a segment `seg` to an existing polyline `cv` at the front. */
    // void operator()(Curve_2& cv, const Subcurve_2& seg) const
    // { Base::Push_front_2::operator()(cv, seg); }

    // /*! Append a segment `seg` to an existing polyline `xcv` at the front. */
    // void operator()(X_monotone_curve_2& xcv,
    //                 const X_monotone_subcurve_2& seg) const
    // { Base::Push_front_2::operator()(xcv, seg); }

    /* Append a point `p` to an existing polyline `cv` at the front. */
    void operator()(Curve_2& cv, const Point_2& p) const
    {
      CGAL_precondition_code
        (
         typedef typename Curve_2::size_type size_type;
         size_type num_seg = cv.number_of_subcurves();
         );
      CGAL_precondition(num_seg > 0);

      const Segment_traits_2* geom_traits =
        this->m_poly_traits.subcurve_traits_2();
      typename Segment_traits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
        geom_traits->compare_endpoints_xy_2_object();

      if (cmp_seg_endpts(cv[0]) == SMALLER) {
        typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
          geom_traits->construct_min_vertex_2_object();
        cv.push_front(Subcurve_2(p, get_min_v(cv[0])));
      }
      else {
        typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
          geom_traits->construct_max_vertex_2_object();
        cv.push_front(Subcurve_2(p, get_max_v(cv[0])));
      }
    }

    /*! Append a point `p` to an existing polyline `xcv` at the front. */
    void operator()(const X_monotone_curve_2& xcv, Point_2& p) const
    {
      const Segment_traits_2* geom_traits =
        this->m_poly_traits.subcurve_traits_2();
      CGAL_precondition_code
        (
         typedef typename X_monotone_curve_2::size_type size_type;
         size_type num_seg = xcv.number_of_subcurves();
         typename Segment_traits_2::Compare_x_2 comp_x =
         geom_traits->compare_x_2_object();
         typename Segment_traits_2::Compare_xy_2 comp_xy =
         geom_traits->compare_xy_2_object();
         typename Base::Is_vertical_2 is_vertical =
           this->m_poly_traits.is_vertical_2_object();
         );
      CGAL_precondition(num_seg > 0);

      if (geom_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) {
        // xcv is oriented left-to-right
        typename Segment_traits_2::Construct_max_vertex_2 get_max_v =
          geom_traits->construct_max_vertex_2_object();
        CGAL_precondition
          (
           (!is_vertical(xcv) &&
            (comp_x(get_max_v(xcv[0]), p) == LARGER)) ||
           (is_vertical(xcv) &&
            (comp_x(get_max_v(xcv[0]), p) == EQUAL) &&
            (comp_xy(get_max_v(xcv[0]), p) == LARGER))
           );
        xcv.push_front(X_monotone_subcurve_2(p, get_max_v(xcv[0])));
      }
      else {
        // xcv is oriented right-to-left
        typename Segment_traits_2::Construct_min_vertex_2 get_min_v =
          geom_traits->construct_min_vertex_2_object();
        CGAL_precondition
          (
           (!is_vertical(xcv) &&
            (comp_x(get_min_v(xcv[0]), p) == SMALLER)) ||
           (is_vertical(xcv) &&
            (comp_x(get_min_v(xcv[0]), p) == EQUAL) &&
            (comp_xy(get_min_v(xcv[0]), p) == SMALLER))
           );
        xcv.push_front(X_monotone_subcurve_2(p, get_min_v(xcv[0])));
      }
    }
  };

  /*! Obtain a Push_front_2 functor object. */
  Push_front_2 push_front_2_object() const { return Push_front_2(*this); }

  /*! Construct a general curve. */
  class Construct_curve_2 : public Base::Construct_curve_2 {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>         Polyline_traits_2;

  public:
    /*! Constructor.
     */
    Construct_curve_2(const Polyline_traits_2& traits) :
      Base::Construct_curve_2(traits)
    {}

    /* Obtain an polyline connecting two given endpoints.
     */
    Curve_2 operator()(const Point_2& p, const Point_2& q) const
    { return Curve_2(Subcurve_2(p, q)); }

    /* Obtain a polyline consists of one given segment.
     */
    Curve_2 operator()(const Subcurve_2& seg) const
    { return Base::Construct_curve_2::operator()(seg); }

    /* Construct a well-oriented polyline from a range of either
     * `Segment_traits::Point_2` or `Segment_traits::Subcurve_2`.
     */
    template <typename ForwardIterator>
    Curve_2 operator()(ForwardIterator begin, ForwardIterator end) const
    {
      typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
      typedef typename boost::is_same<VT, Point_2>::type Is_point;
      // Dispatch the range to the appropriate implementation.
      return constructor_impl(begin, end, Is_point());
    }

    /*! Construction implementation from a range of segments.
     *  Note that the segments in the range are NOT necessarily x-monotone,
     *  thus it is impossible to test (even in precondition) whether the input
     *  forms a continuous and well oriented polyline.
     *  \pre Range should contain at least one segment.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end,
                             boost::false_type) const
    { return Base::Construct_curve_2::operator()(begin, end); }

    /*! Construction of a polyline from a range of points.
     * \pre The range contains at least two points
     * \pre Consecutive points are disjoint.
     * \return Well-oriented polyline connecting the given points. The order
     *         of the vertices is determined by their order in the range.
     *         Furthermore, the orientation of the polyline is induced by
     *         their order.
     */
    template <typename ForwardIterator>
    Curve_2 constructor_impl(ForwardIterator begin, ForwardIterator end,
                             boost::true_type) const
    {
      // The range must contain at least two points.
      CGAL_precondition_msg(std::distance(begin, end) > 1,
                            "Range of points must contain at least 2 points");

      // Container of the segments to be constructed from the range of points
      std::list<Subcurve_2> segs;

      CGAL_precondition_code
        (
         typename Segment_traits_2::Equal_2 equal =
         this->m_poly_traits.subcurve_traits_2()->equal_2_object();
         );

      // Check whether there are no points in the range:
      ForwardIterator next = begin;
      ForwardIterator curr = next++;

      // Construct a segment from each two adjacent points.
      Curve_2 cv;
      while (next != end) {
        CGAL_precondition_msg(!equal(*curr,*next),
                              "Cannot construct a degenerated segment");
        segs.push_back(Subcurve_2(*curr, *next));
        curr = next++;
      }

      return operator()(segs.begin(), segs.end());
    }
  };

  /*! Obtain a Construct_curve_2 functor object. */
  Construct_curve_2 construct_curve_2_object() const
  { return Construct_curve_2(*this); }

  /*! Construct an x-monotone curve. */
  class Construct_x_monotone_curve_2 :
    public Base::Construct_x_monotone_curve_2
  {
  protected:
    typedef Arr_polyline_traits_2<Segment_traits_2>     Polyline_traits_2;

  public:
    /*! Constructor.
     */
    Construct_x_monotone_curve_2(const Polyline_traits_2& traits) :
      Base::Construct_x_monotone_curve_2(traits)
    {}

    /*! Obtain an x-monotone polyline connecting two given endpoints.
     * \param p The first point.
     * \param q The second point.
     * \pre p and q must not be the same.
     * \return A segment connecting p and q.
     */
    X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q) const
    {
      CGAL_precondition_msg
        (!this->m_poly_traits.subcurve_traits_2()->equal_2_object()(p,q),
         "Cannot construct a degenerated segment as a polyline");
      X_monotone_subcurve_2 seg =  this->m_poly_traits.subcurve_traits_2()->
        construct_x_monotone_curve_2_object()(p, q);

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      if (this->m_poly_traits.subcurve_traits_2()->compare_xy_2_object()(p,q) ==
          LARGER)
        seg = this->m_poly_traits.subcurve_traits_2()->
          construct_opposite_2_object()(seg);
#endif

      return X_monotone_curve_2(seg);
    }

    /*! Obtain an x-monotone polyline that consists of one given segment.
     * \param seg input segment.
     * \pre seg is not degenerated.
     * \return An x-monotone polyline with one segment.
     */
    X_monotone_curve_2 operator()(const X_monotone_subcurve_2& seg) const
    { return Base::Construct_x_monotone_curve_2::operator()(seg); }

    /*! Construct an x-monotone polyline from a range of elements.
     * \pre Range should from a continuous well-oriented x-monotone polyline.
     */
    template <typename ForwardIterator>
    X_monotone_curve_2 operator()(ForwardIterator begin, ForwardIterator end)
      const
    {
      typedef typename std::iterator_traits<ForwardIterator>::value_type VT;
      typedef typename boost::is_same<VT, Point_2>::type Is_point;
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
    X_monotone_curve_2 constructor_impl(ForwardIterator begin,
                                        ForwardIterator end,
                                        boost::false_type) const
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
    X_monotone_curve_2 constructor_impl(ForwardIterator begin,
                                        ForwardIterator end,
                                        boost::true_type) const
    {
      // The range must contain at least two points.
      CGAL_precondition_msg(std::distance(begin, end) > 1,
                            "Range of points must contain at least 2 points");

      // Container of the segments to be constructed from the range of points
      std::list<X_monotone_subcurve_2> segs;
      // Make sure the range of points contains at least two points.
      ForwardIterator next = begin;
      ForwardIterator curr = next++;

      CGAL_precondition_code
        (
         const Segment_traits_2* geom_traits =
           this->m_poly_traits.subcurve_traits_2();
         // Initialize two comparison functors
         typename Segment_traits_2::Compare_x_2 compare_x =
           geom_traits->compare_x_2_object();
         typename Segment_traits_2::Compare_xy_2 compare_xy =
           geom_traits->compare_xy_2_object();
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

        segs.push_back(X_monotone_subcurve_2(*curr, *next));
        curr = next++;
      }

#ifdef CGAL_ALWAYS_LEFT_TO_RIGHT
      if (this->m_poly_traits.subcurve_traits_2()->
          compare_endpoints_xy_2_object()(*segs.begin()) == LARGER)
      {
        X_monotone_curve_2 xcv(segs.begin(), segs.end());
        return this->m_poly_traits.construct_opposite_2_object()(xcv);
      }
#endif

      return operator()(segs.begin(), segs.end());
    }
  };

  /*! Obtain a Construct_x_monotone_curve_2 functor object. */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const
  { return Construct_x_monotone_curve_2(*this); }

  /*! Deprecated!
   * Obtain the segment traits.
   * \return the segment traits.
   */
  CGAL_DEPRECATED const Segment_traits_2* segment_traits_2() const
  { return this->subcurve_traits_2(); }
};

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
