// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel    <efif@post.tau.ac.il>
//            Ron Wein     <wein@post.tau.ac.il>
//            (base on old version by: Iddo Hanniel)

#ifndef CGAL_ARR_NON_CACHING_SEGMENT_TRAITS_H
#define CGAL_ARR_NON_CACHING_SEGMENT_TRAITS_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file The non-caching segment traits-class for the arrangement package.
 * This traits class handles general segments. It is a model of the
 * ArrangementTraits_2 concept, a refinement of the ArrangementBasicTraits_2
 * concept. The class is templated by a kernel and inherits from the
 * Arr_non_caching_segment_basic_traits_2 class instanciated with the kernel -
 * a model of the ArrangementBasicTraits_2 concept. It defined a few additional
 * functors required by the concept it models.
 */

#include <boost/variant.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_non_caching_segment_basic_traits_2.h>
#include <CGAL/intersections.h>

namespace CGAL {

/*! \class
 * A model of the ArrangementTraits_2 concept that handles general
 * line segments.
 */
template <typename Kernel_T = Exact_predicates_exact_constructions_kernel>
class Arr_non_caching_segment_traits_2 :
  public Arr_non_caching_segment_basic_traits_2<Kernel_T>
{
public:
  typedef Kernel_T                                         Kernel;

  typedef Arr_non_caching_segment_basic_traits_2<Kernel>   Base;
  typedef typename Base::Segment_assertions                Segment_assertions;
  typedef typename Base::Has_exact_division                Has_exact_division;

  /*! Default constructor */
  Arr_non_caching_segment_traits_2() : Base() {}

  /// \name Types and functors inherited from the base
  //@{

  // Traits types:
  typedef typename Base::Has_left_category           Has_left_category;
  typedef typename Base::Has_do_intersect_category   Has_do_intersect_category;

  typedef typename Base::Left_side_category      Left_side_category;
  typedef typename Base::Bottom_side_category    Bottom_side_category;
  typedef typename Base::Top_side_category       Top_side_category;
  typedef typename Base::Right_side_category     Right_side_category;

  typedef typename Base::Point_2                     Point_2;
  typedef typename Base::X_monotone_curve_2          X_monotone_curve_2;
  typedef typename Base::Multiplicity                Multiplicity;

  /*! Compare the x-coordinates of two points */
  typedef typename Base::Compare_x_2            Compare_x_2;

  /*! Compare two points lexigoraphically; by x, then by y */
  typedef typename Base::Compare_xy_2           Compare_xy_2;

  /*! Obtain the left endpoint of a given segment */
  typedef typename Base::Construct_min_vertex_2 Construct_min_vertex_2;

  /*! Obtain the right endpoint of a given segment */
  typedef typename Base::Construct_max_vertex_2 Construct_max_vertex_2;

  /*! Check whether a given segment is vertical */
  typedef typename Base::Is_vertical_2          Is_vertical_2;

  /*! Return the location of a given point with respect to an input segment */
  typedef typename Base::Compare_y_at_x_2       Compare_y_at_x_2;

  /*! Check if two segments or if two points are identical */
  typedef typename Base::Equal_2                Equal_2;

  /*! Compare the y value of two segments immediately to the left of their
   * intersection point
   */
  typedef typename Base::Compare_y_at_x_left_2  Compare_y_at_x_left_2;

  /*! Compare the y value of two segments immediately to the right of their
   * intersection point
   */
  typedef typename Base::Compare_y_at_x_right_2 Compare_y_at_x_right_2;

  //@}

  /// \name Types and functors introduced here (based on the kernel)
  //@{

  // Categories:
  typedef Tag_true                              Has_merge_category;

  // Traits types:
  typedef X_monotone_curve_2                    Curve_2;

  /*! \class Make_x_monotone_2
   * A functor for subdividing a curve into x-monotone curves.
   */
  class Make_x_monotone_2 {
  public:
    /*! Subdivide a given curve into x-monotone subcurves and insert them into
     * a given output iterator. As segments are always x_monotone, only one
     * x-monotone curve is inserted into the output iterator.
     * \param cv the segment.
     * \param oi the output iterator for the result. Its dereference type is a
     *           variant that wraps a \c Point_2 or an \c X_monotone_curve_2
     *           objects.
     * \return the past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const
    {
      // Wrap the segment with a variant.
      typedef boost::variant<Point_2, X_monotone_curve_2>
        Make_x_monotone_result;
      *oi++ = Make_x_monotone_result(cv);
      return oi;
    }
  };

  /*! Obtain a Make_x_monotone_2 functor object. */
  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(); }

  /*! \class
   * A functor for splitting a segment into two segements.
   */
  class Split_2 {
    typedef Arr_non_caching_segment_traits_2<Kernel_T>    Self;

  public:
    /*! Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint).
     * \param c2 Output: The right resulting subcurve (p is its left endpoint).
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& cv, const Point_2& p,
                    X_monotone_curve_2& c1, X_monotone_curve_2& c2) const
    {
      Base base;

      // Make sure that p lies on the interior of the curve.
      CGAL_precondition_code(auto compare_xy = base.compare_xy_2_object());

      Construct_min_vertex_2 min_vertex = base.construct_min_vertex_2_object();
      Construct_max_vertex_2 max_vertex = base.construct_max_vertex_2_object();

      const Point_2& left = min_vertex(cv);
      const Point_2& right = max_vertex(cv);
      CGAL_precondition
        (Segment_assertions::_assert_is_point_on(p, cv, Has_exact_division()) &&
         compare_xy(left, p) == SMALLER &&
         compare_xy(right, p) == LARGER);

      typename Base::Construct_segment_2 construct_segment =
        base.construct_segment_2_object();

      Self self;
      if (self.compare_endpoints_xy_2_object()(cv) == SMALLER) {
        c1 = construct_segment(left, p);
        c2 = construct_segment(p, right);
      }
      else {
        c1 = construct_segment(p, left);
        c2 = construct_segment(right, p);
      }
    }
  };

  /*! Obtain a Split_2 functor object. */
  Split_2 split_2_object() const { return Split_2(); }

  /*! \class
   * A functor for computing intersections.
   */
  class Intersect_2 {
  protected:
    typedef Arr_non_caching_segment_traits_2<Kernel>        Traits;

    /*! The traits (in case it has state) */
    const Traits& m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Intersect_2(const Traits& traits) : m_traits(traits) {}

    friend class Arr_non_caching_segment_traits_2<Kernel>;

  public:
    /*! Find the intersections of the two given segments and insert them into
     * the given output iterator. As two segments may itersect only once, only
     * a single intersection will be contained in the iterator.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              OutputIterator oi) const
    {
      typedef std::pair<Point_2, Multiplicity>          Intersection_point;
      typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;

      const Kernel& kernel = m_traits;
      auto res = kernel.intersect_2_object()(cv1, cv2);

      // There is no intersection:
      if (! res) return oi;

      // Chack if the intersection is a point:
      const Point_2* p_p = boost::get<Point_2>(&*res);
      if (p_p != nullptr) {
        // Create a pair representing the point with its multiplicity,
        // which is always 1 for line segments for all practical purposes.
        // If the two segments intersect at their endpoints, then the
        // multiplicity is undefined, but we deliberately ignore it for
        // efficieny reasons.
        *oi++ = Intersection_result(Intersection_point(*p_p, 1));
        return oi;
      }

      // The intersection is a segment.
      const X_monotone_curve_2* cv_p = boost::get<X_monotone_curve_2>(&*res);
      CGAL_assertion(cv_p != nullptr);

      Comparison_result cmp1 = m_traits.compare_endpoints_xy_2_object()(cv1);
      Comparison_result cmp2 = m_traits.compare_endpoints_xy_2_object()(cv2);

      if (cmp1 == cmp2) {
        // cv1 and cv2 have the same directions, maintain this direction
        // in the overlap segment
        if (m_traits.compare_endpoints_xy_2_object()(*cv_p) != cmp1) {
          auto ctr_opposite = kernel.construct_opposite_segment_2_object();
          *oi++ = Intersection_result(ctr_opposite(*cv_p));
          return oi;
        }
      }
      *oi++ = Intersection_result(*cv_p);
      return oi;
    }
  };

  /*! Obtain an Intersect_2 functor object. */
  Intersect_2 intersect_2_object() const { return Intersect_2(*this); }

  /*! \class
   * A functor for testing whether two segments are mergeable.
   */
  class Are_mergeable_2 {
  protected:
    typedef Arr_non_caching_segment_traits_2<Kernel>    Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Are_mergeable_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_non_caching_segment_traits_2<Kernel>;

  public:
    /*! Check whether it is possible to merge two given x-monotone curves.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are mergeable, that is, if they are
     *         supported by the same line; (false) otherwise.
     * \pre cv1 and cv2 share a common endpoint.
     */
    bool operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2) const
    {
      const Base* base = m_traits;
      Equal_2 equal = base->equal_2_object();
      Construct_min_vertex_2 min_vertex = base->construct_min_vertex_2_object();
      Construct_max_vertex_2 max_vertex = base->construct_max_vertex_2_object();
      if (! equal(max_vertex(cv1), min_vertex(cv2)) &&
          ! equal(max_vertex(cv2), min_vertex(cv1)))
        return false;

      // Check if the two curves have the same supporting line.
      return (base->compare_slope_2_object()(cv1, cv2) == EQUAL);
    }
  };

  /*! Obtain an Are_mergeable_2 functor object */
  Are_mergeable_2 are_mergeable_2_object() const
  { return Are_mergeable_2(this); }

  /*! \class Merge_2
   * A functor that merges two x-monotone arcs into one.
   */
  class Merge_2 {
  protected:
    typedef Arr_non_caching_segment_traits_2<Kernel>    Traits;

    /*! The traits (in case it has state) */
    const Traits* m_traits;

    /*! Constructor
     * \param traits the traits (in case it has state)
     */
    Merge_2(const Traits* traits) : m_traits(traits) {}

    friend class Arr_non_caching_segment_traits_2<Kernel>;

  public:
    /*! Merge two given segments into a single segment.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The merged curve.
     * \pre The two curves are mergeable.
     */
    void operator()(const X_monotone_curve_2& cv1,
                    const X_monotone_curve_2& cv2,
                    X_monotone_curve_2& c) const
    {
      CGAL_precondition(m_traits->are_mergeable_2_object()(cv2, cv1));

      const Base* base = m_traits;
      Equal_2 equal = base->equal_2_object();
      Construct_min_vertex_2 min_vertex = base->construct_min_vertex_2_object();
      Construct_max_vertex_2 max_vertex = base->construct_max_vertex_2_object();

      // Check which curve extends to the right of the other.
      const Point_2& left1 = min_vertex(cv1);
      const Point_2& right1 = max_vertex(cv1);
      const Point_2& left2 = min_vertex(cv2);
      const Point_2& right2 = max_vertex(cv2);

      if (equal(right1, left2)) {
        // cv2 extends cv1 to the right.
        c = base->construct_segment_2_object()(left1, right2);
        return;
      }
      // cv1 extends cv2 to the right.
      CGAL_precondition(equal(right2, left1));

      c = base->construct_segment_2_object()(left2, right1);
    }
  };

  /*! Obtain a Merge_2 functor object */
  Merge_2 merge_2_object() const { return Merge_2(this); }
  //@}

  //! \name Functor definitions for the Boolean set-operations.
  //@{
  typedef typename Kernel::Construct_opposite_segment_2  Construct_opposite_2;

  /*! Obtain a Construct_opposite_2 functor object */
  Construct_opposite_2 construct_opposite_2_object() const
  { return Construct_opposite_2(); }

  class Compare_endpoints_xy_2 {
  public:
    /*!
     * Compare the two endpoints of a given curve lexigoraphically.
     * \param cv The curve.
     * \return SMALLER if cv is directed from left to right and LARGER
     * otherwise.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv) const
    {
      typedef typename Kernel::Construct_vertex_2     Construct_vertex_2;

      Kernel k;
      Base b;
      Construct_vertex_2 ctr_v = k.construct_vertex_2_object();
      Compare_xy_2 cmp_xy = b.compare_xy_2_object();
      return(cmp_xy(ctr_v(cv,0), ctr_v(cv,1)));
    }
  };

  /*! Obtain a Compare_endpoints_xy_2 functor object */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const
  { return Compare_endpoints_xy_2(); }
  //@}
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
